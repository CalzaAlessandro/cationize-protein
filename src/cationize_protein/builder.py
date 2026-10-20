"""Coordinate-preserving, unparameterized residue assembly for arbitrary proteins.

ParmEd is imported only when assembly is requested. No Amber installation or
executable is needed. Full force-field assignment belongs to the caller's
backend (for example OpenMM). Existing bonded coefficients are deliberately
not reused after a chemical mutation.
"""
from __future__ import annotations

from collections import deque
from collections.abc import Mapping
from copy import copy
from itertools import combinations
import numbers

import numpy as np


def _parmed():
    try:
        import parmed
    except ImportError as exc:
        raise ImportError('Residue assembly requires the optional ParmEd package') from exc
    return parmed


def _named_atoms(atoms, label):
    result = {atom.name: atom for atom in atoms}
    if len(result) != len(atoms) or any(not name for name in result):
        raise ValueError(f'{label} has empty or duplicate atom names')
    return result


def _xyz(atoms):
    try:
        result = np.array([[atom.xx, atom.xy, atom.xz] for atom in atoms], dtype=float)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError('Every atom needs finite three-dimensional coordinates') from exc
    if result.shape != (len(atoms), 3) or not np.isfinite(result).all():
        raise ValueError('Every atom needs finite three-dimensional coordinates')
    return result


def _frame(origin, first, second):
    x = first - origin
    length = np.linalg.norm(x)
    if length < 1e-10:
        return None
    x = x / length
    y = second - origin
    y -= np.dot(y, x) * x
    length = np.linalg.norm(y)
    if length < 1e-10:
        return None
    y /= length
    return np.column_stack((x, y, np.cross(x, y)))


def _placed_coordinates(template, inherited, elements):
    """Fit each added component rigidly about its own retained attachment atom."""
    atoms = {atom.name: atom for atom in template.atoms}
    original = dict(zip(atoms, _xyz(template.atoms)))
    output = {name: xyz.copy() for name, xyz in original.items()}
    retained = {name: _xyz([atom])[0] for name, atom in inherited.items()}
    graph = {name: set() for name in atoms}
    for bond in template.bonds:
        a, b = bond.atom1.name, bond.atom2.name
        if a not in graph or b not in graph or a == b:
            raise ValueError('Invalid template bond graph')
        graph[a].add(b); graph[b].add(a)
    pending = set(atoms) - set(retained)
    while pending:
        component = {min(pending)}
        queue = list(component)
        while queue:
            name = queue.pop()
            for neighbor in graph[name] & pending - component:
                component.add(neighbor); queue.append(neighbor)
        pending -= component
        attachments = {neighbor for name in component for neighbor in graph[name]
                       if neighbor in retained}
        if len(attachments) != 1:
            raise ValueError('Each added component must have exactly one retained attachment atom')
        root = attachments.pop()
        # Local oxygen and carbon anchors define a carbonyl frame. Searching
        # farther along retained bonds also supports new hydrogens on N atoms.
        distance, queue = {root: 0}, deque([root])
        while queue:
            name = queue.popleft()
            for neighbor in sorted(graph[name] & set(retained)):
                if neighbor not in distance:
                    distance[neighbor] = distance[name] + 1; queue.append(neighbor)
        candidates = sorted((name for name in distance if name != root),
                            key=lambda name: (distance[name],
                                              {8: 0, 6: 1}.get(elements[name], 2), name))
        frames = None
        for first, second in combinations(candidates, 2):
            source = _frame(original[root], original[first], original[second])
            target = _frame(retained[root], retained[first], retained[second])
            if source is not None and target is not None:
                frames = source, target
                break
        if frames is None:
            raise ValueError('Cannot define a noncollinear local placement frame')
        source, target = frames
        for name in component:
            output[name] = (original[name] - original[root]) @ source @ target.T + retained[root]
    output.update(retained)
    return output


def replace_residue(structure, residx: int, template,
                    native_to_template_map: Mapping[str, str] | None = None):
    """Return a new ParmEd Structure; ``residx`` is a zero-based residue ordinal.

    The input graph is authoritative: every bond outside the replaced residue,
    including all peptide links, disulfides and other covalent attachments, is
    retained. Removing an externally attached atom is rejected. Template atoms
    supply product charges, types and intraresidue bonds; no bonded force-field
    coefficients are assigned. Source and template objects remain unchanged.

    An explicit map identifies *all* inherited atoms whose coordinates must be
    retained, allowing arbitrary legacy names. Without a map, identical names
    are inherited. Unknown OFF elements on inherited atoms and zero OFF masses
    are resolved from native/elemental metadata without modifying the library.
    Other template atoms are added by rigid local-frame fitting
    per connected new component. Multiple attachment atoms in one new component
    are rejected. Subsequent steric optimization and physical relaxation are
    separate steps; this routine does not claim a minimized product geometry.
    """
    pmd = _parmed()
    if (isinstance(residx, bool) or not isinstance(residx, numbers.Integral)
            or not 0 <= residx < len(structure.residues)):
        raise ValueError('residx must be a valid zero-based residue ordinal')
    if not template.atoms or not template.name:
        raise ValueError('A named, nonempty product template is required')
    _xyz(structure.atoms)
    old_residue = structure.residues[int(residx)]
    old = _named_atoms(old_residue.atoms, 'Input residue')
    product = _named_atoms(template.atoms, 'Product template')
    mapping = ({name: name for name in old.keys() & product.keys()}
               if native_to_template_map is None else dict(native_to_template_map))
    if (not mapping or len(set(mapping.values())) != len(mapping) or
            any(a not in old or b not in product for a, b in mapping.items())):
        raise ValueError('Atom map must uniquely match existing native and template atoms')
    for native, new in mapping.items():
        if product[new].atomic_number > 0 and old[native].atomic_number != product[new].atomic_number:
            raise ValueError(f'Inherited atom element changes: {native} to {new}')
    for bond in structure.bonds:
        for atom, other in ((bond.atom1, bond.atom2), (bond.atom2, bond.atom1)):
            if atom.residue is old_residue and other.residue is not old_residue and atom.name not in mapping:
                raise ValueError(f'Deleted atom is an external covalent attachment: {atom.name}')
    inherited = {new: old[native] for native, new in mapping.items()}
    elements = {}
    for name, atom in product.items():
        z = atom.atomic_number
        if z <= 0 and name in inherited:
            z = inherited[name].atomic_number
        if not 1 <= z < len(pmd.periodic_table.Element):
            raise ValueError(f'Unknown element for product atom {name}; explicit metadata is required')
        elements[name] = z
    placed = _placed_coordinates(template, inherited, elements)
    result = pmd.Structure()
    old_to_new, product_to_new = {}, {}
    for ri, residue in enumerate(structure.residues):
        renamed = pmd.Residue(template.name if ri == residx else residue.name,
                              residue.number, chain=residue.chain,
                              insertion_code=residue.insertion_code, segid=residue.segid,
                              list=result.residues)
        renamed.ter = residue.ter
        result.residues.append(renamed)
        atoms = template.atoms if ri == residx else residue.atoms
        for atom in atoms:
            new = copy(atom)
            if ri == residx:
                new.atomic_number = elements[atom.name]
                native = inherited.get(atom.name)
                if native is not None and native.mass > 0:
                    new.mass = native.mass
                elif not np.isfinite(new.mass) or new.mass <= 0:
                    new.mass = pmd.periodic_table.Mass[pmd.periodic_table.Element[new.atomic_number]]
                new.xx, new.xy, new.xz = placed[atom.name]
                product_to_new[atom.name] = new
            else:
                old_to_new[atom] = new
            renamed.add_atom(new); result.atoms.append(new)
    old_to_new.update({old[native]: product_to_new[new] for native, new in mapping.items()})
    for bond in structure.bonds:
        if bond.atom1.residue is old_residue and bond.atom2.residue is old_residue:
            continue
        result.bonds.append(pmd.Bond(old_to_new[bond.atom1], old_to_new[bond.atom2],
                                     order=bond.order, qualitative_type=bond.qualitative_type))
    for bond in template.bonds:
        result.bonds.append(pmd.Bond(product_to_new[bond.atom1.name], product_to_new[bond.atom2.name],
                                     order=bond.order, qualitative_type=bond.qualitative_type))
    if structure.box is not None:
        result.box = np.array(structure.box, copy=True)
    return result


def export_mutation_pdb(structure) -> tuple[str, dict]:
    """Return PDB text and a backend-neutral, JSON-compatible graph manifest.

    Atom and residue indices in the manifest are zero based. PDB records use
    consecutive one-based indices and include every graph edge as CONECT, plus
    TER at declared or disconnected chain boundaries. Original residue metadata
    remains in the manifest. PDB coordinates have its usual 0.001 Å precision;
    use ``coordinates_A`` to retain the full precision in a simulation backend.
    No files are written and no runtime force-field program is invoked.
    """
    pmd = _parmed()
    xyz = _xyz(structure.atoms)
    if len(structure.atoms) > 99999 or len(structure.residues) > 9999:
        raise ValueError('Structure exceeds PDB index limits; use a direct topology or mmCIF backend')
    serials = {atom: i for i, atom in enumerate(structure.atoms)}
    residues = {residue: i for i, residue in enumerate(structure.residues)}
    bonds, crosslinks, residue_pairs = [], [], set()
    for bond in structure.bonds:
        i, j = serials[bond.atom1], serials[bond.atom2]
        r1, r2 = residues[bond.atom1.residue], residues[bond.atom2.residue]
        bonds.append([i, j, float(bond.order)])
        if r1 != r2:
            crosslinks.append(dict(residue1=r1, atom1=bond.atom1.name,
                                  residue2=r2, atom2=bond.atom2.name, order=float(bond.order)))
            residue_pairs.add(tuple(sorted((r1, r2))))
    lines, atoms_info, residues_info = [], [], []
    for ri, residue in enumerate(structure.residues):
        _named_atoms(residue.atoms, 'Export residue')
        if (len(residue.name) > 3 or len(residue.chain) > 1 or
                len(residue.insertion_code) > 1):
            raise ValueError('Residue name, chain or insertion code exceeds PDB field width')
        residues_info.append(dict(index=ri, name=residue.name, source_number=residue.number,
                                  chain=residue.chain, insertion_code=residue.insertion_code,
                                  segid=residue.segid, declared_ter=bool(residue.ter)))
        for atom in residue.atoms:
            if (not atom.name or len(atom.name) > 4 or atom.altloc or atom.other_locations
                    or not 1 <= atom.atomic_number < len(pmd.periodic_table.Element)):
                raise ValueError('PDB export needs resolved atom names, elements and alternate locations')
            index = serials[atom]
            x, y, z = xyz[index]
            if any(len(f'{value:8.3f}') > 8 for value in (x, y, z)):
                raise ValueError('Coordinates exceed PDB fixed-column limits')
            element = pmd.periodic_table.Element[atom.atomic_number].upper()
            lines.append(f'ATOM  {index+1:5d} {atom.name:>4s} {residue.name:>3s} '
                         f'{residue.chain:1s}{ri+1:4d}{residue.insertion_code:1s}   '
                         f'{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}  \n')
            atoms_info.append(dict(index=index, residue=ri, name=atom.name,
                                   atomic_number=atom.atomic_number, amber_type=str(atom.type),
                                   charge_e=float(atom.charge)))
        next_res = structure.residues[ri+1] if ri+1 < len(structure.residues) else None
        if (residue.ter or next_res is None or next_res.chain != residue.chain
                or (ri, ri+1) not in residue_pairs):
            lines.append('TER\n')
    for i, j, order in bonds:
        lines.append(f'CONECT{i+1:5d}{j+1:5d}\n')
    lines.append('END\n')
    return ''.join(lines), dict(schema_version=1, requires_parameterization=True,
                                coordinates_A=xyz.tolist(), atoms=atoms_info,
                                residues=residues_info, bonds=bonds,
                                interresidue_bonds=crosslinks,
                                box_A_degrees=None if structure.box is None else np.asarray(structure.box).tolist())
