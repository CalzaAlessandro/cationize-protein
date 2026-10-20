"""Topology-derived protein identities and independently targetable carboxyls.

ParmEd is imported only by :func:`load_amber`; inspection also accepts objects
with ParmEd's atoms/residues/bonds/coordinates interface. Chemistry is declared
by residue state names and checked against covalent connectivity. Coordinates,
charges, PB radii, and LJ half-Rmin radii are read, never predicted or changed.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Mapping
from urllib.parse import quote

import numpy as np

from . import catalog


class UnsupportedTopology(ValueError):
    """Input chemistry or covalent identity cannot be interpreted unambiguously."""


LEGACY_ALIASES = {'ASL': 'ASC', 'GLS': 'GLC', 'NSL': 'ASX', 'LCL': 'LET'}
# Native Amber names -> names in the four preserved fitted templates.
_REUSED_MAPS = {
    'ASL': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'HA': 'H2',
            'CB': 'C3', 'CG': 'C4', 'OD1': 'O2'},
    'NSL': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'HA': 'H2',
            'CB': 'C3', 'CG': 'C4', 'OD1': 'O2'},
    'GLS': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'HA': 'H2',
            'CB': 'C3', 'CG': 'C4', 'CD': 'C5', 'OE1': 'O2'},
    'LCL': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'HA': 'H2',
            'CB': 'CB', 'CG': 'CG', 'CD1': 'CD1', 'CD2': 'CD2'},
}
_BACKBONE = frozenset(('N', 'CA', 'C', 'O', 'OXT'))


@dataclass(frozen=True)
class AtomInfo:
    index: int
    name: str
    residue_ordinal: int
    residue_name: str
    atomic_number: int


@dataclass(frozen=True)
class ResidueState:
    ordinal: int
    input_name: str
    canonical_name: str
    parent: str
    context: str
    modified_groups: tuple[str, ...]
    chain_id: str
    chain_ordinal: int
    original_id: str | None
    atom_indices: tuple[int, ...]
    native_atom_indices: dict[str, int]
    native_atom_names: dict[str, str]
    source_spec: dict | None


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    residue_ordinal: int
    group: str
    oxygen_indices: tuple[int, int]
    target_spec: dict


@dataclass(frozen=True)
class ProteinTopology:
    coordinates: np.ndarray
    charges_e: np.ndarray
    pb_radii_A: np.ndarray
    lj_radii_A: np.ndarray
    atoms: tuple[AtomInfo, ...]
    bonds: tuple[tuple[int, int], ...]
    peptide_bonds: tuple[tuple[int, int], ...]
    disulfide_bonds: tuple[tuple[int, int], ...]
    protected_atoms: tuple[int, ...]
    chains: tuple[tuple[int, ...], ...]
    residues: tuple[ResidueState, ...]
    candidates: tuple[Candidate, ...]
    provenance: dict


def _numeric_atoms(atoms, attribute, name, *, nonnegative=False):
    try:
        values = np.array([getattr(atom, attribute) for atom in atoms], dtype=float)
    except (AttributeError, TypeError, ValueError) as exc:
        raise UnsupportedTopology(f'Every atom must supply {name}') from exc
    if not np.isfinite(values).all() or (nonnegative and np.any(values < 0)):
        raise UnsupportedTopology(f'Invalid supplied {name}')
    return values


def _initial_states(residues, atom_indices):
    supported_parents = {spec['parent'] for spec in catalog.templates()} | {'ASH', 'GLH', 'LYN', 'TYM'}
    template_names = {spec['name'] for spec in catalog.templates()}
    rows = []
    for ordinal, residue in enumerate(residues, 1):
        name = str(residue.name).strip()
        canonical = LEGACY_ALIASES.get(name, name)
        spec = catalog.by_name(canonical) if canonical in template_names else None
        if spec is None and canonical not in supported_parents:
            raise UnsupportedTopology(f'Unsupported or nonprotein residue {name} at ordinal {ordinal}')
        mapping = _REUSED_MAPS.get(spec['reused_from'], {}) if spec else {}
        inverse = {actual: native for native, actual in mapping.items()}
        indices, names = {}, {}
        for atom in residue.atoms:
            actual = str(atom.name).strip()
            native = inverse.get(actual, actual)
            if not actual or native in indices:
                raise UnsupportedTopology(f'Ambiguous atom name in residue ordinal {ordinal}: {actual}')
            indices[native] = atom_indices[id(atom)]
            names[native] = actual
        if not _BACKBONE.difference({'OXT'}).issubset(indices):
            raise UnsupportedTopology(f'Incomplete protein backbone in residue ordinal {ordinal}')
        rows.append(dict(input_name=name, canonical_name=canonical,
                         parent=spec['parent'] if spec else {'ASH':'ASP','GLH':'GLU'}.get(canonical, canonical),
                         modified_groups=tuple(spec['modified_groups']) if spec else (),
                         native_atom_indices=indices, native_atom_names=names, source_spec=spec))
    return rows


def _bond_graph(structure, atoms, residues, atom_indices, residue_indices, rows):
    edges, graph = set(), [set() for _ in atoms]
    for bond in structure.bonds:
        try:
            i, j = atom_indices[id(bond.atom1)], atom_indices[id(bond.atom2)]
        except KeyError as exc:
            raise UnsupportedTopology('Bond references an atom outside the structure') from exc
        if i == j:
            raise UnsupportedTopology('Self-bond in topology')
        edge = tuple(sorted((i, j)))
        if edge in edges:
            raise UnsupportedTopology('Duplicate bond in topology')
        edges.add(edge)
        graph[i].add(j)
        graph[j].add(i)
    successors, predecessors = {}, {}
    peptide, disulfides = [], []
    for i, j in sorted(edges):
        ri = residue_indices[id(atoms[i].residue)]
        rj = residue_indices[id(atoms[j].residue)]
        if ri == rj:
            continue
        a, b = rows[ri]['native_atom_indices'], rows[rj]['native_atom_indices']
        if a.get('SG') == i and b.get('SG') == j:
            if rows[ri]['parent'] not in {'CYS', 'CYX'} or rows[rj]['parent'] not in {'CYS', 'CYX'}:
                raise UnsupportedTopology('Disulfide does not connect supported cysteine states')
            disulfides.append((i, j))
            continue
        if a['C'] == i and b['N'] == j:
            source, target = ri, rj
        elif a['N'] == i and b['C'] == j:
            source, target = rj, ri
        else:
            raise UnsupportedTopology('Unsupported interresidue covalent linkage')
        if source in successors or target in predecessors:
            raise UnsupportedTopology('Branched or ambiguous peptide connectivity')
        successors[source] = target
        predecessors[target] = source
        peptide.append((rows[source]['native_atom_indices']['C'], rows[target]['native_atom_indices']['N']))
    for index, residue in enumerate(residues):
        members = {atom_indices[id(atom)] for atom in residue.atoms}
        seed = next(iter(members))
        connected, stack = {seed}, [seed]
        while stack:
            node = stack.pop()
            for neighbor in graph[node] & members - connected:
                connected.add(neighbor)
                stack.append(neighbor)
        if connected != members:
            raise UnsupportedTopology(f'Disconnected atoms within residue ordinal {index + 1}')
        native = rows[index]['native_atom_indices']
        for first, second in [('N', 'CA'), ('CA', 'C'), ('C', 'O')]:
            if native[second] not in graph[native[first]]:
                raise UnsupportedTopology(f'Missing canonical backbone bond in residue ordinal {index + 1}')
    return tuple(sorted(edges)), graph, successors, predecessors, tuple(sorted(peptide)), tuple(disulfides)


def _peptide_chains(residues, successors, predecessors):
    chains, visited = [], set()
    for seed in range(len(residues)):
        if seed in predecessors:
            continue
        chain, current = [], seed
        while current is not None:
            if current in visited:
                raise UnsupportedTopology('Cyclic or ambiguous peptide topology')
            visited.add(current)
            chain.append(current)
            current = successors.get(current)
        if len(chain) < 2:
            raise UnsupportedTopology('Single-residue peptide chains are unsupported')
        explicit = {str(getattr(residues[index], 'chain', '') or '').strip() for index in chain}
        explicit.discard('')
        if len(explicit) > 1:
            raise UnsupportedTopology('Chain labels contradict a continuous peptide chain')
        label = next(iter(explicit)) if explicit else f'chain{len(chains) + 1}'
        chains.append((tuple(chain), label))
    if len(visited) != len(residues):
        raise UnsupportedTopology('Cyclic peptide chains are unsupported')
    return chains


def _disulfide_states(rows, graph, disulfides, atoms):
    endpoints = {index for pair in disulfides for index in pair}
    if len(endpoints) != 2 * len(disulfides):
        raise UnsupportedTopology('A sulfur participates in multiple disulfides')
    for row in rows:
        sg = row['native_atom_indices'].get('SG')
        if sg in endpoints and int(atoms[sg].atomic_number) != 16:
            raise UnsupportedTopology('An SG-SG bond endpoint is not sulfur')
        if sg in endpoints:
            if any(int(getattr(atoms[neighbor], 'atomic_number', 0) or 0) == 1 for neighbor in graph[sg]):
                raise UnsupportedTopology('Disulfide sulfur still has a bonded hydrogen')
            if row['source_spec'] is not None and row['parent'] != 'CYX':
                raise UnsupportedTopology('Modified cysteine template contradicts actual disulfide')
            row['parent'] = 'CYX'
            if row['source_spec'] is None:
                row['canonical_name'] = 'CYX'
        elif row['parent'] == 'CYX':
            raise UnsupportedTopology('CYX residue has no actual SG-SG bond')
    return endpoints


def _oxygen_pair(native, graph, carbon, oxygen_names, ordinal, atoms, *, protonated=False):
    if carbon not in native or any(name not in native for name in oxygen_names):
        raise UnsupportedTopology(f'Incomplete carboxyl group in residue ordinal {ordinal}')
    center = native[carbon]
    oxygens = tuple(native[name] for name in oxygen_names)
    if any(oxygen not in graph[center] for oxygen in oxygens):
        raise UnsupportedTopology(f'Carboxyl oxygen is not bonded to its carbon in residue ordinal {ordinal}')
    if int(atoms[center].atomic_number) != 6 or any(int(atoms[oxygen].atomic_number) != 8 for oxygen in oxygens):
        raise UnsupportedTopology(f'Carboxyl atom elements are inconsistent at residue ordinal {ordinal}')
    extras = [graph[oxygen] - {center} for oxygen in oxygens]
    if protonated:
        if extras[0] or len(extras[1]) != 1 or int(atoms[next(iter(extras[1]))].atomic_number) != 1:
            raise UnsupportedTopology(f'Protonated acid requires exactly one H on {oxygen_names[1]} at residue ordinal {ordinal}')
        hydrogen = next(iter(extras[1]))
        if graph[hydrogen] != {oxygens[1]}:
            raise UnsupportedTopology('An acid hydrogen has inconsistent connectivity')
    elif any(extras):
        raise UnsupportedTopology(f'Carboxyl oxygen has an unsupported substituent or proton at residue ordinal {ordinal}')
    return oxygens


def _original_id(residue):
    number = getattr(residue, 'number', None)
    if number is None:
        return None
    insertion = str(getattr(residue, 'insertion_code', '') or '').strip()
    return f'{number}{insertion}'


def inspect_structure(structure, *, provenance: Mapping | None = None) -> ProteinTopology:
    """Validate a protein-only topology and enumerate remaining COO⁻ groups.

    Peptide chains follow actual interresidue C-N bonds, independently of SG-SG
    links and residue numbering. Each chain must have at least two residues.
    Modified labels must agree with their actual terminal context. ASH/GLH retain
    their explicit acid H and chemical target; ambiguous HIS, caps, solvent, ligands, branched/cyclic peptides, and unknown
    covalent crosslinks stop explicitly. No geometric bond inference is used.

    ``native_atom_indices`` maps canonical native names to global atom indices;
    ``native_atom_names`` maps those names to the input labels. Unmapped added
    atoms retain their input labels. Candidate target specs include all existing
    modified groups, making side and terminal amidation independent operations.
    """
    atoms, residues = list(structure.atoms), list(structure.residues)
    if not atoms or not residues:
        raise UnsupportedTopology('Empty protein topology')
    atom_indices = {id(atom): index for index, atom in enumerate(atoms)}
    residue_indices = {id(residue): index for index, residue in enumerate(residues)}
    if len(atom_indices) != len(atoms) or len(residue_indices) != len(residues):
        raise UnsupportedTopology('Repeated atom or residue objects')
    seen = set()
    for residue in residues:
        if not residue.atoms:
            raise UnsupportedTopology('Empty residue in topology')
        for atom in residue.atoms:
            if id(atom) not in atom_indices or atom.residue is not residue or id(atom) in seen:
                raise UnsupportedTopology('Atom/residue membership is inconsistent')
            seen.add(id(atom))
    if len(seen) != len(atoms):
        raise UnsupportedTopology('Atoms are missing from residue membership')
    coordinates = np.array(structure.coordinates, dtype=float, copy=True)
    if coordinates.shape != (len(atoms), 3) or not np.isfinite(coordinates).all():
        raise UnsupportedTopology('Finite coordinates are required for every topology atom')
    atomic_numbers = _numeric_atoms(atoms, 'atomic_number', 'atomic numbers')
    if np.any(atomic_numbers < 1) or np.any(atomic_numbers > 118) or np.any(atomic_numbers != np.floor(atomic_numbers)):
        raise UnsupportedTopology('Atoms must have explicit physical element identities; virtual sites are unsupported')
    charges = _numeric_atoms(atoms, 'charge', 'topology charges')
    pb_radii = _numeric_atoms(atoms, 'solvent_radius', 'topology PB radii', nonnegative=True)
    if np.any(pb_radii <= 0):
        raise UnsupportedTopology('Protein PB radii must be positive; unparameterized inputs are unsupported')
    lj_radii = _numeric_atoms(atoms, 'rmin', 'topology LJ half-Rmin radii', nonnegative=True)
    rows = _initial_states(residues, atom_indices)
    bonds, graph, successors, predecessors, peptide, disulfides = _bond_graph(
        structure, atoms, residues, atom_indices, residue_indices, rows)
    chains = _peptide_chains(residues, successors, predecessors)
    protected = _disulfide_states(rows, graph, disulfides, atoms)
    chain_lookup = {}
    for chain_number, (chain, label) in enumerate(chains, 1):
        for index in chain:
            context = 'n_terminal' if index == chain[0] else 'c_terminal' if index == chain[-1] else 'internal'
            chain_lookup[index] = (chain_number, label, context)
    states, candidates = [], []
    for index, (residue, row) in enumerate(zip(residues, rows)):
        ordinal = index + 1
        chain_number, chain_id, context = chain_lookup[index]
        native = row['native_atom_indices']
        groups = set(row['modified_groups'])
        if row['source_spec'] and row['source_spec']['context'] != context:
            raise UnsupportedTopology(f'Modified state context contradicts peptide bonds at residue ordinal {ordinal}')
        if 'OXT' in native and (context != 'c_terminal' or 'terminal' in groups):
            raise UnsupportedTopology(f'Unexpected unmodified terminal oxygen at residue ordinal {ordinal}')
        backbone = {native[name] for name in _BACKBONE if name in native}
        protected.update(backbone)
        protected.update(neighbor for atom in backbone for neighbor in graph[atom]
                         if int(getattr(atoms[neighbor], 'atomic_number', 0) or 0) == 1)
        original = _original_id(residue)
        state = ResidueState(
            ordinal=ordinal, input_name=row['input_name'], canonical_name=row['canonical_name'],
            parent=row['parent'], context=context, modified_groups=tuple(sorted(groups)),
            chain_id=chain_id, chain_ordinal=chain_number, original_id=original,
            atom_indices=tuple(atom_indices[id(atom)] for atom in residue.atoms),
            native_atom_indices=dict(native), native_atom_names=dict(row['native_atom_names']),
            source_spec=row['source_spec'],
        )
        states.append(state)
        remaining = []
        if state.parent in {'ASP', 'GLU'}:
            carbon, oxygens = ('CG', ('OD1', 'OD2')) if state.parent == 'ASP' else ('CD', ('OE1', 'OE2'))
            if 'side' not in groups:
                remaining.append(('side', _oxygen_pair(native, graph, carbon, oxygens, ordinal, atoms,
                                                       protonated=state.canonical_name in {'ASH','GLH'})))
            elif oxygens[1] in native:
                raise UnsupportedTopology(f'Modified side-chain label retains an unmodified COO group at residue ordinal {ordinal}')
        if context == 'c_terminal' and 'terminal' not in groups:
            remaining.append(('terminal', _oxygen_pair(native, graph, 'C', ('O', 'OXT'), ordinal, atoms)))
        for group, oxygen_indices in remaining:
            target = catalog.resolve(state.parent, context, groups | {group})
            identity = f'r{ordinal}:chain={quote(chain_id, safe="")}:pdb={quote(original or "unknown", safe="")}:group={group}'
            candidates.append(Candidate(identity, ordinal, group, oxygen_indices, target))
    atom_info = tuple(AtomInfo(index, str(atom.name), residue_indices[id(atom.residue)] + 1,
                               str(atom.residue.name), int(getattr(atom, 'atomic_number', 0) or 0))
                      for index, atom in enumerate(atoms))
    return ProteinTopology(
        coordinates, charges, pb_radii, lj_radii, atom_info, bonds, peptide, disulfides,
        tuple(sorted(protected)), tuple(tuple(index + 1 for index in chain) for chain, _ in chains),
        tuple(states), tuple(candidates), dict(provenance or {}),
    )


def _sha_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def load_amber(topology_path, coordinates_path, *, identity_pdb_path=None) -> ProteinTopology:
    """Load a protein-only Amber topology/restart, optionally restoring PDB IDs.

    Optional PDB identities must match the topology's atom/residue order and
    exact atom labels; coordinates always come from the supplied restart. Input
    files remain untouched and SHA-256 provenance is returned with the result.
    """
    paths = {'topology': Path(topology_path), 'coordinates': Path(coordinates_path)}
    if identity_pdb_path is not None:
        paths['identity_pdb'] = Path(identity_pdb_path)
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    try:
        import parmed
    except ImportError as exc:
        raise ImportError('Loading Amber inputs requires the external ParmEd dependency') from exc
    try:
        structure = parmed.load_file(str(paths['topology']), xyz=str(paths['coordinates']))
        if identity_pdb_path is not None:
            identity = parmed.load_file(str(paths['identity_pdb']))
            if len(identity.residues) != len(structure.residues) or [a.name for a in identity.atoms] != [a.name for a in structure.atoms]:
                raise UnsupportedTopology('Identity PDB atom/residue order differs from topology')
            for source, target in zip(identity.residues, structure.residues):
                if LEGACY_ALIASES.get(source.name, source.name) != LEGACY_ALIASES.get(target.name, target.name) or len(source.atoms) != len(target.atoms):
                    raise UnsupportedTopology('Identity PDB residue labels differ from topology')
                target.chain = source.chain
                target.number = source.number
                target.insertion_code = source.insertion_code
    except UnsupportedTopology:
        raise
    except Exception as exc:
        raise UnsupportedTopology(f'Cannot load Amber topology and coordinates: {exc}') from exc
    provenance = {f'{name}_path': str(path.resolve()) for name, path in paths.items()}
    provenance.update({f'{name}_sha256': _sha_file(path) for name, path in paths.items()})
    return inspect_structure(structure, provenance=provenance)
