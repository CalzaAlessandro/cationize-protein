"""Portable GROMACS include files with isolated protein atom-type names."""
from __future__ import annotations

import copy
import io
from pathlib import Path


def write_gromacs(structure, directory):
    """Write global parameters, molecular definitions, system and coordinates.

    Bonded coefficients remain inline to avoid ambiguity with another molecule's
    bonded types. CMAP and all global nonbonded definitions stay in forcefield.itp.
    All protein chains share one molecule definition without adding any bonds.
    """
    import parmed as pmd

    directory=Path(directory)
    # ParmEd's generic deepcopy omits its symmetry attribute; its public copy
    # constructor also preserves a complete, independently typed Structure.
    isolated=structure.copy(pmd.Structure)
    for atom in isolated.atoms:
        atom.type='CP_'+atom.type
        atom.atom_type=copy.copy(atom.atom_type)
        atom.atom_type.name=atom.type
        if atom.atom_type._bond_type is not None:
            atom.atom_type._bond_type='CP_'+atom.atom_type._bond_type
    topology=pmd.gromacs.GromacsTopologyFile.from_structure(isolated)
    top=io.StringIO()
    molecule=io.StringIO()
    molecule.name='protein.itp'
    topology.write(top,parameters='inline',molfile=molecule,combine='all')
    include='#include "protein.itp"'
    before,separator,after=top.getvalue().partition(include)
    if not separator or '[ defaults ]' not in before:
        raise ValueError('GROMACS writer did not produce the expected molecular include layout')
    # Strip ParmEd's machine/user/command header; keep numeric data unchanged.
    globals_text=before[before.index('[ defaults ]'):]
    (directory/'forcefield.itp').write_text(
        '; cationize-protein: physical global definitions for this protein.\n'
        '; Include exactly once, before all molecule definitions.\n'+globals_text)
    (directory/'protein.itp').write_text(
        '; cationize-protein: protein molecule definitions, including modified residues.\n'
        '; Native atom order and connected molecular components are preserved.\n'+molecule.getvalue())
    (directory/'system.top').write_text(
        '; cationize-protein: standalone topology with relative, portable includes.\n'
        '#include "forcefield.itp"\n#include "protein.itp"\n'+after)
    isolated.save(str(directory/'system.gro'),format='gro',overwrite=False,combine=True)
    (directory/'GROMACS_README.md').write_text('''# GROMACS files

`system.top` includes `forcefield.itp` (global defaults, protein atom types and
CMAP maps) and `protein.itp` (molecule definitions, atoms, bonds and other
intramolecular interactions). `system.gro` and `protein.pdb` have the same atom
order. Keep the files together; all include paths are relative.

For an independently parameterized ligand, retain exactly one compatible
`[ defaults ]` block. Add ligand global definitions before any molecule
definitions, include the ligand molecule ITP, and update `[ molecules ]` and
the complex coordinates consistently. Protein atom types use the `CP_` prefix;
avoid that namespace for ligand types. The ligand must be compatible with the
selected Amber force-field family, combination rule and 1-4 treatment. Adding
an ITP does not parameterize a ligand or create protein-ligand covalent bonds.

The protein ITP groups all input chains in one molecule definition while
preserving the original atom order and covalent graph; no interchain bonds are
created. It does not treat residues as independent molecules. The original
fitted residue libraries and complete Amber parameter bundle remain separate.
No dynamics protocol, ligand placement or small-molecule parameterization is
performed by this export.
''')
    return ['system.top','system.gro','forcefield.itp','protein.itp','GROMACS_README.md']
