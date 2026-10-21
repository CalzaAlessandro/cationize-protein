# Add a noncovalent ligand

Use one cationization checkpoint and a ligand that already has validated
force-field parameters. This example uses ff14SB. Protein and ligand must
share the intended binding pose, coordinate frame and units. Within each
component, coordinates must follow its topology's atom order. A ligand PDB
provides coordinates and atom names; it does not supply charges, atom types
or bonded parameters.

Work in a new directory and retain the original exports. These instructions
assemble a dry complex. Solvation, ions, box selection, equilibration and MD
need a separate protocol. Covalent ligands and metal coordination need
additional bonded models.

## Know which files to use

| File | Role |
|---|---|
| `series/25/exports/protein.pdb` | Protein coordinates and actual modified-residue names; use it to inspect the ligand pose. |
| `system.parm7` + `system.rst7` | Fully parameterized Amber protein and matching coordinates. |
| `protein.lib` + `parameters.frcmod` + `rebuild.in` | System-specific Amber OFF unit, ordinary parameters and atom-type declarations for rebuilding this whole protein. |
| `residue_templates/ASC.lib`, `GLC.lib`, etc. | Individual fitted residue templates for building a protein with matching chemical and terminal states. |
| `forcefield.itp` | Exported GROMACS global defaults, protein atom types and any CMAP maps. |
| `protein.itp` | Exported GROMACS molecule definition, atoms and intramolecular interactions. |

The whole-protein `protein.lib` includes coordinates and connectivity across
residues, including disulfides. Loading it also creates the LEaP variable
`protein`. An individual residue library cannot replace this whole-protein
unit. The `.lib` extension here denotes Amber OFF format. A `.frcmod` supplies
force-field terms; a GROMACS `.itp` supplies topology sections. Renaming either
file does not convert it into the other format.

## Amber: combine the whole-protein unit with a ligand

Copy the checkpoint's `exports/` files to a new working directory. Supply
`ligand.mol2` with the ligand's validated charges, atom types and posed
coordinates, plus `ligand.frcmod` and its required base force field. The
example below assumes a ligand prepared with GAFF2. Use the base force field
that was used for your ligand.

First compare the protein types in `ATOM_TYPES.csv` with the ligand and its
base force field. Resolve conflicting definitions before combining them:
the protein bundle uses system-specific two-character types, and loading a
second definition can overwrite a type's parameters.

Copy `rebuild.in` to `complex.in`. Keep its complete `addAtomTypes` block,
`loadamberparams parameters.frcmod` and `loadoff protein.lib` lines. Put
`source leaprc.gaff2` before that block. Replace its final `saveamberparm`
and `quit` lines with:

```text
loadamberparams ligand.frcmod
ligand = loadmol2 ligand.mol2
complex = combine { protein ligand }
check complex
charge complex
savepdb complex complex.pdb
saveamberparm complex complex.parm7 complex.rst7
quit
```

Run `tleap -s -f complex.in` in that directory. Inspect `tleap.log` and resolve
missing parameters, type conflicts and unexpected close contacts. LEaP's
`combine` retains component coordinates and existing bonds; it adds no
protein–ligand bond.
This recipe uses the protein coordinates stored in `protein.lib`. If you
move the protein during ligand placement, update its coordinate inputs in
the matching atom order before assembly. Opening another PDB does not update
the OFF unit. The official [Amber ligand tutorial](https://ambermd.org/tutorials/basic/tutorial4b/index.php)
explains loading MOL2/OFF libraries, additional parameters and checking units;
the [Amber reference manuals](https://ambermd.org/Manuals.php) document LEaP's
`combine` command.

This LEaP recipe applies to an ff14SB export without CMAP. For ff19SB,
ordinary `.lib`/`.frcmod` reconstruction omits CMAP. Its `rebuild.in` therefore
directs you to the supplied `rebuild.py`; do not replace that guard with
the commands above. Use the existing `system.parm7`/`system.rst7` in the
Python route below to retain CMAP, and use an engine that supports it.

## GROMACS: convert two parameterized Amber components together

This route uses the installed ParmEd dependency and converts both components
in one pass. Supply `protein.parm7`/`protein.rst7` copied from the same export,
and an independently parameterized `ligand.parm7`/`ligand.rst7`. Both restarts
must already contain the intended complex pose. The ligand parameters must
use a compatible Amber nonbonded model and 1–4 treatment. This script does
not fit ligand parameters or place the ligand.

Save this as `make_complex.py` in a new directory with those four files,
activate the cationize-protein environment, and run `python make_complex.py`:

```python
import copy
import io
from pathlib import Path
import parmed as pmd

protein = pmd.load_file("protein.parm7", xyz="protein.rst7").copy(pmd.Structure)
ligand = pmd.load_file("ligand.parm7", xyz="ligand.rst7").copy(pmd.Structure)
for component, prefix in ((protein, "CP_"), (ligand, "LIG_")):
    for atom in component.atoms:
        atom.type = prefix + atom.type
        atom.atom_type = copy.copy(atom.atom_type)
        atom.atom_type.name = atom.type
        if atom.atom_type._bond_type is not None:
            atom.atom_type._bond_type = prefix + atom.atom_type._bond_type

complex_system = protein + ligand
complex_system.box = None
topology = pmd.gromacs.GromacsTopologyFile.from_structure(complex_system)
main, molecules = io.StringIO(), io.StringIO()
molecules.name = "complex_molecules.itp"
topology.write(main, parameters="inline", molfile=molecules, combine="all")
globals_text, separator, tail = main.getvalue().partition(
    '#include "complex_molecules.itp"'
)
assert separator and "[ defaults ]" in globals_text
Path("complex_forcefield.itp").write_text(
    globals_text[globals_text.index("[ defaults ]"):]
)
Path("complex_molecules.itp").write_text(molecules.getvalue())
Path("complex.top").write_text(
    '#include "complex_forcefield.itp"\n#include "complex_molecules.itp"\n' + tail
)
complex_system.save("complex.gro", format="gro", combine="all", precision=6)
```

The resulting `complex.top` includes globals first and molecular interactions
second. Bonded coefficients remain inline to preserve distinct terms that
can share atom-type labels. The script retains atom order, prefixes protein
types with `CP_` and ligand types with `LIG_`, and keeps explicit interactions
and any CMAP maps. `combine="all"` puts the whole complex into one GROMACS
`[ moleculetype ]` without creating new bonds. Keep these four output files
together. The GRO writer adds an enclosing box; choose the physical box and
solvation protocol separately. See [ParmEd's topology writer documentation](https://parmed.github.io/ParmEd/html/api/parmed/parmed.gromacs.gromacstop.html#parmed.gromacs.gromacstop.GromacsTopologyFile.write)
for its output options.

### If you already have a matching ligand ITP

Use the original exported `protein.itp` directly only when the ligand's
parameters already match the exported global defaults. Split any ligand
global definitions into `ligand_types.itp`; keep its molecular definition
in `ligand.itp`. The include order is:

```text
#include "forcefield.itp"
#include "ligand_types.itp"
#include "protein.itp"
#include "ligand.itp"

[ system ]
Cationized protein with ligand

[ molecules ]
system  1
LIG     1
```

Use the actual names from the two `[ moleculetype ]` sections in
`[ molecules ]`. Build complex coordinates in that same molecule order,
with each component's internal atom order unchanged. The exported protein
ITP groups all protein chains in one definition and retains their bonds.
Avoid the reserved `CP_` prefix for ligand types.

Keep exactly one `[ defaults ]` block. All global atom types and other
parameter types must precede **every** `[ moleculetype ]`. Do not include a
second stock `*.ff/forcefield.itp`. The protein export can use explicit 1–4
Lennard-Jones pairs with `gen-pairs = no` and `fudgeLJ = 1`; a ligand expecting
automatically generated pairs and `fudgeLJ = 0.5` needs conversion. Deleting
its `[ defaults ]` block alone does not make it compatible. The [GROMACS
topology specification](https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html)
defines the molecule and coordinate ordering; its [include-order guidance](https://manual.gromacs.org/current/user-guide/run-time-errors.html#invalid-order-for-directive-xxx)
explains global directive placement.

## Check the assembled complex

Check atom count and order, total charge, ligand identity and pose, modified
residues, and the protein bond graph against the original components.
Retain every original disulfide; the supplied HSA series has 17. Check CMAP
counts for ff19SB. For GROMACS, preprocess with your chosen, reviewed MDP:

```sh
gmx grompp -f your_protocol.mdp -c complex.gro -p complex.top \
  -o complex.tpr -pp processed.top
```

Resolve preprocessing warnings and inspect `processed.top`. Compare the
converted charges, Lennard-Jones parameters, exclusions and 1–4 interactions
with the source, then compare single-point energies at identical coordinates.
Successful file parsing alone does not validate a ligand force field or an
MD protocol.
