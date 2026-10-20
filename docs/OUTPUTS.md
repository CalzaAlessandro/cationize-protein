# Structures and reusable parameters

## Directory layout

```text
cationized_series/
  0/                         native baseline
  25/                        one requested percentage (same layout below)
    protein.pdb
    system.xml
    state.json
    STAGE_COMPLETE.json
    protonation/             PROPKA results and decisions
    exports/
      protein.pdb
      system.parm7
      system.rst7
      protein.lib
      parameters.frcmod
      system.xml
      state.json
      system.top
      system.gro
      forcefield.itp
      protein.itp
      PARAMETERS.json
      ATOM_TYPES.csv
      EXPORTS.json
      rebuild.py
      rebuild.in
      PARAMETERS_README.md
      GROMACS_README.md
  50/
  75/
  100/
  internal/                  per-step minimization, fields and restart evidence
  RUN.json
  REFERENCE.json
  RUN_COMPLETE.json
```

Engine-specific files depend on `--formats`. Custom percentages replace the
default percentage directories. `internal/` is managed automatically; it is
not an MD production or analysis directory. JSON files record settings,
completion and hashes; normal operation does not require editing them.

## Protein-only use

Every percentage contains a PDB with hydrogens and the actual modified residue
names. The export PDB uses the same atom order as the numerical topology files.
PDB and ordinary GRO coordinates have finite decimal precision; `state.json`
and the numeric parameter bundle retain full precision.

Amber uses `system.parm7` and `system.rst7`; no original leaprc is needed to
read them. OpenMM uses `system.xml` with the supplied coordinates. GROMACS uses
`system.top` and `system.gro`. These are parameterized protein structures,
without a production-dynamics or solvation protocol.

## GROMACS includes

1. `forcefield.itp` contains one global defaults block, protein atom types and
   CMAP definitions. It is self-contained, not a reference to an installed
   force-field folder.
2. `protein.itp` contains the molecule definition, atoms, bonds, pairs, angles,
   dihedrals, exclusions and CMAP assignments. Bonded coefficients are inline.
3. `system.top` contains relative includes, system description and molecule
   count. Keep the ITP files alongside it when moving the directory.

All protein chains share one molecule definition to preserve atom order.
No bonds are added and existing disulfides remain intact. Residues are not
exported as independent molecules.

An independently parameterized ligand can provide its own global definitions
and molecule ITP. Keep exactly one compatible `[ defaults ]` block, put all
global type definitions before molecule definitions, update `[ molecules ]`,
and supply matching complex coordinates. Avoid the protein `CP_` atom-type
namespace. Combination rules, exclusions and 1-4 conventions must agree with
the selected Amber parameter family. The export does not merge arbitrary
force fields or infer ligand parameters.

## Amber library and parameter reuse

`protein.lib` is a full-system unit named `protein`: coordinates, residue
membership, atom types, charges and connectivity. `parameters.frcmod` prints
ordinary mass, bonded and LJ parameters. `ATOM_TYPES.csv` maps exported type
names to source types. These files are system-specific.

The PDB can serve as a docking receptor while the full-system library remains
available as a parameterized protein unit. A ligand requires its own compatible
parameters. Subsequent complex preparation must preserve protein chemistry
and resolve any naming or atom-type conflicts. A modified PDB alone is not a
force field. No ligand-combination utility is included.

An ordinary frcmod does not contain CMAP. Preserve ff19SB through its complete
numeric topology or `system.xml`/`PARAMETERS.json`; `rebuild.py` reconstructs
Amber files without losing CMAP. Direct LEaP rebuilding is provided for systems
without CMAP. An unmodified loadPDB workflow cannot be assumed to reconstruct
every custom CMAP context.

The fitted residue templates distributed with the package are distinct from
the per-system OFF export. Format conversion never refits their charge sets.
