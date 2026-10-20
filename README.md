# cationize-protein

Build progressively cationized proteins and export their structures and
parameters for Amber, OpenMM or GROMACS.

<p align="center">
  <img src="docs/images/hsa_progressive.gif" width="560" alt="Five aligned minimized HSA structures, with progressively more modified residues shown as sticks">
</p>

*HSA example: five minimized structures, from 0% to 100% cationization.
Modified residues are shown as sticks. This animation switches between
chemical checkpoints.*

The program converts carboxyl groups into `–C(=O)–NH–CH₂–CH₂–NH₃⁺`.
It selects one site at a time from the electrostatic potential, rotates the
new group to reduce clashes, minimizes the structure and checks protonation
with PROPKA. Each checkpoint retains the modifications made earlier.

[Install](#install) · [Prepare a protein](#prepare-a-protein) ·
[Run cationization](#run-cationization) · [Find your files](#find-your-files) ·
[HSA example](#hsa-example) · [Options](#options)

## Install

Python 3.12 is recommended. From a terminal:

```sh
git clone https://github.com/CalzaAlessandro/cationize-protein.git
cd cationize-protein
python3 install.py
source .venv/bin/activate
cationize-protein doctor
```

The installer creates a local environment and installs the Python dependencies,
including OpenMM, PDBFixer and PROPKA. Activate it with
`source .venv/bin/activate` whenever you open a new terminal. After activation,
use `cationize-protein` directly. Use `deactivate` to leave the environment.

If you prefer Conda, replace the last three commands with:

```sh
conda env create -f environment.yml
conda activate cationize-protein
cationize-protein doctor
```

Amber, GROMACS and APBS are separate programs. Installation reports whether
they are present; you can add them later. APBS or Amber PBSA is needed for
Poisson–Boltzmann calculations. GROMACS is only needed to run GROMACS jobs,
not to write its input files. See [installation details](docs/INSTALL.md).

## Prepare a protein

```sh
cationize-protein prepare input.pdb --output prepared
```

This checks the input, uses PDBFixer to add missing atoms within existing
residues, and assigns histidine charge states with PROPKA at pH 7.4.
OpenMM selects a neutral HID/HIE tautomer when the input does not specify one.
The original PDB is preserved.

```text
prepared/
├── original.pdb
├── protein.pdb
├── preparation.log
├── histidines.csv       When histidines are present
└── protonation/         Detailed PROPKA evidence
```

Missing loops, alternate conformations and malformed coordinate records need
to be resolved before proceeding. The terminal explains what must be fixed.
For inspection without changes, use `cationize-protein check input.pdb`.
See [PDB preparation](docs/PDB_PREPARATION.md) and the
[broken-input examples](examples/pdb_checks/README.md).

## Run cationization

```sh
cationize-protein run --pdb prepared/protein.pdb --output series \
  --formats amber openmm gromacs
```

The default series is **0%, 25%, 50%, 75% and 100%**. Percentages refer to
modifiable carboxyl groups and are rounded to whole groups. The baseline 0%
is always saved. Introduced ammonium groups remain +1 throughout the protocol.

A new terminal run asks you to choose the minimization engine, electrostatic
solver and force field. The defaults are OpenMM, APBS and ff14SB. Available
protein force fields are ff14SB, ff19SB and ff99SBildn. Options can also be
supplied on the command line, as shown below.

Run the same command again to resume an interrupted calculation. The program
checks the saved inputs and parameters before continuing. Use a new output
directory when changing the protocol. Execution is serial; Slurm is optional.
The calculation ends after minimization, before equilibration or MD production.

## Find your files

```text
series/
├── residue_templates/    Individual residue libraries: ASC.lib, GLC.lib, …
├── 0/
├── 25/
│   ├── protein.pdb      Minimized structure with the actual residue names
│   └── exports/         Parameters and files for the requested engines
├── 50/
├── 75/
├── 100/
└── internal/            Intermediate steps and recovery records
```

| What you need | Where to find it |
|---|---|
| Structure for inspection or docking | `25/protein.pdb` |
| Ready Amber protein | `25/exports/system.parm7` and `system.rst7` |
| Ready OpenMM protein | `25/exports/system.xml` and `state.json` |
| Ready GROMACS protein | `25/exports/system.top`, `system.gro`, `forcefield.itp`, `protein.itp` |
| One modified residue | `residue_templates/ASC.lib`, `GLC.lib`, etc. |
| Complete Amber protein unit | `25/exports/protein.lib` |
| Printed protein parameters | `25/exports/parameters.frcmod` and `PARAMETERS.json` |

The individual `.lib` files use Amber **OFF format**; they do not need duplicate
`.off` copies. `protein.lib` contains the entire protein, whereas `ASC.lib`
contains one residue template. Match the library to the residue's chemical
and terminal state. See [output and parameter reuse](docs/OUTPUTS.md).

## HSA example

[examples/HSA](examples/HSA/README.md) contains the existing **585-residue HSA
series** at 0%, 25%, 50%, 75% and 100%, with a PDB, Amber topology and restart
for each minimized structure. Tables list the modified sites and protein
charges. The 17 disulfide bonds are retained at every level.

![Aligned views of the five minimized HSA structures](docs/images/hsa_comparison.png)

To try the program without APBS, use the small peptide example:

```sh
cationize-protein run --pdb examples/pdb_checks/reference_peptide.pdb \
  --output peptide_series --engine openmm --solver coulomb \
  --force-field ff14SB --formats amber openmm gromacs
```

This explicitly uses a screened-Coulomb approximation. Use APBS or PBSA for
a Poisson–Boltzmann calculation.

## Options

| Option | Purpose |
|---|---|
| `--percentages 25 28 35` | Save a custom cumulative series |
| `--ph 7.4` | Target pH; use the same value in `prepare` and `run` |
| `--engine openmm` | Minimize with OpenMM; `amber`, `sander` and `pmemd` select installed serial Amber |
| `--solver apbs` | Choose `apbs`, `pbsa` or the `coulomb` approximation |
| `--force-field ff14SB` | Choose ff14SB, ff19SB or ff99SBildn |
| `--formats amber openmm gromacs` | Choose output formats; Amber is the default |
| `--probe-radius 1.369` | Positive probe LJ half-Rmin in Å |
| `--protonation fixed` | Retain supplied protonation states |

Use `cationize-protein --help`, `prepare --help` or `run --help` for the full
command descriptions. Advanced settings, installation troubleshooting and
scientific qualification are documented in
[installation](docs/INSTALL.md), [methods](docs/METHODS.md),
[force fields](docs/FORCE_FIELDS.md) and [validation](docs/VALIDATION.md).

This repository is private during development. See [GitHub instructions](docs/GITHUB.md)
for cloning, saving changes and granting collaborator access.
