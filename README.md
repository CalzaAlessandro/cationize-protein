# cationize-protein

Progressive electrostatic-guided protein cationization, with cumulative
structures and complete parameters for Amber, OpenMM and GROMACS.

The program converts protein carboxyl groups into
`–C(=O)–NH–CH2–CH2–NH3+`. Each step ranks accessible sites, modifies one group,
rotates permitted bonds to reduce clashes, minimizes the structure and audits
native protonation with PROPKA. Introduced ammonium groups remain +1.

This private development version supports **ff14SB, ff19SB and ff99SBildn**.
Unsupported chemical contexts are rejected. See [methods](docs/METHODS.md)
and [validation](docs/VALIDATION.md).

## Install

Without Conda, from this repository (Python 3.12 recommended):

```sh
python3 install.py
.venv/bin/cationize-protein doctor
```

With an existing Conda installation:

```sh
conda env create -f environment.yml
conda activate cationize-protein
cationize-protein doctor
```

Both methods install the same Python package and fitted residue library.
PDBFixer is included for initial heavy-atom preparation.
Optional Amber, GROMACS and APBS programs are detected, not downloaded.
Missing external tools do not fail a complete Python installation and can be
installed later. Nothing is installed globally. See [installation](docs/INSTALL.md).

## Check and run

```sh
.venv/bin/cationize-protein check protein.pdb
.venv/bin/cationize-protein prepare protein.pdb --output prepared
.venv/bin/cationize-protein run --pdb prepared/protein.pdb --output cationized_series \
  --formats amber openmm gromacs
```

`prepare` uses PDBFixer to add missing heavy atoms and terminal oxygen atoms
within existing native residues. It saves `original.pdb`, `protein.pdb` and a
readable `preparation.log` in a new directory. Supplied atom records, coordinates,
residue states and connectivity records are preserved. It does not add whole
residues, choose alternate conformers, replace nonstandard residues or add H.
Missing sequence and ambiguous or malformed inputs stop preparation. Histidine
states must still be specified as described below. See [PDB preparation](docs/PDB_PREPARATION.md).

For an already prepared input, run directly:

```sh
.venv/bin/cationize-protein run --pdb protein.pdb --output cationized_series \
  --formats amber openmm gromacs
```

With Conda activated, use `cationize-protein` directly. New terminal runs offer
numbered choices for minimizer, solver and force field. Defaults are OpenMM,
APBS and ff14SB. Optional flags bypass the corresponding questions;
noninteractive runs use defaults or saved preferences. An unavailable requested
tool stops that calculation without silently changing the model.

The PDB check prints findings on screen. `--json report.json` additionally
writes a report. Supply complete heavy atoms, resolved alternate locations and
explicit histidine states or diagnostic ring hydrogens. `run` does not repair
missing atoms; use `prepare` first when needed. Protein-only Amber topology/restart input can be
supplied instead with `--topology FILE --coordinates FILE`.

## Main options

| Option | Purpose |
|---|---|
| `--percentages 25 50 75 100` | Default cumulative checkpoints; custom values accepted; baseline 0% always included |
| `--force-field ff14SB` | Optional choice: ff14SB, ff19SB or ff99SBildn |
| `--engine openmm` | Optional minimizer: OpenMM or installed serial Amber (`amber`, `sander`, `pmemd`) |
| `--solver apbs` | Optional solver: nonlinear APBS, PBSA, or explicit screened-Coulomb approximation |
| `--probe-radius 1.369` | Positive probe LJ half-Rmin in Å; default from the reference sodium model |
| `--ph 7.4` | Native PROPKA target pH; introduced ammonium stays +1 |
| `--protonation fixed` | Preserve supplied native states instead of per-step PROPKA |
| `--formats amber openmm gromacs` | Requested outputs; Amber default; PDB and XML/JSON always retained |
| `--apbs PATH`, `--pbsa PATH` | Explicit external solver location |
| `--keep-maps` | Retain every potential grid |

`run --help` documents all numerical settings. Probe and PB electrolyte ion
radii are distinct. Execution is serial and performs minimization, not MD;
Slurm is not required.

## Outputs and reuse

Requested percentages are written as `0/`, `25/`, `50/`, etc. Each contains
`protein.pdb`, `system.xml`, `state.json` and `exports/`. The export directory
includes the PDB, full numerical parameters and requested engine files:

| Destination | Files |
|---|---|
| Amber | `system.parm7`, `system.rst7`; full-system `protein.lib`, `parameters.frcmod` |
| OpenMM | `system.xml`, `state.json` |
| GROMACS | `system.top`, `system.gro`, `forcefield.itp`, `protein.itp` |
| All | `protein.pdb`, `PARAMETERS.json`, `ATOM_TYPES.csv`, reconstruction instructions |

The protein-only system already has its parameters. Its PDB can also be used
for docking and subsequent preparation of a complex with independently
parameterized ligands. This program does not place or parameterize ligands.
The GROMACS ITPs are separate and use a `CP_` atom-type namespace. The Amber
OFF export is a complete system unit. See [output details](docs/OUTPUTS.md).

Repeat the same run command to resume. Recorded choices and all input,
implementation, parameter and checkpoint hashes must match. Valid results are
never overwritten. A changed scientific protocol requires a new directory.

## Development

See [working with this private repository](docs/GITHUB.md). Source, tests and
small parameter files belong in Git; environments, binaries, private
calculations and MD trajectories do not. An open-source licence has not yet
been selected for this private development version.
