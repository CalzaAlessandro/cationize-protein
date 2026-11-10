# Quickstart

Amidate the three carboxyl groups of the bundled Asp–Glu–Leu reference peptide.
From the repository root, create a Python 3.12 environment and install the package:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install .
python examples/quickstart/run_quickstart.py --output-dir work/quickstart
```

The script checks [`reference_peptide.pdb`](../examples/pdb_checks/reference_peptide.pdb),
runs progressive amidation and checks the exported files. It uses OpenMM
minimization and the analytical screened-Coulomb solver. Separate APBS, Amber
and GROMACS installations are unnecessary for this example.

To run the two commands directly:

```bash
cationize-protein check examples/pdb_checks/reference_peptide.pdb
cationize-protein run --pdb examples/pdb_checks/reference_peptide.pdb --output work/peptide_cationized --engine openmm --solver coulomb --force-field ff14SB --formats amber openmm gromacs
```

Expected results:

| Checkpoint | Amidated groups | Atoms | Formal charge |
|---|---:|---:|---:|
| `0/` | 0 | 46 | −2 e |
| `25/` | 1 | 57 | 0 e |
| `50/` | 2 | 68 | +2 e |
| `75/` | 2 | 68 | +2 e |
| `100/` | 3 | 79 | +4 e |

Three groups cannot represent every percentage exactly. The 50% and 75%
checkpoints therefore contain the same chemical state.

Each checkpoint contains `protein.pdb`. The final `100/exports/` directory
contains `system.xml` (OpenMM), `system.parm7` and `system.rst7` (Amber),
and `system.top`, `system.gro`, `protein.itp` and `forcefield.itp` (GROMACS).
`RUN.json` records the settings and modification history.

This peptide checks the software workflow. It does not validate protein
structures, amidation kinetics or force-field transferability. See the
[validation notes](VALIDATION.md) before running a scientific study.
