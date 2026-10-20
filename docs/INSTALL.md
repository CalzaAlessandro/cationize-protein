# Installation

Use Python 3.12 if available. The package supports Python 3.10–3.14 subject
to compatible dependency wheels. Choose either the local venv or Conda route.

## Local venv

From the repository folder:

```sh
python3 install.py
source .venv/bin/activate
cationize-protein doctor
```

The installer creates `.venv` and installs this package, its residue library,
NumPy, SciPy, ParmEd, OpenMM, PROPKA and PDBFixer. It does not change the system
Python or shell startup files. Activate the environment in each new terminal;
after activation, use `cationize-protein` directly. Run `deactivate` to leave it.

An unrelated existing `.venv` is preserved. Choose another directory with
`python3 install.py --venv .venv-cationize`, then activate that environment.
To select a Python interpreter, use `--python /path/to/python3.12`.
A venv cannot be moved reliably; recreate it after moving the repository.

## Conda

With Conda already installed, run from the repository folder:

```sh
conda env create -f environment.yml
conda activate cationize-protein
cationize-protein doctor
```

The environment file creates Python 3.12 and installs the same package and
Python dependencies. It does not create a nested `.venv`. A new Conda
environment has not yet been exercised in the current validation; the tested
installation uses the private venv.

## External programs

Python installation completes even if the following programs are absent:

| Program | Needed for | Search locations |
|---|---|---|
| Amber `sander` or `pmemd` | Optional serial Amber minimization | `AMBERHOME/bin`, then `PATH` |
| Amber `pbsa` | Optional Poisson–Boltzmann solver | `AMBERHOME/bin`, then `PATH` |
| APBS | Default Poisson–Boltzmann solver | Project `.tools/apbs/bin/apbs`, then `PATH` |
| GROMACS `gmx` or `gmx_mpi` | Running GROMACS calculations | `PATH` |
| MODELLER | Optional native comparative model from PDB + FASTA | `CATIONIZE_MODELLER_PYTHON`, otherwise active Python |
| ESMFold2 | Optional native prediction from FASTA | `CATIONIZE_ESMFOLD_PYTHON`, otherwise active Python |

OpenMM minimization and all three export formats work without Amber or
GROMACS executables. Gaussian and Amber parameter-development tools are only
needed to develop new fitted templates. They are not required to use the
accepted library.

Install external programs separately, then run `cationize-protein doctor`
again. Detection is repeated, so reinstalling this package is unnecessary.
MODELLER and ESMFold2 are not installed by the standard installer or Conda
environment file. Their licence/model weights are not checked by `doctor`.
See [optional modeling](MODELING.md) for official installation links, separate
Python environments and commands. Their numerical modeling runs are not yet
qualified in this package.
For platforms supported by conda-forge, APBS and GROMACS can be installed with:

```sh
conda install -c conda-forge apbs gromacs
```

Without Conda, follow the official [APBS installation guide](https://apbs.readthedocs.io/en/stable/getting/)
and [GROMACS installation guide](https://manual.gromacs.org/current/install-guide/index.html).
Amber must be obtained through its [official distribution](https://ambermd.org/).

A requested solver or engine must be available before that calculation starts.
The program never substitutes a different electrostatic model automatically.
`--solver coulomb` selects an explicit screened-Coulomb approximation.

## Choose settings when running

A new interactive `run` offers choices for minimizer, solver and force field.
Defaults are OpenMM, APBS and ff14SB. You can set options directly instead:

```sh
cationize-protein run --pdb prepared/protein.pdb --output series \
  --engine openmm --solver apbs --force-field ff14SB
```

Noninteractive runs use supplied options or saved defaults. A resumed run uses
its recorded scientific choices. Probe radius defaults to 1.369 Å.

## Check an installation without changing it

```sh
python3 install.py --check
cationize-protein doctor
```

`install.py --check` does not install packages. `doctor` reports Python
packages and external executables. The installer returns success when the
Python runtime is complete; a missing selected external solver is reported
separately. Optional JSON output is available with `--json`.

The installation environment is in `.venv`; any separately installed private
APBS bundle is under `.tools/apbs`. Inputs and calculation outputs are stored
outside these directories. Preserve unrelated files when removing an
installation.
