# Installation with or without Conda

Both methods install the complete Python package, its runtime dependencies and
the fitted residue library, including PDBFixer for initial PDB preparation.
Amber, GROMACS, APBS and Gaussian are separate
external programs. Missing external programs do not invalidate a successful
Python installation. The program checks their current availability when needed;
installing one later does not require reinstalling this package.

## Without Conda

Use a downloaded or cloned checkout. Python 3.10–3.14 is supported by the current runtime wheels; Python 3.12 is recommended. The installer searches for Python 3.12, 3.13, 3.11, 3.14 and 3.10, then the launching interpreter if compatible. From the project folder:

```sh
python3 install.py
```

The installer creates `.venv` and runs that environment's `python -m pip` to install the project and its Python dependencies: NumPy, SciPy, ParmEd, OpenMM and PROPKA. It does not replace the Python that started it, use Conda, run `sudo`, modify shell startup files, or install global tools. This follows Python's [isolated virtual environment](https://docs.python.org/3/library/venv.html) mechanism and pip's [interpreter-specific invocation](https://pip.pypa.io/en/stable/user_guide/#running-pip).

Installation does not ask for a scientific model. New terminal runs ask for
the minimizer, electrostatic solver and protein force field. The optional
installer flags below only save preferences for later commands:

```sh
python3 install.py --engine openmm --force-field ff14SB
```

The default potential solver is APBS. Engine and force field are separate choices. The registry also exposes `ff19SB` and `ff99SBildn`; the report checks their availability in the selected environment and records the registry's validation scope. Availability of a native protein force field does not validate modified residues or an entire simulation. All custom RESP templates must still pass the workflow's chemistry and content-hash checks. The selected settings are saved in `.venv/cationize-protein-settings.json` for the runtime CLI.

An existing, unrelated `.venv` is preserved and refused. Use a fresh project-local environment, for example `--venv .venv-cationize`, or keep the existing environment and install there manually. A relative environment path is interpreted from this checkout, even if the script is run from a different working directory. Symlinks and paths overlapping the source, parameter library or other project data are rejected. An interrupted creation keeps its partial directory for inspection; remove or rename that specific directory before retrying if its interpreter is incomplete.

To choose an existing Python explicitly:

```sh
python3 install.py --python /path/to/python3.12
```

If Python/venv support is absent or no compatible Python wheels are available, the installer reports the limitation. It does not silently install a new system Python or build a global toolchain. pip errors are summarized because pip output can contain authenticated package-index URLs. No raw pip output or credentials are printed.

## Read-only inspection

```sh
python3 install.py --check
python3 install.py --check --json
```

`--check` never installs or creates an environment and never prompts. If `.venv` exists, it inspects that interpreter; otherwise it reports the missing private environment and shows package discovery in the selected bootstrap Python. For this explicit readiness check, exit code 0 means the selected workflow dependencies were detected, 1 means incomplete/unavailable, and 2 means invalid command syntax. An actual installation returns 0 when the complete Python runtime is installed, even when optional external programs are absent. Its JSON output records this separately as `installation_complete`; `ready` still describes the selected scientific workflow.

From either an activated Conda environment or the private environment, inspect
the currently running installation without creating files:

```sh
cationize-protein doctor
cationize-protein doctor --json
```

Use `.venv/bin/cationize-protein` if the private environment is not activated.
`doctor` returns 0 when the Python runtime is complete; unavailable external
tools remain visible in the report. Executable discovery is not numerical
validation.

## With Conda

From the downloaded or cloned project folder:

```sh
conda env create -f environment.yml
conda activate cationize-protein
cationize-protein doctor
cationize-protein check native.pdb
cationize-protein run --pdb native.pdb --output cationized_series
```

`environment.yml` creates Python 3.12 and pip, then installs this project and
all its Python requirements into that Conda environment. No nested `.venv` is
needed. The command must be run from the project folder because the pip entry
is `.`. Conda itself must already be installed. See the official
[environment documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

Optional open-source external programs can be installed into the same Conda
environment if packages are available for the platform:

```sh
conda install -c conda-forge apbs gromacs
cationize-protein doctor
```

The official conda-forge recipes are
[APBS](https://github.com/conda-forge/apbs-feedstock) and
[GROMACS](https://github.com/conda-forge/gromacs-feedstock). This is an explicit
user choice; the Python installer does not install them automatically. Amber
and Gaussian must be obtained separately under their applicable licences.

A check can use the same selections as an installation, for example:

```sh
python3 install.py --check --engine amber --solver pbsa --force-field ff19SB
```

## Amber and other external software

The checker looks for `pmemd.cuda`, `pmemd`, `sander`, `tleap`, `pbsa`, `antechamber` and `resp` in `AMBERHOME/bin`, then on `PATH`. It also checks `gmx`/`gmx_mpi` and APBS. Their presence is reported separately: a topology builder or PBSA executable is not automatically a minimizer. Amber minimization requires one of the supported serial programs; requested programs are never replaced silently. GROMACS is optional and is not needed to export `.top`/`.gro` files. Executable discovery does not establish a valid licence, GPU availability or successful computation.

No external engine or solver is downloaded during ordinary installation. If an explicitly selected Amber engine or PBSA solver is missing, installation reports it after installing the private Python runtime. Obtain and configure those programs separately through the applicable official distribution. Installing OpenMM does not install Amber.

## APBS

The search checks `.tools/apbs/bin/apbs` in this checkout and then `PATH`. If
APBS is absent, installation reports that fact and still completes the Python
package. Install a compatible APBS separately, through Conda or the official
release, then run `doctor`. A scientific run selecting APBS stops until its
executable is available; it never substitutes Coulomb automatically.

The existing private downloader is retained as a development utility, with
pinned official APBS 3.4.1 assets in `config/install-manifest.json`; ordinary
installation does not call it. It verifies SHA256 before safe extraction,
retains companion libraries, and does not overwrite an existing installation.

Official APBS recommends [release binaries](https://apbs.readthedocs.io/en/stable/getting/) and publishes them on [GitHub releases](https://github.com/Electrostatics/apbs/releases/tag/v3.4.1). The Linux release has documented glibc requirements. Apple Silicon, unsupported platforms or incompatible libraries may require a separate compatible installation or an [official local source build](https://apbs.readthedocs.io/en/stable/getting/source.html). Those cases are reported honestly; the installer does not use an unofficial `apbs-binary` wheel or Homebrew. Please follow APBS's registration and citation requests in its official documentation.

An APBS executable being found is not an end-to-end Poisson–Boltzmann validation. The workflow must check that the selected solver actually produces usable potential data. PBSA is a separate explicitly selected option. `--solver coulomb` requests a screened Coulomb approximation and is labelled as such; it is never an automatic fallback for a missing Poisson–Boltzmann solver.

## Running and removing the private installation

```sh
.venv/bin/python -m cationize_protein --help
```

On Windows, the Python interpreter is `.venv/Scripts/python.exe`. The pinned APBS download currently covers only the platforms listed above.

The private Python environment is confined to `.venv`; an existing private
APBS bundle may occupy `.tools/apbs`. Preserve unrelated files in `.tools`.
Scientific inputs and outputs remain separate. A virtual environment is not
relocatable: after copying the checkout, create a fresh environment there.
