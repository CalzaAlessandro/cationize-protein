# Install optional software

The standard installation already includes OpenMM, PDBFixer, PROPKA, ParmEd
and the ff14SB parameters. It can minimize and rank sites with the included
screened-Coulomb approximation, and export Amber, OpenMM and GROMACS files.
Install the tools below only for the additional calculation you need.

[APBS](#apbs) · [AmberTools / Amber](#ambertools-and-amber) ·
[GROMACS](#gromacs) · [MODELLER](#modeller) · [ESMFold2](#esmfold2)

| Software | Use in this package |
|---|---|
| APBS | Poisson-Boltzmann potential, `--solver apbs` |
| AmberTools | Serial `sander` minimization and `pbsa` potential |
| Amber | Optional serial `pmemd` minimization |
| GROMACS | Run exported GROMACS systems |
| MODELLER | Build a native model from an incomplete PDB and complete FASTA |
| ESMFold2 | Predict a native model from FASTA with local model weights |

These are separate installations. Installing this Python package with `pip`
does not install these programs. The commands below use macOS/Linux shells;
on a cluster, check for an existing software module first.

## Keep the main environment available

Start in the repository folder with the main environment active, and record
its absolute path:

```sh
export CATIONIZE_REPO="$PWD"
mkdir -p "$CATIONIZE_REPO/.tools"
cationize-protein doctor
```

The examples install optional Conda environments by absolute prefix, without
activating them. The modelers are called through their own Python, so the main
environment keeps providing `cationize-protein`.

If you activate another environment while following upstream instructions,
leave it first. For a main venv, use `conda deactivate` or `deactivate` as
appropriate, then `source "$CATIONIZE_REPO/.venv/bin/activate"`. For a main
Conda installation, use `conda activate cationize-protein`. Replace `.venv`
if you chose a different main environment directory.

Exports apply to the current terminal. Repeat the relevant exports in a new
terminal, using the absolute repository path, or add them to your shell
configuration after checking them.

## APBS

**With Conda, on supported Linux/macOS platforms:** install the
[conda-forge APBS package](https://anaconda.org/conda-forge/apbs):

```sh
conda create --prefix "$CATIONIZE_REPO/.tools/apbs-env" \
  -c conda-forge --strict-channel-priority apbs
export PATH="$CATIONIZE_REPO/.tools/apbs-env/bin:$PATH"
```

**Without Conda:** download a binary archive matching your OS and architecture
from the [official APBS releases](https://github.com/Electrostatics/apbs/releases),
then unpack it. Keep its `bin`, `lib` and `share` directories together in a
directory such as `$HOME/.local/apbs`. Follow the
[APBS binary installation notes](https://apbs.readthedocs.io/en/stable/getting/)
for any required library path. Then connect that installation:

```sh
export PATH="$HOME/.local/apbs/bin:$PATH"
```

If no matching binary is available, use the
[official source build](https://apbs.readthedocs.io/en/latest/getting/source.html).
For example, after installing CMake and the compiler dependencies listed there:

```sh
git clone --branch v3.4.1 --depth 1 https://github.com/Electrostatics/apbs.git \
  "$CATIONIZE_REPO/.tools/apbs-src"
cmake -S "$CATIONIZE_REPO/.tools/apbs-src" \
  -B "$CATIONIZE_REPO/.tools/apbs-build" \
  -DCMAKE_INSTALL_PREFIX="$CATIONIZE_REPO/.tools/apbs" \
  -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTS=ON
cmake --build "$CATIONIZE_REPO/.tools/apbs-build" --parallel 4
PATH="$CATIONIZE_REPO/.tools/apbs-build/bin:$PATH" \
  ctest --test-dir "$CATIONIZE_REPO/.tools/apbs-build" --output-on-failure
cmake --install "$CATIONIZE_REPO/.tools/apbs-build"
export PATH="$CATIONIZE_REPO/.tools/apbs/bin:$PATH"
```

Check the installed executable:

```sh
command -v apbs
apbs --version
cationize-protein doctor
```

Use `--solver apbs` on `run`; `--apbs /absolute/path/to/apbs` selects one
installation explicitly. The project also checks `.tools/apbs/bin/apbs`.

## AmberTools and Amber

AmberTools supplies both `sander` and `pbsa`. The
[official AmberTools download page](https://ambermd.org/GetAmber.php) documents
the Conda binary and source distributions.

**With Conda, on supported Linux/macOS platforms:**

```sh
conda create --prefix "$CATIONIZE_REPO/.tools/ambertools-env" \
  -c conda-forge --strict-channel-priority python=3.12 ambertools=26
export AMBERHOME="$CATIONIZE_REPO/.tools/ambertools-env"
source "$AMBERHOME/amber.sh"
```

This uses the official conda-forge AmberTools26 package. It provides serial
executables suitable for this package; it does not include the separately
licensed Amber `pmemd` engine.

**Without Conda:** download AmberTools source from that official page, install
the C/C++/Fortran compilers required for your OS, and follow Chapter 2 of the
[Amber manual](https://ambermd.org/Manuals.php) and the
[OS-specific installation instructions](https://ambermd.org/Installation.php).
After unpacking the source archive, edit its `build/run_cmake` for your
compiler and writable installation prefix, then configure and install:

```sh
cd /absolute/path/to/amber26_src/build
./run_cmake
make install
```

Check the build output, run the serial tests described in the manual, and
source `amber.sh` from the chosen installation prefix:

```sh
export AMBERHOME="/absolute/path/to/amber-installation"
source "$AMBERHOME/amber.sh"
```

For `pmemd`, obtain and build the separately licensed
[Amber distribution](https://ambermd.org/GetAmber.php), or load the cluster's
Amber module. This interface uses serial `sander` or `pmemd`; a `pmemd.cuda`
installation alone does not provide its Amber minimization engine.

From the main environment, check:

```sh
command -v sander
command -v pbsa
cationize-protein doctor
```

Choose `--engine sander` for AmberTools minimization, or `--solver pbsa` for
Poisson-Boltzmann potential. `--engine amber` selects available serial Amber
software automatically. `--engine-executable /absolute/path/to/sander` and
`--pbsa /absolute/path/to/pbsa` select executables explicitly. Discovery checks
`AMBERHOME/bin` before `PATH`.

## GROMACS

**With Conda:** install the
[conda-forge GROMACS package](https://anaconda.org/conda-forge/gromacs):

```sh
conda create --prefix "$CATIONIZE_REPO/.tools/gromacs-env" \
  -c conda-forge --strict-channel-priority gromacs
export PATH="$CATIONIZE_REPO/.tools/gromacs-env/bin:$PATH"
```

**Without Conda, on macOS with Homebrew:** the
[Homebrew formula](https://formulae.brew.sh/formula/gromacs) provides `gmx`:

```sh
brew install gromacs
source "$(brew --prefix)/share/gromacs/GMXRC"
```

**Linux or a custom CPU/GPU build:** download and unpack a release from
[GROMACS](https://manual.gromacs.org/current/download.html), then follow its
[installation guide](https://manual.gromacs.org/current/install-guide/index.html).
For a basic CPU build, with the source unpacked at `.tools/gromacs-src` and
the compiler/CMake requirements from that guide installed:

```sh
cmake -S "$CATIONIZE_REPO/.tools/gromacs-src" \
  -B "$CATIONIZE_REPO/.tools/gromacs-build" \
  -DCMAKE_INSTALL_PREFIX="$CATIONIZE_REPO/.tools/gromacs" \
  -DGMX_GPU=OFF -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
cmake --build "$CATIONIZE_REPO/.tools/gromacs-build" --parallel 4
cmake --build "$CATIONIZE_REPO/.tools/gromacs-build" --target check
cmake --install "$CATIONIZE_REPO/.tools/gromacs-build"
source "$CATIONIZE_REPO/.tools/gromacs/bin/GMXRC"
```

Check with `gmx --version` and `cationize-protein doctor`. An MPI installation
may provide `gmx_mpi`; discovery accepts either name. GPU support depends on
the selected GROMACS build. GROMACS is used to run exported systems; exporting
their TOP/GRO/ITP files is included in the standard Python installation.

## MODELLER

First obtain an individual licence through
[MODELLER registration](https://salilab.org/modeller/registration.html).
Use the official structural modeling distribution. The unrelated PyPI package
named `modeller` does not provide it.

**With Conda:** create a separate Python environment using the official Sali
Lab channel, then configure the licence file as instructed by the installer:

```sh
conda create --prefix "$CATIONIZE_REPO/.tools/modeller-env" \
  -c salilab -c conda-forge python=3.12 modeller
export CATIONIZE_MODELLER_PYTHON="$CATIONIZE_REPO/.tools/modeller-env/bin/python"
```

**Without Conda, on macOS with Homebrew:** use the
[official MODELLER Homebrew instructions](https://salilab.org/modeller/10.8/release.html#mac).
Complete the licence setup printed by Homebrew:

```sh
brew tap salilab/salilab
brew install modeller
export CATIONIZE_MODELLER_PYTHON="$(brew --prefix python@3.14)/bin/python3.14"
```

The current [Sali Lab formula](https://github.com/salilab/homebrew-salilab/blob/master/modeller.rb)
uses Python 3.14; match its Python version if the formula changes.
On Linux, download the official `.deb` or RPM for your architecture and follow
its licence/install command on the
[MODELLER installation page](https://salilab.org/modeller/download_installation.html).
The generic Unix installer is also available for installations without root
access; its Python environment must be configured before use here. For a
native Linux installation, set `CATIONIZE_MODELLER_PYTHON` to the absolute
Python path that imports MODELLER, for example `/usr/bin/python3`.

Verify using the selected Python, then return to the main environment if needed:

```sh
"$CATIONIZE_MODELLER_PYTHON" -c 'from modeller import Environ; Environ(); print("MODELLER licence OK")'
cationize-protein doctor
```

This variable selects a Python executable that imports MODELLER, not the
`mod10.x` launcher. A real modeling command is shown in [MODELING.md](MODELING.md).

## ESMFold2

Use a separate Python 3.12 environment. Install Python 3.12 from
[Python.org](https://www.python.org/downloads/) if needed, or use the Conda
route below. The commands use the package versions from the recorded
short-fragment CUDA validation: PyTorch 2.11.0,
[Transformers 5.20.1](https://pypi.org/project/transformers/5.20.1/) and
[Accelerate 1.15.0](https://pypi.org/project/accelerate/1.15.0/). The adapter
uses the [Transformers ESMFold2 interface](https://huggingface.co/docs/transformers/main/model_doc/esmfold2).

Create the environment using **one** of these routes:

```sh
# Without Conda
python3.12 -m venv "$CATIONIZE_REPO/.tools/esmfold2-env"
```

```sh
# With Conda
conda create --prefix "$CATIONIZE_REPO/.tools/esmfold2-env" python=3.12 pip
```

Point to its Python; activation is unnecessary:

```sh
export CATIONIZE_ESMFOLD_PYTHON="$CATIONIZE_REPO/.tools/esmfold2-env/bin/python"
"$CATIONIZE_ESMFOLD_PYTHON" -m pip install --upgrade pip
```

Install **one** PyTorch build, following the
[official PyTorch version instructions](https://pytorch.org/get-started/previous-versions/#v2110)
for your hardware. For Linux with an NVIDIA GPU and a driver compatible with
CUDA 12.8:

```sh
"$CATIONIZE_ESMFOLD_PYTHON" -m pip install torch==2.11.0 \
  --index-url https://download.pytorch.org/whl/cu128
```

For Linux CPU execution:

```sh
"$CATIONIZE_ESMFOLD_PYTHON" -m pip install torch==2.11.0 \
  --index-url https://download.pytorch.org/whl/cpu
```

For macOS on Apple Silicon, use its standard wheel and select `--device cpu`
when modeling; this adapter does not expose the MPS device:

```sh
"$CATIONIZE_ESMFOLD_PYTHON" -m pip install torch==2.11.0
```

Other platforms need an available compatible PyTorch wheel; this pinned
profile does not cover every OS/architecture. Finish the Python installation:

```sh
"$CATIONIZE_ESMFOLD_PYTHON" -m pip install "transformers==5.20.1" "accelerate==1.15.0"
"$CATIONIZE_ESMFOLD_PYTHON" -m pip check
"$CATIONIZE_ESMFOLD_PYTHON" -c 'import torch; from transformers import EsmFold2Model; print(torch.__version__, "CUDA:", torch.cuda.is_available())'
cationize-protein doctor
```

The import check does not load weights. To download/cache the public
[ESMFold2-Fast checkpoint](https://huggingface.co/biohub/ESMFold2-Fast) and run
local prediction, explicitly request a model from the main environment:

```sh
cationize-protein model --method esmfold --fasta mature.fasta \
  --esmfold-model fast --device cuda --output predicted
```

Use `--device cpu` for a CPU installation. CPU prediction can be slow and
memory intensive; assess the actual sequence length. `--esmfold-model full`
selects [ESMFold2](https://huggingface.co/biohub/ESMFold2). The program records
the checkpoint revision; `--model-revision` can pin a published commit.

`doctor` checks installed code. It does not download weights, validate a
MODELLER licence or prove GPU execution. ESMFold2 weights are fetched only on
the explicit modeling command above; this local route uses no Biohub API key
and sends no sequence to a remote folding service. See [MODELING.md](MODELING.md)
for model checks, output files and current validation limits.

The installation recipes follow the linked upstream instructions. Their
platform coverage extends beyond the environments tested by this project;
the functional validation reported in [MODELING.md](MODELING.md) remains the
evidence for successful model execution.
