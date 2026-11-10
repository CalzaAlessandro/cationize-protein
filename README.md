<h1 align="center">cationize-protein</h1>

Build progressively cationized proteins and export their structures and
parameters for Amber, OpenMM or GROMACS.

<p align="center">
  <a href="https://github.com/CalzaAlessandro/cationize-protein/actions/workflows/tests.yml"><img src="https://img.shields.io/github/actions/workflow/status/CalzaAlessandro/cationize-protein/tests.yml?branch=main&amp;style=flat&amp;label=Tests" alt="Software tests"></a>
  <img src="https://img.shields.io/badge/Status-Active%20development-9A6700?style=flat" alt="Status: active development">
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/Python-3.10%E2%80%933.14-58616B?style=flat" alt="Python: 3.10 to 3.14"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache--2.0-607D8B?style=flat" alt="License: Apache-2.0"></a>
</p>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/banner_dark.svg">
  <img src="docs/images/banner.svg" width="100%" alt="A schematic folded chain with negatively charged sites progressively converted into positive sites; the aminoethyl-amide chemical scheme is shown below">
</picture>

> [!NOTE]
> Parameters, charges and topologies are preliminary.
> See the [validation notes](docs/VALIDATION.md) before using them in simulations.

### Quick navigation

- **Getting started:** [Quickstart](docs/QUICKSTART.md) · [How it works](#how-it-works) · [Install](#install) · [Optional external software](#optional-external-software)
- **Workflow:** [Get a structure or sequence](#get-a-structure-or-sequence) · [Prepare a protein](#prepare-a-protein) · [Run cationization](#run-cationization)
- **Outputs & ligands:** [Find your files](#find-your-files) · [Add a ligand](#add-a-ligand)
- **Examples:** [HSA example](#hsa-example) · [Lysozyme walkthrough](#lysozyme-walkthrough)
- **Reference:** [CLI options](#options) · [Software and references](#software-and-references) · [Residue parameters](#residue-parameters) · [Acknowledgements](#acknowledgements) · [License](#license)

## How it works

Chemical cationization reacts carboxylic acid groups (`Asp` and `Glu` side
chains and the biological C-terminus) with ethylenediamine to form positively
charged aminoethyl-amides (`-C(=O)-NH-CH₂-CH₂-NH₃⁺`), shifting each modified
site from `-1` to `+1` (`Δq = +2 e`).

At each step, `cationize-protein` repeats four operations until all requested
modification percentages (by default `0%`, `25%`, `50%`, `75%` and `100%`) are
reached:

1. **Electrostatic ranking:** computes the electrostatic potential around every
   unmodified carboxyl group and selects the site with the most negative
   accessible potential.
2. **Amidation & steric search:** replaces the selected carboxylate with the
   corresponding validated RESP aminoethyl-amide template and rotates torsional
   bonds to remove steric clashes while preserving the protein backbone and
   disulfides.
3. **Energy minimization:** relaxes the modified structure in implicit solvent
   (GBn2) with positional restraints on the native backbone.
4. **Protonation check:** re-evaluates native titratable residues (`His`, `Asp`,
   `Glu`, `Cys`, `Lys`, `Tyr`) with PROPKA at the target pH while keeping
   introduced ammonium groups protonated (`+1`).

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/hsa_workflow_dark.gif">
  <img src="docs/images/hsa_workflow.gif" width="100%" alt="Text-only colored cycle at left; real HSA structures and electrostatic surfaces for each percentage, alternating with fixed illustrations of amidation, minimization and protonation">
</picture>

*HSA structures and electrostatic surfaces are shown for each cationization
level. Amidation, minimization and protonation use fixed illustrative examples.
[View the image details](examples/HSA/README.md#figure).*

---

## Install

The core installation is self-contained and includes **OpenMM**, **PDBFixer**,
**PROPKA**, **ParmEd**, the screened-Coulomb electrostatic solver, and the
validated **ff14SB / ff19SB / ff99SBildn** residue libraries. Python 3.12 is
recommended ([full installation guide](docs/INSTALL.md)).

```sh
# Option A: Python with a local venv
git clone https://github.com/CalzaAlessandro/cationize-protein.git
cd cationize-protein
python3 install.py
source .venv/bin/activate
cationize-protein doctor
```

<details>
<summary><b>Option B: Install with Conda</b></summary>

```sh
git clone https://github.com/CalzaAlessandro/cationize-protein.git
cd cationize-protein
conda env create -f environment.yml
conda activate cationize-protein
cationize-protein doctor
```

Activate with `conda activate cationize-protein` in each new terminal (`conda deactivate` to leave).

</details>

---

## Optional external software

You do **not** need any external software to run the default workflow or to
export Amber, OpenMM and GROMACS topologies. Install external programs only if
you specifically need Poisson-Boltzmann electrostatics (`APBS` or `PBSA`),
serial Amber minimization (`sander` / `pmemd`), or missing-residue modeling
(`MODELLER` / `ESMFold2`).

Step-by-step installation and connection instructions for each tool are in
[**docs/OPTIONAL_SOFTWARE.md**](docs/OPTIONAL_SOFTWARE.md):

| Program | When you need it | How to select it in `cationize-protein` | Setup guide |
|---|---|---|---|
| **APBS** | Nonlinear Poisson-Boltzmann site ranking | `run --solver apbs` | [Install APBS](docs/OPTIONAL_SOFTWARE.md#apbs) |
| **AmberTools / Amber** | PBSA electrostatics or serial `sander`/`pmemd` minimization | `run --solver pbsa` or `run --engine amber` | [Install AmberTools](docs/OPTIONAL_SOFTWARE.md#ambertools-and-amber) |
| **GROMACS** | Run MD simulations on exported `.top`/`.gro` files | Exports work without GROMACS installed | [Install GROMACS](docs/OPTIONAL_SOFTWARE.md#gromacs) |
| **MODELLER** | Rebuild missing loops/residues from a template PDB + FASTA | `model --method modeller` | [Install MODELLER](docs/OPTIONAL_SOFTWARE.md#modeller) |
| **ESMFold2** | Predict a 3D structure directly from a FASTA sequence | `model --method esmfold` | [Install ESMFold2](docs/OPTIONAL_SOFTWARE.md#esmfold2) |

---

## Get a structure or sequence

Skip this step if you already have your input files. `fetch` downloads records
from RCSB PDB and UniProt ([details](docs/FETCH.md)):

- **Have only a PDB ID?** `fetch --pdb 1AKI` downloads `1AKI.pdb`, extracts
  any `SEQADV` sequence mutations/differences into `pdb_mutations.csv` and
  `DBREF` UniProt cross-references into `uniprot_candidates.csv`, and—when the
  PDB maps to a single UniProt entry—automatically downloads `uniprot.json`,
  `sequence.fasta`, `mature.fasta` and `disulfides.csv`.
- **Have only a UniProt ID?** `fetch --uniprot P02768` downloads the UniProt
  sequence, disulfide annotations (`uniprot.json`, `disulfides.csv`) and a list
  of all experimental PDB structures linked to that protein (`pdb_candidates.csv`,
  including method, resolution, chain bounds, engineered mutations and ligands)
  so you can pick which PDB to download.

```sh
# Download from PDB ID (auto-fetches unique mapped UniProt record)
cationize-protein fetch --pdb 1AKI --output downloaded

# Or specify only UniProt ID (lists available PDBs in pdb_candidates.csv), or both
cationize-protein fetch --pdb 1N5U --uniprot P02768 --output downloaded_hsa
```

## Prepare a protein

```sh
cationize-protein prepare input.pdb --output prepared
```

By default, `prepare` checks the PDB, repairs missing heavy atoms in existing
residues with PDBFixer, and assigns histidine states with PROPKA at pH 7.4
while preserving `original.pdb` ([details](docs/PDB_PREPARATION.md)). Use
`cationize-protein check input.pdb` for a read-only check.

<details>
<summary><b>Verify disulfides from UniProt or reconstruct missing residues (MODELLER / ESMFold2)</b></summary>

#### Disulfide verification from UniProt

Map the PDB sequence to a downloaded `uniprot.json` record and write verified
pairs as `CYX`, `SSBOND` and `CONECT` records:

```sh
cationize-protein prepare input.pdb --uniprot-record downloaded/uniprot.json \
  --output prepared_with_disulfides
```

#### Optional structure modeling for missing loops/residues

If whole residues are missing, build a complete single-chain native structure
first, then pass its `protein.pdb` to `prepare` ([details](docs/MODELING.md)):

```sh
# MODELLER: incomplete native PDB + complete FASTA sequence
cationize-protein model --method modeller --pdb incomplete.pdb \
  --fasta mature.fasta --output modeled

# ESMFold2-Fast: predict from complete FASTA sequence (--esmfold-model full for ESMFold2)
cationize-protein model --method esmfold --fasta mature.fasta --output predicted

# Prepare the resulting model
cationize-protein prepare modeled/protein.pdb --output prepared
```

</details>

## Run cationization

```sh
cationize-protein run --pdb prepared/protein.pdb --output series
```

| Setting | Default | Alternatives |
|---|---|---|
| Saved levels | `0%`, `25%`, `50%`, `75%`, `100%` | `--percentages 10 25 60 100` |
| Minimization engine | `openmm` (GBn2) | `--engine amber` (`sander` / `pmemd`) |
| Electrostatic solver | `coulomb` (screened Coulomb) | `--solver apbs`, `--solver pbsa` |
| Protein force field | `ff14SB` | `--force-field ff19SB`, `ff99SBildn` |
| Protonation | `propka` at pH `7.4` | `--protonation fixed`, `--ph <value>` |
| Export formats | `amber` (plus PDB & OpenMM XML/JSON) | `--formats amber openmm gromacs` |

```sh
# Non-interactive run with explicit settings and all three export formats
cationize-protein run --pdb prepared/protein.pdb --output series \
  --engine openmm --solver coulomb --force-field ff14SB \
  --formats amber openmm gromacs
```

Re-running the exact same command resumes an interrupted calculation after
verifying all saved checkpoints. If a remaining carboxyl site is sterically
buried at the selected probe radius, a warning is logged and the site is
ranked via the buried-site fallback.

## Find your files

```text
series/
├── residue_templates/    Individual residue libraries: ASC.lib, GLC.lib, …
├── 0/ … 100/
│   ├── protein.pdb       Minimized structure with modified residue names
│   ├── system.xml        OpenMM physical system
│   ├── state.json        Coordinates and chemical state
│   └── exports/          parm7/rst7, top/gro/itp, protein.lib, parameters.frcmod, PARAMETERS.json
└── internal/             Intermediate steps and recovery records
```

<details>
<summary><b>File lookup table & ligand complexation (Amber / GROMACS)</b></summary>

| What you need | Where to find it |
|---|---|
| Structure for inspection or docking | `25/protein.pdb` |
| Ready Amber protein | `25/exports/system.parm7` and `system.rst7` |
| Ready OpenMM protein | `25/exports/system.xml` and `state.json` |
| Ready GROMACS protein | `25/exports/system.top`, `system.gro`, `forcefield.itp`, `protein.itp` |
| Single modified residue template | `residue_templates/ASC.lib`, `GLC.lib`, etc. (Amber OFF format) |
| Complete Amber protein unit & parameters | `25/exports/protein.lib`, `parameters.frcmod`, `PARAMETERS.json` |

See [output and parameter reuse](docs/OUTPUTS.md) for details.

### Add a ligand

Start with the cationized `protein.pdb` and an independently parameterized
ligand in the same coordinate frame:

- **Amber:** copy one checkpoint's `exports/` folder, keep the atom-type
  declarations from `rebuild.in`, load `parameters.frcmod` + `protein.lib` and
  the ligand parameters, then `combine` in LEaP.
- **GROMACS:** follow the [complex tutorial](docs/COMPLEXES.md) to merge
  protein and ligand topologies into a consistent `TOP`/`GRO`/`ITP` set.

</details>

## Examples

### HSA example

[examples/HSA](examples/HSA/README.md) contains the **585-residue HSA series**
at 0%, 25%, 50%, 75% and 100% (17 disulfides preserved), with PDBs, Amber
topologies/coordinates and modification tables. For a fast test on a 3-residue
peptide:

```sh
cationize-protein run --pdb examples/pdb_checks/reference_peptide.pdb \
  --output peptide_series --engine openmm --solver coulomb \
  --force-field ff14SB --formats amber openmm gromacs
```

### Lysozyme walkthrough

Hen egg-white lysozyme (`1AKI`, UniProt `P00698`; 129 mature residues, 4
disulfides, 1 His, 10 carboxyl groups) exercises retrieval, optional MODELLER
reconstruction, PDBFixer repair, UniProt disulfide mapping, PROPKA protonation
and progressive cationization (including the buried-site fallback on `Asp66`).

<details>
<summary><b>Show complete Lysozyme commands (1AKI / P00698)</b></summary>

```sh
# 1. Download deposited PDB and automatically retrieve its mapped UniProt record (P00698)
cationize-protein fetch --pdb 1AKI --output lysozyme_fetch

# 2. Keep protein ATOM/TER/END records (or use an incomplete template)
grep -E '^(ATOM  |TER|END)' lysozyme_fetch/1AKI.pdb > lysozyme_protein.pdb

# 3. Optional: if whole residues/loops are missing, reconstruct with MODELLER
cationize-protein model --method modeller --pdb lysozyme_protein.pdb \
  --fasta lysozyme_fetch/mature.fasta --output lysozyme_modeled

# 4. Repair missing heavy atoms (PDBFixer), map UniProt disulfides and assign His (PROPKA)
cationize-protein prepare lysozyme_protein.pdb \
  --uniprot-record lysozyme_fetch/uniprot.json \
  --output lysozyme_prepared --ph 7.4

# 5. Inspect the prepared structure
cationize-protein check lysozyme_prepared/protein.pdb

# 6. Run 0%, 25%, 50%, 75% and 100% cationization and export Amber, OpenMM and GROMACS
cationize-protein run --pdb lysozyme_prepared/protein.pdb \
  --output lysozyme_series --engine openmm --solver coulomb \
  --force-field ff14SB --formats amber openmm gromacs
```

If whole residues were rebuilt in step 3, pass `lysozyme_modeled/protein.pdb`
to `prepare` in step 4.

</details>

## Options

Use `cationize-protein --help` or `cationize-protein COMMAND --help` at any
time. Further documentation: [installation](docs/INSTALL.md),
[methods](docs/METHODS.md), [force fields](docs/FORCE_FIELDS.md) and
[validation](docs/VALIDATION.md).

<details>
<summary><b>Cationization flags (<code>cationize-protein run</code>)</b></summary>

| Flag | Purpose, accepted values and requirements | Default |
|---|---|---|
| `--pdb` | Native protein PDB; use this or `--topology`, not both. | Required unless `--topology` is supplied |
| `--topology` | Start from an existing protein-only topology instead of a PDB: Amber (`.parm7`/`.prmtop`), GROMACS (`.top`), CHARMM (`.psf`) or OpenMM (`.xml`) (requires `--coordinates`). | Not used |
| `--coordinates` | Coordinates matching `--topology`: Amber (`.rst7`/`.inpcrd`), GROMACS (`.gro`), CHARMM (`.crd`/`.pdb`) or OpenMM (`state.json`/`.pdb`). | Required with `--topology` |
| `--output` | New output directory, or the same verified run to resume. | Required |
| `--percentages` | Saved cumulative levels (0–100); baseline 0% is always saved. | `25 50 75 100` |
| `--engine` | Minimizer: `openmm`, `amber`, `sander` or `pmemd`. | `openmm`* |
| `--engine-executable` | Path to serial `sander` or `pmemd`. | Search `AMBERHOME` and `PATH` |
| `--solver` | Potential: `coulomb`, `apbs` or `pbsa`. | `coulomb`* |
| `--apbs` / `--pbsa` | Path to installed `apbs` or `pbsa` executable. | Private install, `PATH`, then saved path |
| `--force-field` | Protein parameters: `ff14SB`, `ff19SB` or `ff99SBildn`. | `ff14SB`* |
| `--formats` | Parameter exports: `amber`, `openmm`, `gromacs`. | `amber` |
| `--protonation` | `propka` (re-check native protonation with PROPKA at each step) or `fixed` (keep all residue protonation states frozen as supplied in the input). | `propka` |
| `--ph` | Target pH (0–14); use the same value in `prepare` and `run`. | `7.4` |
| `--probe-radius` | Positive probe Lennard-Jones half-Rmin in Å. | `1.369` Å |
| `--parameters` | Custom validated residue library directory. | Bundled library |

*Omitted engine, solver and force field are prompted in interactive terminals.

</details>

<details>
<summary><b>Preparation, retrieval and modeling flags (<code>prepare</code>, <code>fetch</code>, <code>model</code>)</b></summary>

| Flag | Command & purpose | Default |
|---|---|---|
| `--uniprot-record` | `prepare`: local `uniprot.json` from `fetch` to map and verify intrachain disulfides. | Not used |
| `--protonation` | `prepare`: `propka` assigns histidine states (`HID`/`HIE`/`HIP`) at `--ph`; `fixed` keeps the residue names already present in the input PDB without running PROPKA. | `propka` |
| `--ph` | `prepare`: pH for initial histidine assignment (0–14). | `7.4` |
| `--pdb` | `fetch`: download deposited PDB coordinates (and auto-fetch unique `DBREF` UniProt record). `model`: incomplete native template for MODELLER. | Optional for `fetch`; required for `modeller` |
| `--uniprot` | `fetch`: download UniProt sequence, annotations, disulfides and linked PDB list (`pdb_candidates.csv`). | Optional |
| `--method` | `model`: `modeller` or `esmfold`. | Required |
| `--fasta` | `model`: complete native sequence of one chain. | Required |
| `--python` | `model`: path to the external Python interpreter where MODELLER or ESMFold2 is installed (when installed in a separate venv/Conda environment). | Env variable, then active Python / `PATH` |
| `--esmfold-model` | `model`: ESMFold2 checkpoint (`fast` or `full`). | `fast` |
| `--model-revision` | `model`: Hugging Face commit/branch revision. | `main` |
| `--device` | `model`: `auto`, `cpu` or `cuda`. | `auto` |
| `--output` | New output directory for `prepare`, `fetch` or `model`. | Required |

</details>

<details>
<summary><b>Advanced physical, grid and minimization flags (<code>run</code>)</b></summary>

| Flag | Purpose | Default |
|---|---|---|
| `--platform` | OpenMM platform: `CPU`, `Reference`, `CUDA` or `OpenCL`. | `CPU` |
| `--temperature` | Temperature of the electrostatic calculation (K). | `298` K |
| `--solute-dielectric` | Protein interior dielectric for PB. | `2` |
| `--solvent-dielectric` | Solvent dielectric for electrostatics. | `78.5` |
| `--salt` | Monovalent salt concentration for potential ranking (M). | `0.15` M |
| `--pb-ion-radius` | Electrolyte ion radius for PB (Å). | `2` Å |
| `--pb-spacing` | Maximum fine-grid spacing for PB (Å). | `0.75` Å |
| `--sampling-spacing` / `--grid-spacing` | Spacing of the local candidate probe lattice (Å). | `1` Å |
| `--sphere-radius` | Radius of candidate sampling spheres around carboxyl oxygens (Å). | `5` Å |
| `--pbsa-nonlinear` | Use nonlinear PBSA (requires `--solver pbsa`). | Off (linear ICCG) |
| `--keep-maps` | Keep potential DX maps produced at each APBS step. | Off |
| `--beam-width` | Number of steric-search candidates retained. | `8` |
| `--max-iterations` | Maximum iterations per minimization. | `20000` |
| `--restraint` | Native backbone/disulfide positional restraint strength, kcal/(mol Å²). | `10` kcal/(mol Å²) |
| `--tolerance` | Minimization RMS force target, kcal/(mol Å). | `0.1` kcal/(mol Å) |
| `--verbose` | Print additional diagnostic messages. | Off |

</details>

## Software and references

| Software | Role | Installation & connection | Primary reference |
|---|---|---|---|
| [OpenMM](https://openmm.org/) | Minimization, molecular systems and hydrogen placement | Included in the Python runtime | Eastman et al., *J. Phys. Chem. B* 2024 ([10.1021/acs.jpcb.3c06662](https://doi.org/10.1021/acs.jpcb.3c06662)) |
| [PDBFixer](https://github.com/openmm/pdbfixer) | Default repair of missing atoms in native residues | Included | Eastman et al., *PLOS Comput. Biol.* 2017 ([10.1371/journal.pcbi.1005659](https://doi.org/10.1371/journal.pcbi.1005659)) |
| [PROPKA](https://github.com/jensengroup/propka) | Initial histidine states and protonation checks after each modification | Included | Olsson et al., *J. Chem. Theory Comput.* 2011 ([10.1021/ct100578z](https://doi.org/10.1021/ct100578z)); Søndergaard et al., *J. Chem. Theory Comput.* 2011 ([10.1021/ct200133y](https://doi.org/10.1021/ct200133y)) |
| [ParmEd](https://parmed.github.io/ParmEd/html/index.html) | Topology handling and export | Included | Shirts et al., *J. Comput.-Aided Mol. Des.* 2017 ([10.1007/s10822-016-9977-1](https://doi.org/10.1007/s10822-016-9977-1)) |
| [Amber / AmberTools](https://ambermd.org/) | Optional Amber minimization and PBSA electrostatics | Separate — [setup guide](docs/OPTIONAL_SOFTWARE.md#ambertools-and-amber) | Case et al., *J. Chem. Inf. Model.* 2023 ([10.1021/acs.jcim.3c01153](https://doi.org/10.1021/acs.jcim.3c01153)) |
| [GROMACS](https://www.gromacs.org/) | Run exported GROMACS systems | Separate — [setup guide](docs/OPTIONAL_SOFTWARE.md#gromacs) | Abraham et al., *SoftwareX* 2015 ([10.1016/j.softx.2015.06.001](https://doi.org/10.1016/j.softx.2015.06.001)) |
| [APBS](https://www.poissonboltzmann.org/) | Poisson-Boltzmann electrostatics | Separate — [setup guide](docs/OPTIONAL_SOFTWARE.md#apbs) | Jurrus et al., *Protein Sci.* 2018 ([10.1002/pro.3280](https://doi.org/10.1002/pro.3280)) |
| [MODELLER](https://salilab.org/modeller/) | Optional template-based native modeling | Separate (licence required) — [setup guide](docs/OPTIONAL_SOFTWARE.md#modeller) | Webb & Sali, *Curr. Protoc. Bioinform.* 2016 ([10.1002/cpbi.3](https://doi.org/10.1002/cpbi.3)) |
| [ESMFold2 / ESMFold2-Fast](https://huggingface.co/biohub/ESMFold2) | Optional prediction of native structures from sequence | Separate ([PyTorch](https://pytorch.org/) + [Transformers](https://huggingface.co/docs/transformers/main/model_doc/esmfold2)) — [setup guide](docs/OPTIONAL_SOFTWARE.md#esmfold2) | Lin et al., *Science* 2023 ([10.1126/science.ade2574](https://doi.org/10.1126/science.ade2574)) |

<details>
<summary><b>Force fields and physical model references</b></summary>

- **Amber ff14SB**: Maier et al., *J. Chem. Theory Comput.* 2015 ([10.1021/acs.jctc.5b00255](https://doi.org/10.1021/acs.jctc.5b00255)).
- **Amber ff19SB**: Tian et al., *J. Chem. Theory Comput.* 2020 ([10.1021/acs.jctc.9b00591](https://doi.org/10.1021/acs.jctc.9b00591)).
- **Amber ff99SB-ILDN**: Lindorff-Larsen et al., *Proteins* 2010 ([10.1002/prot.22711](https://doi.org/10.1002/prot.22711)).
- **GBn2 implicit solvent**: Nguyen, Roe & Simmerling, *J. Chem. Theory Comput.* 2013 ([10.1021/ct3010485](https://doi.org/10.1021/ct3010485)).
- **RESP charge fitting**: Bayly et al., *J. Phys. Chem.* 1993 ([10.1021/j100142a004](https://doi.org/10.1021/j100142a004)).

</details>

## Residue parameters

The bundled library contains **32 chemical states** covering modified Asp/Glu
side chains and C termini across terminal, histidine and cysteine contexts
([validation summary](parameters/validation_summary.csv), [parameter guide](docs/PARAMETERS.md)).

<details>
<summary><b>Parameterization protocol (QM / RESP / force-field transfer)</b></summary>

Four accepted HSA templates (`ASC`, `GLC`, `ASX`, `LET`) retain their original
fitted charges; the other 28 states were developed with capped amino-acid models:

| Component | How it was obtained |
|---|---|
| Quantum model | ACE/NME caps reproduce the peptide context; disulfide variants include an explicit partner |
| Geometry and electrostatic potential | Two constrained backbone conformations optimized at HF/6-31G(d), with Merz-Kollman ESP samples |
| Partial charges | Joint two-stage RESP fitting, with fixed cap charges and the specified integer residue charge |
| Bonded and Lennard-Jones terms | Amber protein atom types and transfer by chemical analogy |

The introduced ammonium remains +1. The `ff14SB`, `ff19SB` and `ff99SBildn`
profiles share the custom RESP charges and use their corresponding native
protein contexts. To extend the library, follow [docs/PARAMETERS.md](docs/PARAMETERS.md)
and pass the custom directory with `--parameters`.

</details>

## Acknowledgements

This work was developed within the **NanoBio Interface Lab** at the Department
of Chemistry "Giacomo Ciamician," Alma Mater Studiorum – University of Bologna,
and the Preclinical & Translational Research in Oncology lab (PRO), IRCCS
Azienda Ospedaliero–Universitaria di Bologna.

## License

This software, its bundled residue parameter libraries, example structures and
documentation are released under the
[Apache License, Version 2.0 (Apache-2.0)](LICENSE).

See [CITATION.cff](CITATION.cff) for citation metadata. Third-party software and
reference force fields retain their original terms and attribution.
