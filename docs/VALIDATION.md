# Implementation validation — 21 October 2026

Tests use short peptides, synthetic systems and the initial 585-residue HSA model. No HSA production was submitted
or restarted. Original structures and accepted fitted RESP libraries were
preserved. A machine-readable summary is in [validation/latest.json](../validation/latest.json).

## Current standalone package

A fresh private Python 3.12 installation on macOS installed the complete
runtime, including PDBFixer 1.12.0, without Amber, GROMACS or APBS executables.
The installation completed and reported missing external programs. Fresh
installations select OpenMM, screened Coulomb and ff14SB. APBS-based
calculations still require an installed solver; a requested solver is never
replaced automatically. Conda and venv recipes install the same project, but creation
of a new Conda environment was not executed in this validation.

The installed package passed **268 tests, with no skips**, including initial
PDB repair, PROPKA histidine assignment, residue-library export, retrieval and
optional modeling interfaces. The earlier Linux run passed 202 tests before PDBFixer
preparation was added; the new preparation tests have currently run on macOS.
Source files and installed Python modules were compared byte for byte.

PDB preparation was tested on actual structures with missing side-chain,
backbone and terminal atoms. It retained supplied heavy coordinates and declared disulfides.
Histidine labels and ring H are updated when initial PROPKA is selected. The repaired peptide passed parameter assignment.
Malformed records, alternate locations, missing SEQRES residues, cationized
residues and existing output directories were rejected. A separate console
example demonstrates one missing CA repaired and a truncated record rejected
without an output directory.

A real ASP–ASP–LEU demonstration completed cumulative 0/25/50/75/100% checkpoints
with OpenMM minimization, PROPKA and all three export formats. A second complete
run started from the PDBFixer-repaired peptide. Both demonstrations selected
screened-Coulomb ranking explicitly; they do not test a new APBS installation.
Three available carboxyl groups produce rounded counts of 0, 1, 2, 2 and 3.

A third complete peptide run used only `--pdb` and `--output`. The fresh-install
defaults completed all five checkpoints with OpenMM, screened Coulomb and
PROPKA, without an APBS executable. Checkpoint hashes and exported artifacts
were verified. Existing APBS/PBSA run records retained their solver choices,
implementation fingerprints and unchanged `RUN.json` bytes during full resume
validation. No APBS/PBSA calculation was performed for that compatibility check.

## Current PDB and split GROMACS exports

PDB atom identity/order and coordinates, moved relative ITP includes, the
`CP_` atom-type namespace, disconnected chains, CMAP retention and composition
with an independently defined synthetic ligand are covered by regression
tests. No ligand-combination utility is provided.

Three previously generated peptide systems were exported again using the
current split ITP writer. ParmEd reloaded these files and OpenMM compared the
physical graph, charges, LJ, exclusions, 1–4 terms, CMAP, energy and forces at
identical original coordinates. All checks passed:

| Profile | GROMACS round-trip energy error, kJ/mol | CMAP terms |
|---|---:|---:|
| ff14SB | 1.139 × 10⁻⁵ | 0 |
| ff19SB | 1.132 × 10⁻⁵ | 1 |
| ff99SBildn | 2.689 × 10⁻⁶ | 0 |

These are serialization and physical-model round trips, not new GROMACS binary
executions. The current split include layout has not yet been rerun in an
external GROMACS executable.

The [complex-assembly tutorial](COMPLEXES.md) was checked with an ff14SB
peptide, an ff19SB peptide and the 25% HSA example, each combined with an
independently parameterized synthetic ligand. The documented ParmEd script
retained atom order, bond graphs, all 17 HSA disulfides and the ff19SB CMAP.
Coordinate round-trip errors were at most 5.1 × 10⁻⁶ Å. Two-coordinate
energy checks for the peptide complexes differed by at most
1.96 × 10⁻⁵ kJ/mol. LEaP and `gmx grompp` examples are documented but were
not executed locally; ligand parameter fitting and production MD are outside
these checks.

## Earlier independent engine qualification

Before the standalone rename and ITP split, Amber 24/AmberTools 25 and GROMACS
2025.4 evaluated the same physical peptide representations. These were vacuum
CPU single points, not MD. Coordinates were identical; GROMACS coordinates used
nine decimal places to remove ordinary GRO rounding.

| Case | Absolute Amber difference, kJ/mol | Absolute GROMACS difference, kJ/mol |
|---|---:|---:|
| ff14_apbs | 0.029605 | 0.00036694 |
| ff19SB_coulomb | 0.04594 | 0.00053912 |
| ff99SBildn_coulomb | 0.041589 | 8.3444 × 10⁻⁵ |

The ff19SB XML/JSON reconstruction retained its 24 × 24 CMAP and passed energy
and force comparison. A conventional frcmod alone does not encode CMAP.
Complete exports include actual 1–4 scaling derived from physical exceptions;
these can differ from bare ParmEd torsion defaults.

Earlier serial integration also covered APBS/PROPKA, Amber minimization,
linear PBSA, native PDB input and verified resume. The capped-template matrix
had 64 passing assignments for ff14SB and 62 each for ff19SB and ff99SBildn,
including vacuum/GBn2 checks. CMT was explicitly unsupported under the latter
two profiles because their distributed native CCYM template is absent. All
32 chemical states were supported under ff14SB. Accepted custom RESP charges
were not refitted or altered for this repository work.

## Qualification limits

These checks establish tested behavior and parameter conversion. They do not
establish transferability to arbitrary proteins or replace experimental
validation. No additional QM torsion fitting was performed for alternate
profiles. Nonlinear PBSA remains unqualified; APBS is the earlier tested
nonlinear reference. Old engines without CMAP support cannot run ff19SB
correctly. Other engine releases require separate checks.

The optional modeling changes received an independent source review. A reported
invalid-sulfur error path was corrected and covered by a regression test.
Earlier export changes received author review and the numerical checks above;
an independent review of the whole package is not claimed.
Private environments, calculation outputs and development logs stay outside
the source repository.

## Initial histidines and the HSA example

`prepare` was run on the full native HSA starting PDB at pH 7.4. PROPKA assigned
14 HIE and two HIP states. The predicted pKas were 8.12 for His67 and 7.42 for
His247. His247 is close to the selected pH and is flagged as borderline; its
binary assignment represents the slightly more populated state from this
single-structure estimate. The workflow checks protonation again after minimization.
The existing minimized campaign structures were not changed.

The 585-residue HSA example includes the five existing minimized PDB/topology/
restart sets. All source hashes, atom identity/order and coordinate agreement
were checked. The modification sets are cumulative, and all 17 disulfides are
retained. Modified partial charges match the accepted libraries within
2.7 × 10⁻⁹ e. The preview uses these actual structures; its five frames show
chemical checkpoints, not MD motion. VMD used STRIDE and ignored .vmdrc.

A new serial peptide run and resume passed with all three export formats and
a root residue_templates directory containing 32 hash-verified libraries.
An activated venv resolves the short cationize-protein command.

## Retrieval and optional modelers

The actual `fetch` command downloaded public RCSB 1N5U and UniProt P02768
records. It preserved the deposited PDB and canonical 609-residue sequence,
extracted the uniquely annotated mature chain at positions 25–609 (585
residues), and listed 189 linked PDB candidates. This macOS Python had no
configured CA file; the network check used the trusted system bundle via
`SSL_CERT_FILE=/etc/ssl/cert.pem`. TLS verification remained enabled.

PDBFixer remains the default preparation. MODELLER and ESMFold2 are discovered
without importing their engines, testing a licence or downloading weights.
Actual CLI tests confirmed that both missing installations produce a clear
error and create no model output. Interface tests cover separate interpreters,
bad probe output, FASTA validation, template mapping, DOPE candidate selection,
the ESMFold2 API contract, sequence/geometry rejection and failed-run diagnostics.
The standalone worker tests use external API fixtures.

MODELLER 10.8 from the existing Homebrew installation completed two real
native-fragment tests: one deleted residue was rebuilt in HSA fragments of
30 and 16 residues, with five scored models per run. The latter preserved a
disulfide and passed PDBFixer preparation and subsequent ff14SB assignment
(16 residues, 217 atoms, one SG–SG topology bond). The actual test identified
and corrected a truncated SSBOND serialization error. A regression test now
checks that written disulfide records pass the preparation checker.

Optional modeling software can be discovered in another Python on PATH,
including Homebrew's Python 3.15, without changing the main Python 3.12
runtime. Explicit interpreter overrides remain authoritative.

The isolated cluster environment contains PyTorch 2.11.0+cu128, Transformers
5.20.1 and Accelerate 1.15.0, with the public ESMFold2-Fast checkpoint revision
`45fe8656f5b3ef493c17fcf9abe9a2968902e712`. A first CUDA run failed because
bfloat16 model weights conflicted with float32 LayerNorm activations. The
adapter now uses float32, covered by a regression test. The second real CLI
run completed on an NVIDIA L40S (Slurm 1268030, 17 seconds) and produced a
30-residue native HSA-fragment candidate that passed the sequence/geometry
checks and included per-residue confidence.

These short-fragment runs establish interface operation, not prediction
accuracy for full HSA or arbitrary missing loops. The full ESMFold2 checkpoint,
MSA conditioning and mixed-precision inference remain untested. Automatic
local grafting of predicted loops is not implemented; the optional command
generates an explicitly identified full native candidate.

## UniProt disulfide assignment

The installed `prepare --uniprot-record` command was tested on the native
585-residue HSA model with all CYX labels changed to CYS and the supplied
SSBOND/CONECT records removed in a private test copy. Exact sequence matching
mapped the chain to canonical UniProt P02768 positions 25–609. All 17 annotated
pairs were recovered, including canonical Cys77–Cys86 mapped to mature
Cys53–Cys62. The output contained 34 CYX residues, retained every supplied heavy
coordinate and preserved the original input and reference bytes. Default
PROPKA assignment covered 16 histidines; subsequent ff14SB assignment passed.

Eighteen focused tests cover sequence ambiguity/mismatch, partial constructs,
insertion codes, homomers, matching and conflicting supplied SSBOND records,
sulfur geometry, thiol H removal, uncertain/qualified annotations, malformed
records, duplicate atom serials and noncontiguous chain blocks. Final graph
verification checks sulfur partners and actual topology-chain identity.
Independent review findings were reproduced and fixed before the installed
full suite was run.

The previously checked README command examples parse with the CLI, and 20 documented
`run` defaults were checked against argparse. The HSA preview uses
five cumulative minimized structures, below a compact
light/dark vector banner. A shared orthographic camera uses yaw 180° and
pitch 30°. A camera roll of −15.8436° aligns the projected native IIB Cα axis
(residues 292–383) vertically. The roll changes the camera Updir only;
scene coordinates, normals, materials and lights are preserved. The same
correction is used for the full-protein cartoons and potential surfaces.
The Structure-colored NewCartoon has opacity 0.20 in the light
version and 0.45 in the dark version. Dark-theme colors are brighter;
modified-residue sticks retain opacity 1.00 and atom-name colors in both
versions. AOChalky is retained, and
`.vmdrc` is bypassed.

Each theme has five downloadable PNGs at the genuine 3840 × 2160 render grid. Tachyon uses
24 antialiasing samples and the recorded lighting/transparency settings.
The 23 October README animation keeps a text-only, colored cycle on the left
and shows the molecular images successively on the right. There are no step
numbers or central reaction formula in the cycle. Each 1600 × 900 GIF has
24 frames and a 67 s loop. HSA structures and potential surfaces remain
specific to each percentage. Three fixed operation illustrations are reused
for amidation, minimization and PROPKA, with the heading "Method illustration"
and no percentage-specific result labels. The 0% reference has no amidation
page. The cycle uses equal calculated clearances at arrow endpoints and a
centered three-line "Amidation / + / rotation" label.

Five screened-Coulomb fields were calculated with the current function and
the respective checkpoint coordinates and topology charges, using dielectric
78.5, 0.15 M salt and 298 K. All fields share a 3 Å grid and a display range
of −5 to +5 kBT/e. The fields and molecular geometry are unchanged by the new
orange–gray–violet display scale. The colorbar is read directly from VMD's
actual 1024-color lookup table. The scale's neutral gray improves visibility
against white. Both themes use the same scale, and 65 scale colors are
reserved in the GIF palette.

The operation illustrations reuse one archived example: Asp249 modification
and local minimization from step 99, and an HIE imidazole fragment. These
images explain operations throughout the cycle. The minimization overlay
retains its rigid local fit and original displacement scale. Percentage-specific
residue identities and pKa values are omitted from these illustrative pages.

The percentage-specific molecular images were rendered at 3840 × 2160 with
VMD/Tachyon. Protein crops are downsampled without upsampling. GIF frames are
composited on the matching light/dark background with one palette per theme;
the companion static PNG retains alpha transparency. The underlying render
alpha comes from independent geometry-transmittance passes that retain the
actual representation opacities. White molecular features are preserved.
The final diagrams and sampled GIF frames were inspected on light and dark
backgrounds. Source structures and accepted parameters are unchanged. Text has
no outline, and the README selects the matching browser color scheme. Render settings,
alpha checks and file hashes are recorded in
[HSA_PREVIEW.json](images/HSA_PREVIEW.json).

The [optional software tutorial](OPTIONAL_SOFTWARE.md) documents installation
and discovery for APBS, AmberTools/Amber, GROMACS, MODELLER and ESMFold2.
Shell blocks passed Bash and Zsh syntax checks; Python snippets and CLI
options were checked separately. Installation sources were consulted, but
these programs were not reinstalled for this documentation update.
