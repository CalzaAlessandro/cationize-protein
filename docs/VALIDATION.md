# Implementation validation — 20 October 2026

Tests use short peptides, synthetic systems and the initial 585-residue HSA model. No HSA production was submitted
or restarted. Original structures and accepted fitted RESP libraries were
preserved. A machine-readable summary is in [validation/latest.json](../validation/latest.json).

## Current standalone package

A fresh private Python 3.12 installation on macOS installed the complete
runtime, including PDBFixer 1.12.0, without Amber, GROMACS or APBS executables.
The installation completed and reported missing external programs. APBS-based
calculations still require an installed solver; no alternative model is
silently selected. Conda and venv recipes install the same project, but creation
of a new Conda environment was not executed in this validation.

The installed package passed **242 tests, with no skips**, including initial
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
single-structure estimate. The later minimized workflow audit remains separate.
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

Neither a licensed MODELLER calculation nor a real ESMFold2 prediction was
executed. These modelers and their numerical results remain unqualified in
this package. No MODELLER key, folding environment or model weights were
installed. Automatic local grafting of predicted loops is not implemented;
the optional command generates an explicitly identified full native candidate.
