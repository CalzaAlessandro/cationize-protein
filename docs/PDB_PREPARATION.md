# Prepare the initial PDB

Activate your environment, then run:

```sh
source .venv/bin/activate
cationize-protein prepare input.pdb --output prepared
```

For Conda, use `conda activate cationize-protein` instead. `prepare` keeps the
source unchanged and writes a new directory. To inspect a file without
preparing it, use `cationize-protein check input.pdb`.

## Files

| File | Contents |
|---|---|
| `original.pdb` | Exact copy of the input |
| `protein.pdb` | Repaired heavy atoms and assigned histidine states |
| `preparation.log` | Changes and remaining structural issues |
| `histidines.csv` | pKa, pH, state and neutral-tautomer choice for each histidine |
| `protonation/` | PROPKA input, mapping, logs and predictions |

The last two entries are created when histidines are present and PROPKA is
selected. Results and failure reasons are also printed on screen. An existing
output directory is never overwritten.

## Histidines

PROPKA chooses the initial charge state at pH 7.4 by default. A histidine is
assigned HIP when its predicted pKa exceeds the target pH; otherwise it is
neutral. Use `--ph` to select another pH:

```sh
cationize-protein prepare input.pdb --output prepared --ph 7.4
cationize-protein run --pdb prepared/protein.pdb --output series --ph 7.4
```

For neutral histidines, supplied HID/HIE labels or a diagnostic ring hydrogen
are retained. Otherwise OpenMM selects HID or HIE from local hydrogen-bond
geometry. PROPKA predicts ionization; it does not select a neutral tautomer.
The CSV flags pKa values within 0.5 of the chosen pH for review.

The program updates histidine residue labels and adds or removes the required
ring hydrogens. Existing heavy-atom coordinates and identifiers are retained.
Other native side-chain and terminal states remain as supplied during this
step. The cationization workflow audits supported native protonation changes
after minimization and after each chemical modification. Introduced ammonium
groups stay +1.

The initial histidine calculation uses a temporary ff14SB protein model to
establish the molecular graph and place hydrogens. PROPKA receives the heavy
structure, not force-field atomic charges. The force field selected for the
subsequent cationization run is assigned separately.

To retain all supplied states during preparation:

```sh
cationize-protein prepare input.pdb --output prepared --protonation fixed
```

In fixed mode, an unlabelled HIS still needs HID/HIE/HIP or diagnostic ring
hydrogens before `run`.

## Missing atoms and structural issues

PDBFixer builds missing side-chain atoms, isolated backbone atoms in an
existing residue and terminal OXT atoms. A fixed seed is used for this step.
Only added heavy-atom records are inserted into the input text. Subsequent
force-field minimization is needed before interpreting their conformations.

| Finding | Required action |
|---|---|
| Missing residues identified from SEQRES | Supply a reviewed sequence/loop model, optionally generated with `model` |
| A peptide gap or incompatible disulfide declaration | Correct or resolve the structural model |
| Multiple models or alternate locations | Select one coherent model or conformer |
| Ligands, solvent, ions, caps or nonstandard residues | Prepare a separate protein-only input |
| Truncated records or non-finite coordinates | Recover a correctly formatted coordinate file |

PDBFixer is used only on native input. Cationized residues are rejected rather
than replaced by standard amino acids. Missing sequence cannot be inferred
reliably from coordinate numbering alone. Without SEQRES or an external
sequence, completeness beyond the supplied residues remains undetermined.

Declared disulfides and CYX labels are retained. The checker verifies sulfur
identities and geometry; the workflow then verifies the prepared bond graph.

OpenMM, used by PDBFixer and the parameter-assignment step, can identify S–S
bonds from cysteine sulfur proximity. During parameter assignment, bonded
cysteines are represented as CYX and the declared SSBOND/CYX identities are
checked against actual SG–SG topology bonds. Proximity alone is not an
experimental assignment of oxidation state; review inferred pairs, especially
in a newly predicted model. `pdb4amber` can still be used externally for an
Amber preparation route, but is not required for the PDBFixer/OpenMM route.
See [OpenMM disulfide detection](https://docs.openmm.org/latest/api-python/generated/openmm.app.topology.Topology.html#openmm.app.topology.Topology.createDisulfideBonds).

A FASTA file identifies cysteine positions but does not specify disulfide
partners. UniProt annotations or an experimentally assigned structure can
supply known pairs, provided their residue numbering is mapped to the input
construct. A CYX label alone also does not identify its partner.

As a geometry check, native HSA was loaded after changing all CYX labels to
CYS and removing SSBOND/CONECT records from an in-memory copy. OpenMM recovered
the same 17 SG–SG pairs as the supplied native topology. This check does not
establish pairing in a distorted or newly predicted structure. Sulfur-bound
HG hydrogens can also prevent automatic recognition and must be consistent
with the intended chemical state.

## Disulfides from UniProt

```sh
cationize-protein fetch --uniprot P02768 --output hsa_reference
cationize-protein prepare input.pdb --uniprot-record hsa_reference/uniprot.json \
  --output prepared
```

The reference is an ordinary local UniProt JSON record. `prepare` makes no
network request. It maps each continuous PDB chain to a unique, exact
contiguous region of the canonical sequence. PDB residue numbers and insertion
codes are retained. A mature HSA chain maps to UniProt positions 25–609, so
the canonical Cys77–Cys86 annotation becomes Cys53–Cys62 in a PDB numbered
from the mature N terminus. No fixed offset is assumed by the program.

For every mapped pair, both cysteine SG atoms must be present after atom
repair, and their distance must be between 1.7 and 2.5 Å. The prepared copy
uses CYX labels and explicit SSBOND/CONECT records. Thiol HG hydrogens on
these cysteines, including associated connectivity records, are removed;
the log reports their number. Supplied heavy coordinates are unchanged.
OpenMM reads the resulting PDB, and each annotated sulfur must have exactly
its expected SG partner in the bond graph.

| Additional output | Contents |
|---|---|
| `uniprot.json` | Unchanged reference record, whose SHA256 is recorded in the log |
| `uniprot_mapping.csv` | Canonical UniProt position mapped to PDB chain, residue number and insertion code |
| `disulfides.csv` | Canonical and PDB pairs, SG distance, verification status and annotation evidence |
| `preparation.log` | Sequence ranges, mapped pairs, removed thiol H and graph verification |

Existing matching SSBOND records are retained. A different declared partner,
ambiguous sequence match, mutation, unsupported interchain/uncertain annotation,
or incompatible geometry stops preparation without creating an output
directory. Automatic assignment accepts unqualified annotations or the label
`Intrachain`; free-text qualifications such as alternative partners or
condition-dependent bonds require review. Atom serials must be globally
unique so CONECT records and hydrogen removal are unambiguous. Separate
chain blocks are never merged merely because they share a chain identifier.
A pair entirely outside the supplied construct is recorded as
`outside_construct`; a pair with only one partner present requires review.
Multiple identical chains can receive separate intrachain pairs. Heteromers,
gapped sequence alignments, engineered variants and interchain disulfides
need a separately reviewed assignment; this option does not infer them or
move cysteines together to repair a distorted bridge.

With no `--uniprot-record`, the existing preparation behavior is unchanged.
Database annotations establish the source of an assignment; they do not
determine the oxidation state under a user's experimental conditions.

Source: [UniProt HSA annotation](https://rest.uniprot.org/uniprotkb/P02768.json),
[PDB SSBOND format](https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#SSBOND).

## Optional complete models

Preserve any experimentally assigned disulfide pairs before modeling. The
UniProt option above can verify known pairs on the resulting candidate.

MODELLER and ESMFold2 are separately selected with `cationize-protein model`;
they do not replace the default PDBFixer preparation or run automatically when
it finds a gap. Their outputs are full structural candidates and must be
checked and passed to `prepare` before cationization. MODELLER requires its
licence; ESMFold2 needs a compatible separate installation and public model
weights. See [commands and validation status](MODELING.md).

## Broken-input examples

```sh
cationize-protein prepare examples/pdb_checks/missing_backbone_atom.pdb \
  --output repaired_peptide
cationize-protein prepare examples/pdb_checks/wrong_columns.pdb \
  --output invalid_example
```

The first command adds one missing CA to a synthetic ASP–ASP–LEU peptide. The
second reports a truncated coordinate record and creates no output directory.
See [the examples](../examples/pdb_checks/README.md) for the files and logs.

Methods: [PDBFixer](https://github.com/openmm/pdbfixer),
[PROPKA](https://propka.readthedocs.io/en/latest/quickstart.html) and
[OpenMM hydrogen placement](https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addHydrogens).
