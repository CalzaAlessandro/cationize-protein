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
| Missing residues identified from SEQRES | Supply a reviewed sequence/loop model |
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
