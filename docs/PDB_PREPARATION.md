# Initial PDB preparation

`check` inspects the original text before a structure reader can discard
ambiguity. `prepare` uses [PDBFixer](https://github.com/openmm/pdbfixer) to add
missing heavy atoms to existing native protein residues. This step is separate
from cationization and never modifies the source file.

```sh
cationize-protein check input.pdb
cationize-protein prepare input.pdb --output prepared
cationize-protein check prepared/protein.pdb
cationize-protein run --pdb prepared/protein.pdb --output series
```

With the private installer, replace `cationize-protein` with
`.venv/bin/cationize-protein`. `prepare --help` describes its two arguments.

## Files

```text
prepared/
  original.pdb       Exact copy of the input
  protein.pdb        Input records plus newly built heavy atoms
  preparation.log   Atom additions, checks and remaining preparation decisions
```

The output directory must be new. Both successful preparation and failure
reasons are printed on screen. Persistent output is created only after the
repaired candidate passes the PDB checks; no JSON file is required.

## What is repaired

Missing side-chain atoms, isolated backbone atoms in an existing residue and
terminal OXT atoms can be built. PDBFixer generates candidate coordinates with
a fixed random seed. Only the new heavy-atom records are inserted into the
original text. Existing atom coordinates, names, residue numbering, chain
identifiers, protonation labels and SSBOND/CONECT records remain unchanged.
Backbone or side-chain reconstruction still requires the subsequent force-field
minimization and structural review. Preparation is not an experimental
reconstruction of an unresolved conformation.

PDBFixer is used only on native input. Cationized residues and other nonstandard
groups are rejected rather than converted to standard amino acids. No CCD
template is downloaded during preparation.

## What stops preparation

| Finding | Required action |
|---|---|
| Missing residues identified from SEQRES | Supply a separately reviewed sequence/loop model |
| A peptide gap or incompatible disulfide declaration | Correct or explicitly resolve the structural model |
| Multiple models or alternate locations | Select one coherent model/conformer before preparation |
| Ligands, solvent, ions, caps or nonstandard residues | Prepare a separate protein-only input; preserve the source |
| Truncated records or non-finite coordinates | Recover a correctly formatted coordinate file |
| Existing output directory | Choose a new directory |

If SEQRES is absent, a coordinate file cannot establish completeness against
the intended sequence. Numbering gaps alone do not prove missing residues.
Large coordinate gaps are flagged by the independent peptide geometry check.

## Protonation and disulfides

This command does not add or remove hydrogens and does not use PDBFixer's
default pH assignment. The existing cationization protocol adds hydrogens,
minimizes and runs PROPKA at the target pH, including after each chemical step.
Specify HID/HIE/HIP or supply diagnostic histidine ring hydrogens before `run`;
PROPKA does not assign a neutral histidine tautomer.

Declared disulfide records and CYX labels are preserved. The PDB checker
verifies declared sulfur identities and geometry; the subsequent workflow
verifies that these bonds are present in the prepared molecular graph.

## Reproducible broken-input examples

```sh
cationize-protein prepare examples/pdb_checks/missing_backbone_atom.pdb \
  --output repaired_peptide
cationize-protein prepare examples/pdb_checks/wrong_columns.pdb \
  --output invalid_example
```

The first command adds one CA to the synthetic ASP–ASP–LEU peptide; the second
fails before creating `invalid_example`. See the example folder for the
unmodified source files and diagnostic output.

PDBFixer's broader capabilities are described in its
[official manual](https://github.com/openmm/pdbfixer/blob/master/Manual.html).
This interface selects only heavy-atom preparation.
