# PDB check examples

These files are synthetic test inputs, not HSA structures. Run the commands
from the repository root after activating the environment:

```sh
cationize-protein check examples/pdb_checks/reference_peptide.pdb
cationize-protein check examples/pdb_checks/missing_backbone_atom.pdb
cationize-protein check examples/pdb_checks/wrong_columns.pdb
```

`reference_peptide.pdb` is the complete three-residue peptide used as the
reference for the examples. It passes. The other two inputs deliberately fail:

- `missing_backbone_atom.pdb` has the first CA removed. The checker identifies
  residue 1 and the missing atom, then requests a reviewed reconstruction.
- `wrong_columns.pdb` contains an ATOM record without the required fixed-column
  coordinate layout. The checker reports a truncated atom record.

The corresponding `.log` files contain the actual console output. Exit codes
are 0 for PASS and 1 for structural FAIL. The input files remain unchanged;
no JSON report or new scientific output is created unless requested.

To see that `run` also refuses the broken input:

```sh
cationize-protein run \
  --pdb examples/pdb_checks/missing_backbone_atom.pdb --output rejected_example \
  --engine openmm --solver coulomb --force-field ff14SB
```

It exits with an error and does not create `rejected_example`. The explicit
Coulomb model makes this input-validation example independent of APBS; it is
not a statement that Coulomb and Poisson-Boltzmann are equivalent.

An optional complete small-peptide demonstration is:

```sh
cationize-protein run \
  --pdb examples/pdb_checks/reference_peptide.pdb --output peptide_example \
  --engine openmm --solver coulomb --force-field ff14SB \
  --formats amber openmm gromacs
```

This runs serial minimization and checks protonation with PROPKA. The peptide has
three candidate carboxyl groups. Rounded checkpoint counts are 0, 1, 2, 2 and
3 at 0/25/50/75/100%; 50% and 75% therefore share a chemical state in this tiny
example. Larger proteins generally have finer percentage resolution.

## Repair with PDBFixer

```sh
cationize-protein prepare examples/pdb_checks/missing_backbone_atom.pdb \
  --output repaired_peptide
```

This adds the missing CA and preserves every supplied atom record. The checked
example result is in `prepared/`; `preparation_success.log` records console
output. The same command with `wrong_columns.pdb` fails before creating output;
see `preparation_failure.log`. Preparation adds heavy atoms only. Whole missing
residues and protonation choices require separate handling as described in
[the preparation guide](../../docs/PDB_PREPARATION.md).
