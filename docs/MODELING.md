# Optional native structure modeling

PDBFixer is the default for `prepare`. MODELLER and ESMFold2 are optional
components, discovered again when `doctor` or `model` is requested. A missing
modeler never prevents installation or ordinary PDBFixer preparation.

## Install the optional software separately

Follow the [optional software installation tutorial](OPTIONAL_SOFTWARE.md#modeller)
for MODELLER, or its [ESMFold2 recipe](OPTIONAL_SOFTWARE.md#esmfold2) for a
separate Python 3.12 environment with CPU/CUDA choices. Both Conda and
non-Conda routes are provided. MODELLER requires an individual licence from
the official distribution; the unrelated PyPI package named `modeller` does
not provide the structural modeling software. This package does not
distribute or request a MODELLER key.

Separate environments are supported, with or without Conda. Point the program
to the Python that can import the modeler:

```sh
export CATIONIZE_MODELLER_PYTHON=/path/to/modeller-environment/bin/python
export CATIONIZE_ESMFOLD_PYTHON=/path/to/esmfold-environment/bin/python
cationize-protein doctor
```

These are paths to Python interpreters, not installation directories. If the
variable is unset, detection checks the active cationize-protein Python, then
other Python interpreters found on `PATH`. This supports separate Homebrew
Python installations without changing the main environment. An
explicit invalid path produces a missing-interpreter message; it never falls
back to another environment. For a single modeling command, `--python` can
select the interpreter without setting an environment variable.

`doctor` checks files belonging to the installed Python packages. It does not
import MODELLER, validate its licence, initialize CUDA or download model
weights. Detection of code therefore does not certify a successful model run.
Licence, API compatibility and hardware errors are reported when modeling
actually starts.

## Choose the complete sequence

Provide one FASTA record containing the 20 canonical amino acids. Use the
sequence of the intended mature protein or construct. Include a signal peptide
only when it belongs to the intended construct. This interface accepts
one native chain; complexes, noncanonical chemistry and modified residues
must be prepared separately.

### MODELLER

```sh
cationize-protein model --method modeller --pdb incomplete.pdb \
  --fasta mature.fasta --output modeled
```

The supplied PDB provides the structural template. Its observed native
sequence must have one unambiguous mapping onto the complete FASTA, allowing
missing residues but no substitutions. Repeated sequences that permit
multiple gap placements require a separately reviewed alignment or model.
The adapter stops when the residue mapping is ambiguous.

The adapter normalizes a copy of the template, writes its alignment, builds
five comparative models and selects the lowest finite DOPE score among
successful models. Explicit template disulfide declarations and paired CYX
states are mapped to target positions and retained. The result is a complete
comparative model; supplied coordinates are not frozen. This is not a local
loop-only refinement protocol. The MODELLER seed is −8123 and is retained in
the provenance record.

### ESMFold2

```sh
cationize-protein model --method esmfold --fasta mature.fasta --output predicted
```

The default checkpoint is `biohub/ESMFold2-Fast`. To use the full checkpoint:

```sh
cationize-protein model --method esmfold --fasta mature.fasta \
  --esmfold-model full --device cuda --output predicted_full
```

Both options use single-sequence prediction in this adapter; MSA input is not
yet exposed. The result is the highest-pTM diffusion sample selected by the
upstream PDB writer. The seed is 0; checkpoint revision, resolved commit when
available, device and software versions are recorded. Set `--model-revision`
to reuse a particular published checkpoint revision.

`--device auto` selects CUDA if available, otherwise CPU. GPU execution is
recommended for practical use; CPU time and memory must be assessed for the
actual protein. The adapter uses float32. A real CUDA test of Transformers
5.20.1 found a LayerNorm dtype mismatch when all weights were loaded as
bfloat16; mixed precision is therefore not enabled in this profile. A requested
but unavailable CUDA device is an error.

Public weights are downloaded or loaded from the Hugging Face cache only
after an explicit `model --method esmfold` request. Neither the installer,
`doctor`, `prepare` nor cationization downloads them. No Biohub inference API
or API key is used, and the sequence is not submitted to a remote folding
service.

## Outputs and next step

The original PDB, when supplied, and FASTA are preserved. A successful model
contains `protein.pdb`, `MODEL.json`, `modeling.log` and the adapter's source
records. MODELLER also saves candidate scores and a source-to-target residue
map. ESMFold2 saves `confidence.csv`; pLDDT follows its documented 0-1 scale.
Model residue numbering can differ from the source PDB.

The result is checked for exact target sequence, heavy-atom completeness,
finite coordinates, one continuous chain, peptide C-N geometry and mapped
template disulfides. These checks do not establish experimental accuracy,
complete stereochemical quality or force-field readiness. A failed calculation
keeps `MODELING_FAILED.json` and diagnostic records, but no successful
`MODEL.json` or prepared `protein.pdb`. Existing output is never overwritten.

Inspect the candidate, then prepare it with the default PDBFixer/PROPKA step:

```sh
cationize-protein prepare modeled/protein.pdb --output prepared --ph 7.4
cationize-protein run --pdb prepared/protein.pdb --output series --ph 7.4
```

ESMFold2 predictions are not grafted automatically into the experimental PDB.
That would require separate alignment, junction optimization and validation.
Modeling is performed before cationization, never after each chemical step.

## Validation status

Discovery, missing-dependency errors, input validation, process adapters,
candidate selection and output checks are covered by automated interface tests.
External API fixtures are used for adapter regression tests. Actual MODELLER
10.8 runs reconstructed one removed residue in native HSA fragments of 30 and
16 residues. Each run produced five scored models. The 16-residue example
retained its disulfide and passed subsequent PDBFixer preparation and ff14SB
parameter assignment.

An isolated Linux Python 3.12 environment installed PyTorch 2.11.0+cu128,
Transformers 5.20.1 and Accelerate 1.15.0, and downloaded the public
ESMFold2-Fast checkpoint at revision
`45fe8656f5b3ef493c17fcf9abe9a2968902e712`. A real CUDA run on an NVIDIA L40S
completed the full CLI workflow for a 30-residue HSA fragment in float32;
the Slurm allocation lasted 17 seconds. The output passed sequence, heavy-atom,
coordinate and peptide-geometry checks and included per-residue confidence.

These are functional tests on short fragments. Full-length HSA prediction,
long-gap accuracy, alternative checkpoint comparisons and mixed-precision
execution have not been qualified.
