# Retrieve a structure or sequence

After activating the environment, use `cationize-protein fetch` with a PDB ID,
a UniProt accession or entry name, or both. At least one identifier is required.
Each example uses a new output directory.

## PDB only

```sh
cationize-protein fetch --pdb 1N5U --output downloaded_pdb
```

This saves the deposited coordinates as `1N5U.pdb`, including all supplied
chains, models and ligands. It also writes `FETCH.json` and `fetch.log`.
It does not retrieve a UniProt sequence.

## UniProt only

```sh
cationize-protein fetch --uniprot P02768 --output downloaded_sequence
```

This saves the UniProt record, its canonical sequence, sequence annotations
and a table of linked PDB entries. It saves `mature.fasta` only when the record
annotates one mature chain with exact boundaries. Linked PDB entries are
listed for review; their coordinates are not downloaded automatically.

## PDB and UniProt together

```sh
cationize-protein fetch --pdb 1N5U --uniprot P02768 --output downloaded
```

This combines the PDB and UniProt outputs. The identifiers are independent:
the command retrieves each requested record but does not verify that they
describe the same protein, chain or construct. Check that correspondence
before using the sequence with the structure.

## Output files

The command downloads from the official RCSB PDB and UniProt endpoints and
retains the original coordinate and JSON response bytes. Existing output
directories are never overwritten.

| Output | Created for | Meaning |
|---|---|---|
| `1N5U.pdb` | PDB | Deposited coordinates; filename follows the requested PDB ID |
| `uniprot.json` | UniProt | Complete UniProt response |
| `sequence.fasta` | UniProt | Canonical sequence, which may include a signal peptide or propeptide |
| `mature.fasta` | UniProt, when uniquely annotated | Separately extracted chain with exact annotated boundaries |
| `sequence_regions.csv` | UniProt | Annotated chains, signal peptides and propeptides |
| `disulfides.csv` | UniProt | Annotated disulfide pairs in canonical numbering, with location qualifiers and source evidence |
| `pdb_candidates.csv` | UniProt | Linked experimental PDB entries, methods, resolutions and chain mappings |
| `FETCH.json` | Every fetch | Requested identifiers, source URLs, SHA256 hashes and retrieval messages |
| `fetch.log` | Every fetch | Human-readable retrieval messages |

UniProt retrieval does not silently choose an experimental PDB or predicted
model. If there are multiple mature chains, no `mature.fasta` is chosen.
Check the desired construct, chain, sequence coverage, ligands and biological
assembly before preparing a native protein-only PDB.

Retrieval does not rebuild missing residues or normalize an experimental
structure. Use `prepare` for missing atoms, or the separately selected
[optional modeling workflow](MODELING.md) for a full structural candidate.

To apply known intrachain disulfides, pass `uniprot.json` to
`prepare --uniprot-record downloaded/uniprot.json`. Preparation maps the
sequence to the PDB and checks sulfur geometry; downloading both records
alone does not establish that correspondence. See
[disulfide preparation](PDB_PREPARATION.md#disulfides-from-uniprot).

## HTTPS certificates

HTTPS uses Python's trusted certificate configuration. If a Python.org macOS
installation reports `CERTIFICATE_VERIFY_FAILED`, configure its certificates
using the Python installation instructions or point `SSL_CERT_FILE` to a
trusted CA bundle. Certificate verification is never disabled. The local HSA
download check used the macOS bundle `/etc/ssl/cert.pem`.

Data sources: [RCSB PDB](https://www.rcsb.org/),
[RCSB usage policy](https://www.rcsb.org/pages/usage-policy),
[UniProt REST API](https://www.uniprot.org/help/api_queries) and
[UniProt licence](https://www.uniprot.org/help/license).
