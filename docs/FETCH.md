# Retrieve a structure or sequence

```sh
cationize-protein fetch --pdb 1N5U --uniprot P02768 --output downloaded
```

Either identifier can be used alone. Existing directories are never
overwritten. The command downloads from the official RCSB PDB and UniProt
endpoints and retains the original response bytes.

| Output | Meaning |
|---|---|
| `1N5U.pdb` | Deposited coordinates, including all supplied chains, models and heterogens |
| `uniprot.json` | Complete UniProt response |
| `sequence.fasta` | Canonical sequence, which may include a signal peptide or propeptide |
| `mature.fasta` | Separately extracted chain only when one unique, exactly bounded Chain feature is annotated |
| `sequence_regions.csv` | Annotated chains, signal peptides and propeptides |
| `pdb_candidates.csv` | Linked experimental PDB entries, methods, resolutions and chain mappings |
| `FETCH.json`, `fetch.log` | Sources, SHA256 hashes and retrieval messages |

UniProt retrieval does not silently choose an experimental PDB or predicted
model. If there are multiple mature chains, no `mature.fasta` is chosen.
Check the desired construct, chain, sequence coverage, ligands and biological
assembly before preparing a native protein-only PDB.

HTTPS uses Python's trusted certificate configuration. If a Python.org macOS
installation reports `CERTIFICATE_VERIFY_FAILED`, configure its certificates
using the Python installation instructions or point `SSL_CERT_FILE` to a
trusted CA bundle. Certificate verification is never disabled. The local HSA
download check used the macOS bundle `/etc/ssl/cert.pem`.

Retrieval does not rebuild missing residues or normalize an experimental
structure. Use `prepare` for missing atoms, or the separately selected
[optional modeling workflow](MODELING.md) for a full structural candidate.

Data sources: [RCSB PDB](https://www.rcsb.org/),
[RCSB usage policy](https://www.rcsb.org/pages/usage-policy),
[UniProt REST API](https://www.uniprot.org/help/api_queries) and
[UniProt licence](https://www.uniprot.org/help/license).
