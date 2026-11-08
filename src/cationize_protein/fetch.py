"""Retrieve deposited structures and sequence records from public databases."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tempfile
from urllib.request import Request, urlopen


def download(url):
    request = Request(url, headers={'User-Agent': 'cationize-protein/0.1'})
    with urlopen(request, timeout=30) as response:
        content = response.read(25 * 1024 * 1024 + 1)
    if len(content) > 25 * 1024 * 1024:
        raise ValueError('Database response exceeds the 25 MiB download limit')
    return content


def _csv(rows, fields):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def _pdb_uniprot_candidates(pdb_text):
    """Extract UniProt cross-references from PDB DBREF / DBREF1 / DBREF2 records."""
    rows = []
    lines = pdb_text.splitlines()
    for index, line in enumerate(lines):
        if line.startswith('DBREF ') and len(line) >= 42 and line[26:32].strip() == 'UNP':
            rows.append(dict(
                chain=line[12:13].strip(),
                pdb_start=line[14:18].strip(),
                pdb_end=line[20:24].strip(),
                uniprot_accession=line[33:41].strip(),
                uniprot_id=line[42:54].strip(),
                uniprot_start=line[55:60].strip(),
                uniprot_end=line[62:67].strip(),
            ))
        elif line.startswith('DBREF1') and len(line) >= 47 and line[26:32].strip() == 'UNP':
            acc = line[47:67].strip()
            unp_start, unp_end = '', ''
            if index + 1 < len(lines) and lines[index + 1].startswith('DBREF2'):
                nxt = lines[index + 1]
                acc = nxt[18:40].strip() or acc
                unp_start = nxt[45:55].strip()
                unp_end = nxt[57:67].strip()
            rows.append(dict(
                chain=line[12:13].strip(),
                pdb_start=line[14:18].strip(),
                pdb_end=line[20:24].strip(),
                uniprot_accession=acc,
                uniprot_id=line[33:47].strip(),
                uniprot_start=unp_start,
                uniprot_end=unp_end,
            ))
    return rows


def fetch(output, *, pdb_id=None, uniprot_id=None):
    """Save downloaded bytes without choosing a model or altering a sequence."""
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise ValueError('Fetch requires a new directory; existing output is preserved')
    if not pdb_id and not uniprot_id:
        raise ValueError('Provide a PDB ID, a UniProt accession, or both')
    if pdb_id and not re.fullmatch(r'[0-9][A-Za-z0-9]{3}', pdb_id):
        raise ValueError('PDB ID must contain four characters, starting with a digit')
    if uniprot_id and not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{1,39}', uniprot_id):
        raise ValueError('Invalid UniProt accession or entry name')
    pdb_id = pdb_id.upper() if pdb_id else None
    uniprot_id = uniprot_id.upper() if uniprot_id else None
    files, sources = {}, {}
    messages = ['Database retrieval; original records are preserved.']
    if pdb_id:
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        raw = download(url)
        text = raw.decode('utf-8')
        atoms = [line for line in text.splitlines()
                 if line.startswith(('ATOM  ', 'HETATM'))]
        try:
            if not atoms or any(len(line) < 54 or not all(math.isfinite(float(line[a:b]))
                    for a, b in ((30, 38), (38, 46), (46, 54))) for line in atoms):
                raise ValueError
        except ValueError:
            raise ValueError('The downloaded record is not a valid PDB coordinate file') from None
        name = f'{pdb_id}.pdb'
        files[name], sources[name] = raw, url
        messages.append(f'PDB: {pdb_id}; deposited chains, models and heterogens retained.')
        unp_rows = _pdb_uniprot_candidates(text)
        if unp_rows:
            files['uniprot_candidates.csv'] = _csv(
                unp_rows,
                ['chain', 'pdb_start', 'pdb_end', 'uniprot_accession', 'uniprot_id', 'uniprot_start', 'uniprot_end'],
            )
            accessions = sorted({r['uniprot_accession'] for r in unp_rows if r['uniprot_accession']})
            messages.append(f'PDB DBREF UniProt mappings: {", ".join(accessions)} (saved in uniprot_candidates.csv).')
            if not uniprot_id and len(accessions) == 1:
                uniprot_id = accessions[0].upper()
                messages.append(f'Automatically retrieving unique mapped UniProt record: {uniprot_id}.')
    if uniprot_id:
        url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.json'
        raw = download(url)
        try:
            record = json.loads(raw)
            sequence = record['sequence']['value'].upper()
            if not re.fullmatch('[A-Z]+', sequence):
                raise ValueError
        except (ValueError, KeyError, TypeError, AttributeError):
            raise ValueError('The UniProt response does not contain a valid protein sequence') from None
        accession = record.get('primaryAccession', uniprot_id)
        files['uniprot.json'], sources['uniprot.json'] = raw, url
        from .uniprot_disulfides import annotation_rows, ANNOTATION_FIELDS
        disulfides = annotation_rows(record)
        files['disulfides.csv'] = _csv(disulfides, ANNOTATION_FIELDS)
        messages.append(f'UniProt disulfides: {len(disulfides)} annotations saved with canonical numbering and evidence.')
        files['sequence.fasta'] = f'>{accession} canonical sequence\n{sequence}\n'.encode()
        regions = []
        for feature in record.get('features', []):
            kind = feature.get('type')
            if kind not in {'Chain', 'Signal', 'Propeptide'}:
                continue
            location = feature.get('location', {})
            start, end = location.get('start', {}), location.get('end', {})
            regions.append(dict(type=kind, start=start.get('value'), end=end.get('value'),
                start_modifier=start.get('modifier', 'EXACT'), end_modifier=end.get('modifier', 'EXACT'),
                description=feature.get('description', '')))
        files['sequence_regions.csv'] = _csv(regions,
            ['type', 'start', 'end', 'start_modifier', 'end_modifier', 'description'])
        chains = [r for r in regions if r['type'] == 'Chain']
        if len(chains) == 1:
            chain = chains[0]
            start, end = chain['start'], chain['end']
            if (isinstance(start, int) and isinstance(end, int) and 1 <= start <= end <= len(sequence)
                    and chain['start_modifier'] == chain['end_modifier'] == 'EXACT'):
                files['mature.fasta'] = (f'>{accession} annotated chain {start}-{end}\n'
                                        f'{sequence[start-1:end]}\n').encode()
                messages.append(f'Unique annotated mature chain: {start}-{end}; saved separately.')
        candidates = []
        for entry in record.get('uniProtKBCrossReferences', []):
            if entry.get('database') != 'PDB':
                continue
            properties = {p.get('key'): p.get('value') for p in entry.get('properties', [])}
            candidates.append(dict(pdb_id=entry.get('id', ''), method=properties.get('Method', ''),
                resolution=properties.get('Resolution', ''), chains=properties.get('Chains', '')))
        files['pdb_candidates.csv'] = _csv(candidates, ['pdb_id', 'method', 'resolution', 'chains'])
        messages.append(f'UniProt: {accession}; canonical sequence {len(sequence)} residues; '
                        f'{len(candidates)} PDB candidates. No structure was selected automatically.')
    record = dict(schema_version=1, pdb_id=pdb_id, uniprot_id=uniprot_id,
                  selected_structure=f'{pdb_id}.pdb' if pdb_id else None,
                  sources=sources, files_sha256={name: hashlib.sha256(raw).hexdigest()
                                                for name, raw in files.items()}, messages=messages)
    files['FETCH.json'] = (json.dumps(record, indent=2) + '\n').encode()
    files['fetch.log'] = ('\n'.join(messages) + '\n').encode()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.cationize-fetch-', dir=output.parent) as tmp:
        staging = Path(tmp)/'record'
        staging.mkdir()
        for name, content in files.items():
            (staging/name).write_bytes(content)
        if output.exists() or output.is_symlink():
            raise ValueError('Fetch output was created by another process; existing files are preserved')
        staging.rename(output)
    return record
