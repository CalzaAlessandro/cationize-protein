"""Map explicit UniProt intrachain disulfides to an unambiguous native PDB."""
from __future__ import annotations

import io
import json
import math

from .preflight import _NATIVE_PARENT, check_pdb

AA = dict(zip(('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL').split(),
              'ARNDCQEGHILKMFPSTWYV'))
ANNOTATION_FIELDS = ['accession', 'uniprot_start', 'uniprot_end', 'start_modifier',
                     'end_modifier', 'description', 'evidence']
BOND_FIELDS = ANNOTATION_FIELDS + ['pdb_chain', 'pdb_residue_1', 'pdb_residue_2',
                                  'distance_A', 'status']
MAPPING_FIELDS = ['accession', 'uniprot_position', 'pdb_chain', 'pdb_residue',
                  'residue_ordinal', 'residue_name']


def annotation_rows(record):
    """Retain database numbering and evidence, without a structural assignment."""
    rows = []
    features = record.get('features', [])
    if not isinstance(features, list):
        raise ValueError('UniProt annotation features must be a list')
    for feature in features:
        if not isinstance(feature, dict):
            raise ValueError('Each UniProt annotation must be a JSON object')
        if feature.get('type') != 'Disulfide bond':
            continue
        location = feature.get('location', {})
        if not isinstance(location, dict) or not isinstance(feature.get('description', ''), str):
            raise ValueError('UniProt disulfide annotation location/description has an invalid format')
        start, end = location.get('start', {}), location.get('end', {})
        if not isinstance(start, dict) or not isinstance(end, dict):
            raise ValueError('UniProt disulfide annotation positions must be JSON objects')
        rows.append(dict(accession=record.get('primaryAccession', ''),
            uniprot_start=start.get('value'), uniprot_end=end.get('value'),
            start_modifier=start.get('modifier', 'EXACT'),
            end_modifier=end.get('modifier', 'EXACT'),
            description=feature.get('description', ''),
            evidence=json.dumps(feature.get('evidences', []), separators=(',', ':'))))
    return rows


def _reference_pairs(record):
    try:
        sequence = record['sequence']['value'].upper()
        if not sequence or not sequence.isascii() or not sequence.isalpha():
            raise ValueError
    except (KeyError, TypeError, AttributeError, ValueError):
        raise ValueError('The UniProt JSON record must contain a valid canonical protein sequence') from None
    rows, partners = annotation_rows(record), {}
    for row in rows:
        a, b = row['uniprot_start'], row['uniprot_end']
        if (row['start_modifier'] != 'EXACT' or row['end_modifier'] != 'EXACT'
                or type(a) is not int or type(b) is not int):
            raise ValueError('UniProt disulfide positions must be exact; uncertain annotations require review')
        if 'interchain' in row['description'].lower() or a == b:
            raise ValueError('Interchain UniProt disulfides require an explicit assembly/partner assignment')
        if row['description'].strip().casefold() not in {'', 'intrachain'}:
            raise ValueError('Qualified UniProt disulfide annotations require review before assignment: '
                             + row['description'])
        if not 1 <= a < b <= len(sequence) or sequence[a-1] != 'C' or sequence[b-1] != 'C':
            raise ValueError(f'UniProt disulfide {a}–{b} does not identify two canonical cysteines')
        for source, target in ((a, b), (b, a)):
            if source in partners and partners[source] != target:
                raise ValueError(f'Conflicting UniProt disulfide partners for cysteine {source}')
            partners[source] = target
    return sequence, rows


def assign_disulfides(path, record):
    """Return PDB text and audit; retain heavy coordinates and require a real graph.

    Each continuous PDB chain must match a unique, exact contiguous region of
    the canonical sequence. No mutation, gap, isoform or assembly is inferred.
    Supplied SSBOND pairs must agree at any shared cysteine.
    """
    report = check_pdb(path)
    if not report['ok']:
        raise ValueError('UniProt disulfides require an otherwise valid, complete native PDB')
    original = path.read_text()
    lines = original.splitlines(keepends=True)
    serials = [int(line[6:11]) for line in lines if line.startswith(('ATOM  ', 'HETATM'))]
    if len(serials) != len(set(serials)):
        raise ValueError('UniProt disulfide assignment requires unique atom serials; duplicate atom serials make CONECT records ambiguous')
    sequence, annotations = _reference_pairs(record)
    groups = []
    for residue in report['residues']:
        key = (residue['model'], residue['segment'], residue['chain_id'])
        if not groups or key != groups[-1][0]:
            groups.append((key, []))
        groups[-1][1].append(residue)
    mapping, rows, pairs, messages = [], [], {}, []
    accession = record.get('primaryAccession', '(unspecified)')
    for (_, _, chain), residues in groups:
        observed = ''.join(AA.get(_NATIVE_PARENT.get(r['residue_name'], r['residue_name']), '?')
                           for r in residues)
        start = sequence.find(observed)
        if not observed or start < 0:
            raise ValueError(f'PDB chain {chain or "(blank)"} sequence does not match a contiguous region of UniProt {accession}; review construct, gaps or mutations')
        if sequence.find(observed, start+1) >= 0:
            raise ValueError(f'PDB chain {chain or "(blank)"} has an ambiguous UniProt sequence mapping')
        positions = {start+i+1: r for i, r in enumerate(residues)}
        mapping.extend(dict(accession=accession, uniprot_position=position, pdb_chain=chain,
            pdb_residue=r['original_id'], residue_ordinal=r['ordinal'], residue_name=r['residue_name'])
            for position, r in positions.items())
        messages.append(f'UniProt sequence: {accession} {start+1}–{start+len(residues)} -> '
                        f'chain {chain or "(blank)"}, {len(residues)} residues; exact unique match.')
        for annotation in annotations:
            a, b = annotation['uniprot_start'], annotation['uniprot_end']
            present = (a in positions, b in positions)
            if any(present) and not all(present):
                raise ValueError(f'UniProt disulfide {a}–{b}: one partner is missing from PDB chain {chain or "(blank)"}; review the partial construct')
            row = dict(annotation, pdb_chain=chain, pdb_residue_1='', pdb_residue_2='',
                       distance_A='', status='outside_construct')
            if all(present):
                first, second = positions[a], positions[b]
                for residue in (first, second):
                    if residue['residue_name'] not in {'CYS', 'CYX'}:
                        raise ValueError('A UniProt disulfide contradicts the supplied cysteine chemical state')
                sulfurs = [next(atom for atom in r['atoms'] if atom['name'] == 'SG')
                           for r in (first, second)]
                distance = math.dist(*(atom['coordinates_A'] for atom in sulfurs))
                if not 1.7 <= distance <= 2.5:
                    raise ValueError(f'UniProt disulfide {a}–{b} maps to {chain}:{first["original_id"]}–{second["original_id"]}, but SG–SG distance is {distance:.3f} Å; repair the geometry before preparation')
                pair = tuple(sorted((first['ordinal'], second['ordinal'])))
                pairs[pair] = (first, second)
                row.update(pdb_residue_1=first['original_id'], pdb_residue_2=second['original_id'],
                           distance_A=f'{distance:.3f}', status='verified')
                messages.append(f'UniProt disulfide: {accession} {a}–{b} -> '
                    f'{chain or "(blank)"}:{first["original_id"]}–{second["original_id"]}; SG–SG {distance:.3f} Å.')
            rows.append(row)
    expected = {}
    for pair in pairs:
        for a, b in (pair, pair[::-1]):
            expected[a] = b
    for declaration in report['ssbond_declarations']:
        a, b = declaration['residue_ordinals']
        if ((a in expected and expected[a] != b) or (b in expected and expected[b] != a)):
            raise ValueError('A supplied SSBOND conflicts with the mapped UniProt disulfide partner')

    bonded = {r['ordinal']: r for rs in pairs.values() for r in rs}
    atoms = {a['line']: r for r in bonded.values() for a in r['atoms']}
    removed = {int(lines[a['line']-1][6:11]) for r in bonded.values()
               for a in r['atoms'] if a['name'] == 'HG'}
    identities = {(r['chain_id'], r['original_id']) for r in bonded.values()}
    from .initial_histidines import _remove_connectivity_atoms
    updated = []
    for number, line in enumerate(lines, 1):
        if line.startswith(('ATOM  ', 'HETATM', 'ANISOU')):
            if int(line[6:11]) in removed:
                continue
            identity = (line[21].strip(), line[22:26].strip()+line[26].strip())
            if number in atoms or identity in identities:
                line = line[:17]+'CYX'+line[20:]
        updated.append(_remove_connectivity_atoms(line, removed))
    declarations = []
    existing = {tuple(sorted(d['residue_ordinals'])) for d in report['ssbond_declarations']}
    serial = max((int(line[7:10].strip() or 0) for line in lines if line.startswith('SSBOND')), default=0)
    connections = set()
    for line in updated:
        if line.startswith('CONECT'):
            numbers = [int(line[i:i+5]) for i in range(6, len(line.rstrip()), 5) if line[i:i+5].strip()]
            if numbers:connections.update((numbers[0], other) for other in numbers[1:])
    added_connections = []
    for pair, residues in pairs.items():
        if pair not in existing:
            serial += 1
            if serial > 999:raise ValueError('Disulfide declarations exceed the PDB SSBOND serial limit')
            line = list(' '*80);line[:6] = 'SSBOND';line[7:10] = f'{serial:3d}'
            for residue, begin in zip(residues, (11, 25)):
                atom_line = lines[residue['atoms'][0]['line']-1]
                line[begin:begin+3] = 'CYS';line[begin+4] = atom_line[21]
                line[begin+6:begin+11] = atom_line[22:27]
            line[59:65] = '  1555';line[66:72] = '  1555'
            line[73:78] = f'{math.dist(*(next(a["coordinates_A"] for a in r["atoms"] if a["name"]=="SG") for r in residues)):5.2f}'
            declarations.append(''.join(line)+'\n')
        sg = [int(lines[next(a['line'] for a in r['atoms'] if a['name']=='SG')-1][6:11]) for r in residues]
        for a, b in (sg, sg[::-1]):
            if (a, b) not in connections:
                added_connections.append(f'CONECT{a:5d}{b:5d}\n')
    insertion = next((i for i, line in enumerate(updated) if line[:6].strip()=='END'), len(updated))
    if insertion and not updated[insertion-1].endswith('\n'):updated[insertion-1] += '\n'
    updated[insertion:insertion] = added_connections
    prepared = ''.join(declarations+updated)
    from openmm.app import PDBFile
    pdb = PDBFile(io.StringIO(prepared))
    chains = {r.index+1: r.chain for r in pdb.topology.residues()}
    if any(chains[a] is not chains[b] for a, b in expected.items()):
        raise ValueError('A mapped intrachain UniProt disulfide crosses distinct topology chains')
    graph = {}
    for a, b in pdb.topology.bonds():
        if a.name == b.name == 'SG':
            graph.setdefault(a.residue.index+1, set()).add(b.residue.index+1)
            graph.setdefault(b.residue.index+1, set()).add(a.residue.index+1)
    if any(graph.get(a) != {b} for a, b in expected.items()):
        raise ValueError('The mapped UniProt disulfides conflict with the actual SG–SG bond graph')
    messages.append(f'UniProt disulfides: {len(pairs)} pairs verified in the bond graph; '
                    f'{len(removed)} thiol HG hydrogens removed; heavy coordinates retained.')
    outside = sum(row['status']=='outside_construct' for row in rows)
    if outside:messages.append(f'UniProt disulfides: {outside} annotations lie entirely outside the supplied construct.')
    return prepared, dict(rows=rows, mapping=mapping, messages=messages)
