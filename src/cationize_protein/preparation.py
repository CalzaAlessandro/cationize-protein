"""Conservative PDBFixer preparation of a native, protein-only input.

Missing heavy atoms are inserted without moving supplied heavy coordinates.
PROPKA assigns initial histidine charge states; ring H and labels are updated.
Missing sequence, alternate conformers and nonstandard chemistry are not guessed.
"""
from __future__ import annotations

import importlib.metadata
import io
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

from .preflight import _NATIVE_PARENT, _SIDE, check_pdb


def prepare_pdb(source, output, *, pH=7.4, protonation='propka', uniprot_record=None):
    """Add missing atoms in existing native residues to a new output directory.

    PDBFixer builds candidate coordinates; only added heavy atoms are copied
    back into the source text. Initial histidine states and ring H are assigned
    by default; fixed mode preserves supplied states. The workflow subsequently
    minimizes and audits all supported native protonation changes with PROPKA.
    """
    source, output = Path(source), Path(output)
    if not math.isfinite(pH) or not 0<=pH<=14:
        raise ValueError('pH must be finite and between 0 and 14')
    if protonation not in {'propka','fixed'}:
        raise ValueError('Choose propka or fixed initial protonation')
    if output.exists() or output.is_symlink():
        raise ValueError('Preparation requires a new directory; existing output cannot be overwritten')
    raw = source.read_bytes()
    uniprot_raw = Path(uniprot_record).read_bytes() if uniprot_record is not None else None
    reference = json.loads(uniprot_raw) if uniprot_raw is not None else None
    if uniprot_raw is not None and not isinstance(reference, dict):
        raise ValueError('Provide a UniProt JSON object downloaded with fetch')
    before = check_pdb(source)
    for residue in before['residues']:
        if _NATIVE_PARENT.get(residue['residue_name'], residue['residue_name']) not in _SIDE:
            raise ValueError('PDBFixer preparation accepts native protein residues only; '
                             'cationized residues, ligands and other nonstandard groups are not replaced')
    blockers = [issue for issue in before['issues']
                if issue['severity'] == 'error' and issue['code'] != 'missing_heavy_atoms']
    if blockers:
        raise ValueError('Initial PDB cannot be repaired automatically: ' +
                         ' '.join(issue['message'] + ' ' + issue['action'] for issue in blockers))

    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    lines = raw.decode('utf-8').splitlines(keepends=True)
    # OpenMM uses canonical names internally. Original labels are never written
    # back from its topology: HID/HIE/HIP, ASH/GLH, CYX etc. stay as supplied.
    normalized = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            name = line[17:20].strip()
            line = line[:17] + f'{_NATIVE_PARENT.get(name, name):>3s}' + line[20:]
        normalized.append(line)
    fixer = PDBFixer(pdbfile=io.StringIO(''.join(normalized)))
    fixer.findMissingResidues()
    if fixer.missingResidues:
        missing = '; '.join(f'chain index {chain}, position {index}: {", ".join(names)}'
                            for (chain, index), names in sorted(fixer.missingResidues.items()))
        raise ValueError('Missing residues require a separately reviewed sequence/loop model: ' + missing)
    # No replacements, heterogen removal, H addition, solvation or CCD downloads.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    if fixer.missingAtoms or fixer.missingTerminals:
        fixer.addMissingAtoms(seed=0)
    residues = list(fixer.topology.residues())
    if len(residues) != len(before['residues']):
        raise ValueError('PDBFixer changed the residue count; preparation was not written')
    rendered = io.StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, rendered, keepIds=True)
    generated = [line for line in rendered.getvalue().splitlines()
                 if line.startswith(('ATOM  ', 'HETATM'))]
    atoms = list(fixer.topology.atoms())
    if len(generated) != len(atoms):
        raise ValueError('PDBFixer atom serialization is inconsistent; preparation was not written')
    serials = []
    for line in lines:
        if line[:6].strip() in {'ATOM', 'HETATM', 'TER'}:
            if line[:6].strip() == 'TER' and not line[6:11].strip():
                continue
            try:
                serials.append(int(line[6:11]))
            except ValueError:
                raise ValueError('Preparation requires valid numeric PDB atom serials') from None
    serial = max(serials, default=0)
    additions, added_atoms = {}, []
    present = [{a['name'] for a in residue['atoms']} for residue in before['residues']]
    res_coords = {(a.residue.index, a.name): (float(rec[30:38]), float(rec[38:46]), float(rec[46:54]))
                  for a, rec in zip(atoms, generated)}
    for atom, record in zip(atoms, generated):
        index = atom.residue.index
        if atom.name in present[index] or atom.element is None or atom.element.atomic_number == 1:
            continue
        if atom.name == 'OXT' and all((index, k) in res_coords for k in ('CA', 'C', 'O')):
            ca, c, o = (res_coords[index, k] for k in ('CA', 'C', 'O'))
            v_ca = tuple(ca[i] - c[i] for i in range(3))
            v_o = tuple(o[i] - c[i] for i in range(3))
            n_ca = math.sqrt(sum(x * x for x in v_ca))
            n_o = math.sqrt(sum(x * x for x in v_o))
            if n_ca > 1e-6 and n_o > 1e-6:
                bis = tuple(v_ca[i] / n_ca + v_o[i] / n_o for i in range(3))
                n_bis = math.sqrt(sum(x * x for x in bis))
                if n_bis > 1e-6:
                    oxt = tuple(c[i] - 1.25 * bis[i] / n_bis for i in range(3))
                    record = record[:30] + f'{oxt[0]:8.3f}{oxt[1]:8.3f}{oxt[2]:8.3f}' + record[54:]
        residue = before['residues'][index]
        original = lines[residue['atoms'][0]['line'] - 1]
        serial += 1
        if serial > 99999:
            raise ValueError('Added atoms exceed the PDB serial limit; use a separately prepared input')
        # Preserve the source residue identity and all supplied coordinate lines.
        record = (record[:6] + f'{serial:5d}' + record[11:17] + original[17:27] + record[27:] + '\n')
        last_line = residue['atoms'][-1]['line']
        additions.setdefault(last_line, []).append(record)
        added_atoms.append(dict(residue_ordinal=residue['ordinal'], chain_id=residue['chain_id'],
                                residue_id=residue['original_id'], residue_name=residue['residue_name'],
                                atom_name=atom.name))
    prepared = ''.join(line + ('' if line.endswith(('\n', '\r')) or i not in additions else '\n') +
                       ''.join(additions.get(i, [])) for i, line in enumerate(lines, 1))
    # Verify before creating persistent output; invalid inputs leave no directory.
    histidines=[];audit_files={};histidine_csv=None;disulfides=None
    with tempfile.TemporaryDirectory(prefix='cationize-prepare-') as tmp:
        candidate = Path(tmp)/'protein.pdb'
        candidate.write_bytes(prepared.encode('utf-8'))
        after = check_pdb(candidate)
        if after['ok'] and reference is not None:
            from .uniprot_disulfides import assign_disulfides
            prepared,disulfides=assign_disulfides(candidate,reference)
            candidate.write_text(prepared)
            after=check_pdb(candidate)
        if after['ok'] and protonation=='propka':
            from .initial_histidines import assign_histidines
            prepared,histidines=assign_histidines(candidate,Path(tmp)/'histidine_assignment',pH=pH)
            if histidines:
                original_heavy={(r['ordinal'],a['name']):a['coordinates_A'] for r in after['residues']
                                for a in r['atoms'] if not a['hydrogen']}
                candidate.write_text(prepared)
                after=check_pdb(candidate)
                actual_heavy={(r['ordinal'],a['name']):a['coordinates_A'] for r in after['residues']
                              for a in r['atoms'] if not a['hydrogen']}
                if original_heavy!=actual_heavy:
                    raise ValueError('Initial histidine assignment changed heavy coordinates or identities')
                audit=Path(tmp)/'histidine_assignment'
                histidine_csv=(audit/'histidines.csv').read_bytes()
                audit_files={path.name:path.read_bytes() for path in (audit/'audit').iterdir() if path.is_file()}
    if not after['ok']:
        raise ValueError('Repaired PDB did not pass validation: ' +
                         ' '.join(i['message'] for i in after['issues'] if i['severity'] == 'error'))
    if source.read_bytes() != raw:
        raise ValueError('The source PDB changed during preparation; retry from a stable input')
    log = ['PDB preparation: PASS', f'PDBFixer: {importlib.metadata.version("pdbfixer")}',
           f'Source SHA256: {before["source_sha256"]}',
           f'Residues: {before["residue_count"]}; added heavy atoms: {len(added_atoms)}',
           'Supplied heavy-atom coordinates and residue identifiers are preserved.',
           'Missing loops/residues, alternate conformers and nonstandard chemistry are not rebuilt.',
           ('Initial histidine charge states were assigned with PROPKA; neutral tautomers use supplied evidence or OpenMM H-bond geometry.'
            if histidines else 'No histidine H changes. The cationization workflow adds H and checks protonation with PROPKA.'),
           'Other native side-chain and terminal states remain supplied; the workflow checks them after minimization.',
           'Coordinate/SEQRES checks do not establish completeness against an external reference sequence.',
           'Added atoms require the subsequent force-field minimization before structural interpretation.']
    log += [f'ADDED: {a["chain_id"] or "(blank)"}:{a["residue_id"]} '
            f'{a["residue_name"]} {a["atom_name"]}' for a in added_atoms]
    if disulfides is not None:
        log += [f'UniProt record SHA256: {hashlib.sha256(uniprot_raw).hexdigest()}']
        log += disulfides['messages']
    log += [f'HISTIDINE: {r["chain_id"] or "(blank)"}:{r["residue_id"]} '
            f'pKa={r["pKa"]:.2f}, pH={pH:g}, {r["source_state"]} -> {r["state"]}; '
            f'{r["neutral_tautomer_method"]}'+(' (pKa within 0.5 of pH)' if r['borderline'] else '')
            for r in histidines]
    log += [f'{issue["severity"].upper()}: {issue["message"]} {issue["action"]}'
            for issue in after['issues']]
    output.mkdir(parents=True, exist_ok=False)
    (output/'original.pdb').write_bytes(raw)
    (output/'protein.pdb').write_bytes(prepared.encode('utf-8'))
    (output/'preparation.log').write_text('\n'.join(log) + '\n')
    if disulfides is not None:
        from .uniprot_disulfides import BOND_FIELDS, MAPPING_FIELDS
        (output/'uniprot.json').write_bytes(uniprot_raw)
        for filename,rows,fields in [('disulfides.csv',disulfides['rows'],BOND_FIELDS),
                                     ('uniprot_mapping.csv',disulfides['mapping'],MAPPING_FIELDS)]:
            with (output/filename).open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    if histidine_csv is not None:
        (output/'histidines.csv').write_bytes(histidine_csv)
        (output/'protonation').mkdir()
        for name,content in audit_files.items():(output/'protonation'/name).write_bytes(content)
    return dict(pdb=str(output/'protein.pdb'), added_atoms=added_atoms, messages=log,
                source_sha256=before['source_sha256'], histidines=histidines, disulfides=disulfides, ok=True)
