"""Conservative PDBFixer preparation of a native, protein-only input.

Only new heavy-atom records are inserted. Supplied atom records, coordinates,
residue states, identifiers and connectivity records are retained verbatim.
Missing sequence, alternate conformers and nonstandard chemistry are not guessed.
"""
from __future__ import annotations

import importlib.metadata
import io
from pathlib import Path
import tempfile

from .preflight import _NATIVE_PARENT, _SIDE, check_pdb


def prepare_pdb(source, output):
    """Add missing atoms in existing native residues to a new output directory.

    PDBFixer builds candidate coordinates; only added heavy atoms are copied
    back into the source text. No hydrogens, residues, solvent or ions are added.
    The workflow subsequently minimizes and audits protonation with PROPKA.
    """
    source, output = Path(source), Path(output)
    if output.exists() or output.is_symlink():
        raise ValueError('Preparation requires a new directory; existing output cannot be overwritten')
    raw = source.read_bytes()
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
    for atom, record in zip(atoms, generated):
        index = atom.residue.index
        if atom.name in present[index] or atom.element is None or atom.element.atomic_number == 1:
            continue
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
    with tempfile.TemporaryDirectory(prefix='cationize-prepare-') as tmp:
        candidate = Path(tmp)/'protein.pdb'
        candidate.write_bytes(prepared.encode('utf-8'))
        after = check_pdb(candidate)
    if not after['ok']:
        raise ValueError('Repaired PDB did not pass validation: ' +
                         ' '.join(i['message'] for i in after['issues'] if i['severity'] == 'error'))
    if source.read_bytes() != raw:
        raise ValueError('The source PDB changed during preparation; retry from a stable input')
    log = ['PDB preparation: PASS', f'PDBFixer: {importlib.metadata.version("pdbfixer")}',
           f'Source SHA256: {before["source_sha256"]}',
           f'Residues: {before["residue_count"]}; added heavy atoms: {len(added_atoms)}',
           'Supplied atom records, residue states, identifiers and connectivity records are preserved.',
           'Missing loops/residues, alternate conformers and nonstandard chemistry are not rebuilt.',
           'No hydrogens were added. The cationization workflow adds H and audits protonation with PROPKA.',
           'Histidines require explicit HID/HIE/HIP states or diagnostic ring H before run.',
           'Coordinate/SEQRES checks do not establish completeness against an external reference sequence.',
           'Added atoms require the subsequent force-field minimization before structural interpretation.']
    log += [f'ADDED: {a["chain_id"] or "(blank)"}:{a["residue_id"]} '
            f'{a["residue_name"]} {a["atom_name"]}' for a in added_atoms]
    log += [f'{issue["severity"].upper()}: {issue["message"]} {issue["action"]}'
            for issue in after['issues']]
    output.mkdir(parents=True, exist_ok=False)
    (output/'original.pdb').write_bytes(raw)
    (output/'protein.pdb').write_bytes(prepared.encode('utf-8'))
    (output/'preparation.log').write_text('\n'.join(log) + '\n')
    return dict(pdb=str(output/'protein.pdb'), added_atoms=added_atoms, messages=log,
                source_sha256=before['source_sha256'], ok=True)
