"""Explicit, optional full native models; conservative preparation stays separate."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import tempfile

from .optional_modelers import inspect_modelers, SOFTWARE
from .preflight import check_pdb, _NATIVE_PARENT

AA = dict(zip(('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL').split(),
              'ARNDCQEGHILKMFPSTWYV'))
CHECKPOINTS = {'fast': 'biohub/ESMFold2-Fast', 'full': 'biohub/ESMFold2'}


def read_fasta(path):
    rows = Path(path).read_text().splitlines()
    if sum(row.startswith('>') for row in rows) != 1 or not rows or not rows[0].startswith('>'):
        raise ValueError('Provide a single-chain FASTA with one header and one sequence')
    sequence = ''.join(row.strip() for row in rows[1:]).upper()
    if not sequence or any(aa not in AA.values() for aa in sequence):
        raise ValueError('The sequence must contain only the 20 canonical amino acids')
    if len(sequence) > 9999:
        raise ValueError('The sequence exceeds this PDB modeling interface residue limit')
    return sequence


def template_mapping(observed, target):
    """Map only missing residues, rejecting mutations and ambiguous placements."""
    early, cursor = [], 0
    for aa in observed:
        cursor = target.find(aa, cursor)
        if cursor < 0:
            raise ValueError('Template sequence is not a subsequence of the supplied complete sequence')
        early.append(cursor)
        cursor += 1
    late, cursor = [], len(target)
    for aa in reversed(observed):
        cursor = target.rfind(aa, 0, cursor)
        late.append(cursor)
    if early != list(reversed(late)):
        raise ValueError('Template-to-FASTA sequence mapping is ambiguous; supply a separately aligned model')
    return early


def _template(path, sequence):
    raw = Path(path).read_bytes()
    report = check_pdb(path)
    # A gap and missing atoms are precisely the reasons to request modeling.
    allowed = {'peptide_chain_gap', 'missing_heavy_atoms'}
    errors = [i['message'] for i in report['issues'] if i['severity']=='error' and i['code'] not in allowed]
    residues = report['residues']
    if errors or not residues or len({r['chain_id'] for r in residues}) != 1:
        raise ValueError('MODELLER template must be one unambiguous native chain: ' + ' '.join(errors))
    parents = [_NATIVE_PARENT.get(r['residue_name'], r['residue_name']) for r in residues]
    if any(name not in AA for name in parents):
        raise ValueError('MODELLER templates must contain native protein residues only')
    mapping = template_mapping(''.join(AA[name] for name in parents), sequence)
    records, rows = [], raw.decode().splitlines()
    for index, (residue, parent) in enumerate(zip(residues, parents), 1):
        for atom in residue['atoms']:
            if atom['hydrogen'] or atom['name'] == 'OXT':
                continue
            line = rows[atom['line']-1]
            records.append(f'ATOM  {len(records)+1:5d}' + line[11:17] + f'{parent:>3s} A{index:4d} ' + line[27:])
    identity = {(r['chain_id'], r['original_id']): i for i, r in enumerate(residues)}
    disulfides = []
    for declaration in report['ssbond_declarations']:
        partners = [identity[(p['chain_id'], p['original_id'])] for p in declaration['endpoints']]
        disulfides.append([mapping[i] for i in partners])
    # Explicit CYX pairs are retained; generic CYS proximity is not an oxidation assignment.
    for candidate in report['disulfide_candidates']:
        i, j = [n-1 for n in candidate['residue_ordinals']]
        if residues[i]['residue_name'] == residues[j]['residue_name'] == 'CYX':
            pair = [mapping[i], mapping[j]]
            if pair not in disulfides:
                disulfides.append(pair)
    aligned = ['-'] * len(sequence)
    for parent, position in zip(parents, mapping):
        aligned[position] = AA[parent]
    stream = io.StringIO(newline='')
    table = csv.writer(stream)
    table.writerow(['source_chain', 'source_residue', 'source_name', 'target_position'])
    for residue, position in zip(residues, mapping):
        table.writerow([residue['chain_id'], residue['original_id'], residue['residue_name'], position+1])
    pir = (f'>P1;template\nstructureX:template:1:A:{len(residues)}:A::::\n'
           f'{"".join(aligned)}*\n>P1;target\nsequence:target::::::::\n{sequence}*\n')
    return {'template.pdb': ('\n'.join(records)+'\nTER\nEND\n').encode(),
            'original.pdb': raw, 'alignment.ali': pir.encode(), 'residue_mapping.csv': stream.getvalue().encode()}, disulfides


def _run_worker(python, request, directory):
    # This intentionally runs in the optional software's own Python environment.
    with (directory/'modeling.log').open('w') as log:
        return subprocess.run([str(python), '-I', '-B', str(Path(__file__).with_name('modeling_worker.py')),
            str(request)], cwd=directory, stdout=log, stderr=subprocess.STDOUT, text=True)


def _validate_candidate(path, sequence, disulfides):
    report = check_pdb(path)
    errors = [i['message'] for i in report['issues'] if i['severity']=='error']
    residues = report['residues']
    observed = ''.join(AA.get(_NATIVE_PARENT.get(r['residue_name'], r['residue_name']), '?') for r in residues)
    if observed != sequence:
        errors.append('Modeled sequence does not exactly match the supplied complete sequence')
    if len({(r['model'], r['segment'], r['chain_id']) for r in residues}) != 1:
        errors.append('The modeled protein must be one continuous chain')
    for i, j in disulfides:
        sulfurs = [[a['coordinates_A'] for a in residues[n]['atoms'] if a['name']=='SG']
                   if n < len(residues) else [] for n in (i,j)]
        if any(len(s)!=1 or s[0] is None for s in sulfurs) or not 1.7 <= math.dist(sulfurs[0][0], sulfurs[1][0]) <= 2.5:
            errors.append(f'Template disulfide at target positions {i+1}, {j+1} was not retained')
    if errors:
        raise ValueError('Modeled structure failed sequence/geometry validation: ' + ' '.join(errors))
    raw = path.read_text()
    if disulfides:
        declarations = []
        for number, pair in enumerate(disulfides, 1):
            line = list(' '*80); line[:6] = 'SSBOND'; line[7:10] = f'{number:3d}'
            for index, start in zip(pair, (11,25)):
                residue = residues[index]
                line[start:start+3] = 'CYS'; line[start+4] = residue['chain_id'] or ' '
                line[start+6:start+10] = f'{int(residue["original_id"]):4d}'
            declarations.append(''.join(line).rstrip())
        raw = '\n'.join(declarations) + '\n' + raw
    return raw, report


def build_model(fasta, output, *, method, template=None, esmfold_model='fast', device='auto',
                python=None, revision='main'):
    """Create a full structural candidate, never silently replace an input PDB."""
    output, fasta = Path(output).absolute(), Path(fasta)
    if output.exists() or output.is_symlink():
        raise ValueError('Modeling requires a new directory; existing output is preserved')
    if method not in SOFTWARE:
        raise ValueError('Choose modeling method modeller or esmfold')
    if esmfold_model not in CHECKPOINTS or device not in {'auto','cpu','cuda'}:
        raise ValueError('Unknown ESMFold2 model or compute device')
    sequence = read_fasta(fasta)
    inputs = {'sequence.fasta': fasta.read_bytes()}
    disulfides = []
    if method == 'modeller':
        if template is None:
            raise ValueError('MODELLER requires a native template PDB in addition to the FASTA')
        generated, disulfides = _template(template, sequence)
        inputs.update(generated)
    elif template is not None:
        raise ValueError('ESMFold2 predicts from sequence; a template is only used by MODELLER')
    if python:
        import os
        environment = dict(os.environ)
        environment[SOFTWARE[method]['environment']] = str(python)
        capability = inspect_modelers(environ=environment)[method]
    else:
        capability = inspect_modelers()[method]
    if not capability['available']:
        raise ValueError(f'{SOFTWARE[method]["name"]}: {capability["reason"]}. '
                         f'Installation: {SOFTWARE[method]["url"]}; PDBFixer preparation remains available.')
    output.parent.mkdir(parents=True, exist_ok=True)
    messages = [f'Optional modeling: {SOFTWARE[method]["name"]}; {len(sequence)} residues.',
                'A full native structural candidate is generated; supplied coordinates are not frozen.',
                'Next: cationize-protein prepare MODEL_DIRECTORY/protein.pdb --output prepared']
    with tempfile.TemporaryDirectory(prefix='.cationize-model-', dir=output.parent) as tmp:
        staging = Path(tmp)/'model'; staging.mkdir()
        for name, raw in inputs.items():
            (staging/name).write_bytes(raw)
        request = dict(method=method, sequence=sequence, checkpoint=CHECKPOINTS[esmfold_model] if method=='esmfold' else None,
                       device=device, seed=-8123 if method=='modeller' else 0, revision=revision, disulfides=disulfides)
        (staging/'request.json').write_text(json.dumps(request, indent=2)+'\n')
        try:
            result = _run_worker(capability['python'], staging/'request.json', staging)
            # In production stdout goes straight to the log; fixtures may return it.
            if result.stdout:
                (staging/'modeling.log').write_text(result.stdout)
            if result.returncode:
                raise ValueError('Optional modeling failed; see modeling.log for licence, weights, API or hardware errors')
            raw, checked = _validate_candidate(staging/'candidate.pdb', sequence, disulfides)
            (staging/'protein.pdb').write_text(raw)
            rows = [line for line in raw.splitlines() if line.startswith(('ATOM  ','HETATM')) and line[12:16].strip()=='CA']
            if method == 'esmfold':
                confidence = io.StringIO(newline=''); table = csv.writer(confidence)
                table.writerow(['residue_position','residue_name','plddt_0_to_1'])
                for index, line in enumerate(rows, 1):
                    value = float(line[60:66])
                    if not math.isfinite(value) or not 0 <= value <= 1:
                        raise ValueError('ESMFold2 pLDDT is outside the documented 0-1 scale')
                    table.writerow([index,line[17:20].strip(),value])
                (staging/'confidence.csv').write_text(confidence.getvalue())
            if fasta.read_bytes() != inputs['sequence.fasta'] or (template and Path(template).read_bytes() != inputs['original.pdb']):
                raise ValueError('An input changed during modeling; retry from stable inputs')
            metadata = json.loads((staging/'worker_metadata.json').read_text()) if (staging/'worker_metadata.json').is_file() else {}
            record = dict(schema_version=1, status='structural_candidate', method=method,
                sequence=sequence, seed=request['seed'], checkpoint=request['checkpoint'], revision=revision,
                runtime=metadata, retained_template_disulfides=disulfides,
                residue_count=checked['residue_count'], messages=messages,
                files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in staging.iterdir() if p.is_file()})
            (staging/'MODEL.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
        except (ValueError, OSError, RuntimeError) as exc:
            (staging/'protein.pdb').unlink(missing_ok=True)
            (staging/'MODELING_FAILED.json').write_text(json.dumps(dict(status='failed',method=method,error=str(exc)))+'\n')
            if output.exists() or output.is_symlink():
                raise ValueError('Model output appeared during execution; existing files were preserved') from exc
            staging.rename(output)
            raise ValueError(f'{exc}; diagnostics: {output}') from exc
        if output.exists() or output.is_symlink():
            raise ValueError('Model output appeared during execution; existing files were preserved')
        staging.rename(output)
    return record
