"""Bounded Slurm parameterization after a validated pilot; machine-readable status."""
import argparse
import csv
import fcntl
import io
import json
import os
import subprocess
from pathlib import Path
from .catalog import by_name
from .qm import digest
from .parameter_provenance import known_resp_source


def template_state(project, name):
    project = Path(project)
    library = project/'parameters'/(name+'.lib')
    fit_marker = project/'work/fits'/name/'FIT_COMPLETE.json'
    validation = fit_marker.with_name('VALIDATED_TEMPLATE.json')
    if not fit_marker.exists():
        return dict(name=name, status='incomplete_or_invalid' if library.exists() else 'awaiting_QM_RESP_validation',
                    reason='Fit marker absent')
    try:
        fit = json.loads(fit_marker.read_text())
        model = project/'work/models'/name
        if fit['library_sha256'] != digest(library):
            raise ValueError('Library hash differs')
        if fit['model_sha256'] != digest(model/'atom_mapping.json'):
            raise ValueError('Model hash differs')
        if not known_resp_source(project,fit['fit_code_sha256']):
            raise ValueError('Original fit source version is unavailable or changed')
        if abs(fit['charge']-by_name(name)['formal_charge']) > 3e-5:
            raise ValueError('Net charge differs from chemical state')
        for label, sha in fit['qm_input_sha256'].items():
            if digest(model/label/'opt_esp.gjf') != sha:
                raise ValueError('QM input changed')
        v = json.loads(validation.read_text())
        if v['fit_marker_sha256'] != digest(fit_marker) or v['library_sha256'] != digest(library):
            raise ValueError('Exported-template validation does not match this fit')
        if v['validation_code_sha256'] != digest(Path(__file__).with_name('parameter_validation.py')):
            raise ValueError('Export validation code changed')
        if v['backend_code_sha256'] != digest(Path(__file__).with_name('openmm_backend.py')):
            raise ValueError('OpenMM backend changed')
        return dict(name=name, status=v['status'], charge_e=fit['charge'],
                    library_sha256=digest(library), metrics=fit['metrics'])
    except (OSError, ValueError, KeyError) as exc:
        return dict(name=name, status='incomplete_or_invalid', reason=str(exc))


def collect(project):
    project = Path(project)
    target = project/'parameters/manifest.json'
    manifest = json.loads(target.read_text())
    for spec in manifest['templates']:
        if spec.get('reused_from'):
            try:
                if digest(project/'parameters'/(spec['name']+'.lib')) != spec['library_sha256']:
                    raise ValueError('Reused library changed')
                state = dict(status='reused_validated_HSA_template')
            except (OSError, ValueError, KeyError) as exc:
                state = dict(status='incomplete_or_invalid', reason=str(exc))
        else:
            state = template_state(project, spec['name'])
        spec.update(state)
        if 'reason' not in state:
            spec.pop('reason', None)
    temporary = target.with_suffix('.json.new')
    temporary.write_text(json.dumps(manifest, indent=2)+'\n')
    os.replace(temporary, target)
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=['name','parent','context','formal_charge','status','reason'])
    writer.writeheader()
    for spec in manifest['templates']:
        writer.writerow({k: spec.get(k, '') for k in writer.fieldnames})
    (project/'parameters/status.csv').write_text(output.getvalue())
    return manifest


def submit(args, project):
    result = subprocess.run(['sbatch', '--parsable', *map(str, args)], cwd=project,
                            text=True, capture_output=True, check=True)
    job = result.stdout.strip().split(';')[0]
    if not job.isdigit():
        raise ValueError('Unexpected sbatch job identifier: '+result.stdout)
    return job


def dispatch(project, pilot='GLX'):
    from .parameter_validation import validate
    project = Path(project).resolve()
    work = project/'work'
    work.mkdir(exist_ok=True)
    with (work/'parameter_dispatch.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        marker = work/'parameter_batch.json'
        if marker.exists():
            record = json.loads(marker.read_text())
            print('[QUEUE] Existing parameter batch: '+json.dumps(record), flush=True)
            return record
        validate(project, pilot)
        state = template_state(project, pilot)
        if state['status'] != 'RESP_capped_Amber_and_OpenMM_validated':
            raise RuntimeError('Pilot validation failed: '+str(state))
        manifest = collect(project)
        remaining = [s['name'] for s in manifest['templates'] if not s.get('reused_from')
                     and s['status'] != 'RESP_capped_Amber_and_OpenMM_validated']
        invalid = [s for s in manifest['templates'] if s['status'] == 'incomplete_or_invalid']
        if invalid:
            raise RuntimeError('Resolve incomplete or changed parameter exports before submitting: '+str(invalid))
        if not remaining:
            record = dict(status='complete', templates=[])
        else:
            array_file = work/'array_templates.json'
            if array_file.exists() and json.loads(array_file.read_text()) != remaining:
                raise RuntimeError('Existing array mapping differs; preserve it before a new batch')
            array_file.write_text(json.dumps(remaining, indent=2)+'\n')
            array = submit(['--array', f'0-{len(remaining)-1}%2', '--output', work/'parameter-%A_%a.log',
                            project/'config/parameterize.slurm', project], project)
            # Record the submitted array before requesting its dependent collector.
            record = dict(status='submitted', pilot=pilot, array_job_id=array,
                          templates=remaining, max_parallel=2)
            with marker.open('x') as stream:
                json.dump(record, stream, indent=2)
                stream.write('\n')
            collector = submit(['--dependency', 'afterany:'+array, '--output', work/'parameter_collect-%j.log',
                                project/'config/parameter_control.slurm', project, 'collect'], project)
            record['collector_job_id'] = collector
        marker.write_text(json.dumps(record, indent=2)+'\n')
        print('[QUEUE] '+json.dumps(record), flush=True)
        return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project', type=Path)
    parser.add_argument('action', choices=['pilot','collect'])
    args = parser.parse_args()
    result = dispatch(args.project) if args.action == 'pilot' else collect(args.project)
    print(json.dumps(result, indent=2))
