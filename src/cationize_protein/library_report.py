"""Export numerical acceptance/charge tables for the local fitted library."""
import argparse
import csv
import json
import math
from pathlib import Path
from .off import read_units
from .qm import digest


def report(project):
    project = Path(project).resolve()
    manifest = json.loads((project/'parameters/manifest.json').read_text())
    rows, records = [], []
    accepted = {'reused_validated_HSA_template', 'RESP_capped_Amber_and_OpenMM_validated'}
    for spec in manifest['templates']:
        row = {k: spec[k] for k in ('name','parent','context','formal_charge','status')}
        row['modified_groups'] = '+'.join(spec['modified_groups'])
        if spec['status'] in accepted:
            lib = project/'parameters'/(spec['name']+'.lib')
            if digest(lib) != spec['library_sha256']:
                raise ValueError('Accepted library checksum differs: '+spec['name'])
            atoms = read_units(lib)[spec['name']]['atoms']
            charge = math.fsum(a['initial_template_charge_e'] for a in atoms)
            if abs(charge-spec['formal_charge']) > 3e-5:
                raise ValueError('Accepted template net charge differs: '+spec['name'])
            row.update(atoms=len(atoms), fitted_charge_e=charge,
                       charge_error_e=charge-spec['formal_charge'])
            record = dict(spec)
            fit = project/'work/fits'/spec['name']/'FIT_COMPLETE.json'
            if fit.is_file():
                record['fit'] = json.loads(fit.read_text())
                validation = fit.with_name('VALIDATED_TEMPLATE.json')
                record['export_validation'] = json.loads(validation.read_text())
                for metric in record['fit']['metrics']:
                    row['ESP_relative_RMS_'+metric['conformation']] = metric['relative_RMS']
            records.append(record)
        rows.append(row)
    columns = ['name','parent','context','modified_groups','formal_charge','fitted_charge_e',
               'charge_error_e','atoms','ESP_relative_RMS_alpha','ESP_relative_RMS_extended','status']
    with (project/'parameters/validation_summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    evidence = dict(catalog_templates=len(rows), accepted_templates=len(records), templates=records)
    (project/'parameters/validation_summary.json').write_text(json.dumps(evidence, indent=2)+'\n')
    if len(records) == len(rows):
        from .openmm_backend import template_xml
        (project/'parameters/cationized_residues.xml').write_text(template_xml(project/'parameters')+'\n')
    return dict(catalog_templates=len(rows), accepted_templates=len(records))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project', type=Path)
    args = parser.parse_args()
    print(json.dumps(report(args.project)))
