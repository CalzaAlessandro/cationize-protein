"""Command line entry point for portable progressive protein cationization."""
from __future__ import annotations

import argparse
import importlib.util
import json
import logging
from pathlib import Path
import shutil
import sys

from .electrostatics import APBSConfig
from .workflow import WorkflowConfig, run
from .resources import default_parameters

LOG=logging.getLogger(__name__)
DEFAULT_PROBE_RADIUS_A = 1.369  # Na+ LJ half-Rmin used in the HSA reference protocol.
_CHOICES = {'engine': ('openmm', 'amber', 'sander', 'pmemd'),
            'solver': ('apbs', 'pbsa', 'coulomb'),
            'force_field': ('ff14SB', 'ff19SB', 'ff99SBildn')}


def parser():
    command = argparse.ArgumentParser(
        prog='cationize-protein',
        description='Serial progressive carboxyl amidation, steric rotation, restrained minimization and PROPKA auditing. '
                    'Physical models and export formats are explicit; Slurm is optional and external.')
    subcommands = command.add_subparsers(dest='command', required=True)
    retrieve = subcommands.add_parser('fetch', help='Download a deposited PDB and/or UniProt sequence record')
    retrieve.add_argument('--pdb', help='Four-character PDB ID; deposited file is kept unchanged')
    retrieve.add_argument('--uniprot', help='UniProt accession or entry name; no structure is chosen automatically')
    retrieve.add_argument('--output', type=Path, required=True, help='New directory for downloaded records')
    modeling = subcommands.add_parser('model', help='Optional full native model with installed MODELLER or ESMFold2; then use prepare')
    modeling.add_argument('--method', choices=['modeller','esmfold'], required=True, help='Optional modeler; esmfold uses ESMFold2')
    modeling.add_argument('--fasta', type=Path, required=True, help='One complete native chain, 20 canonical amino acids')
    modeling.add_argument('--pdb', type=Path, help='Native incomplete template PDB; required for MODELLER only')
    modeling.add_argument('--output', type=Path, required=True, help='New directory for model, provenance and logs')
    modeling.add_argument('--python', help='Python of the optional modeler; otherwise its environment variable or current Python')
    modeling.add_argument('--esmfold-model', choices=['fast','full'], default='fast', help='ESMFold2 checkpoint; default: fast, single sequence')
    modeling.add_argument('--model-revision', default='main', help='Hugging Face checkpoint revision; resolved revision is recorded when available')
    modeling.add_argument('--device', choices=['auto','cpu','cuda'], default='auto', help='ESMFold2 device; auto selects CUDA when available, otherwise CPU')
    fitting = subcommands.add_parser('parameterize', help='Develop new templates sequentially without Slurm')
    fitting.add_argument('--project', type=Path, required=True, help='Parameter-development project containing work/models')
    fitting.add_argument('--names', nargs='+', help='Template names in execution order; default: all 28 new-fit states')
    check=subcommands.add_parser('check',help='Read-only protein PDB adequacy check')
    check.add_argument('pdb',type=Path)
    check.add_argument('--json',type=Path,help='Write a new machine-readable report')
    prepare=subcommands.add_parser('prepare', help='Repair initial heavy atoms with PDBFixer and assign histidines with PROPKA')
    prepare.add_argument('pdb',type=Path,help='Native protein-only PDB; supplied coordinates and residue states are preserved')
    prepare.add_argument('--output',type=Path,required=True,help='New directory for original.pdb, protein.pdb and preparation.log')
    prepare.add_argument('--ph',type=float,default=7.4,help='Initial histidine PROPKA pH (default: 7.4)')
    prepare.add_argument('--protonation',choices=['propka','fixed'],default='propka',help='Default: PROPKA assigns histidine charge; fixed retains supplied states')
    prepare.add_argument('--uniprot-record',type=Path,help='Optional local uniprot.json from fetch: map exact intrachain disulfides by sequence, verify SG geometry and write CYX/SSBOND/CONECT')
    doctor=subcommands.add_parser('doctor',help='Report Python dependencies and optional external programs')
    doctor.add_argument('--json',action='store_true',help='Print the dependency report as JSON; write no files')
    execute = subcommands.add_parser('run', help='Create or automatically resume a verified progressive run')
    source = execute.add_mutually_exclusive_group(required=True)
    source.add_argument('--topology', type=Path, help='Prepared protein-only Amber parm7/prmtop')
    source.add_argument('--pdb', type=Path, help='Native PDB with complete heavy atoms and explicit histidine states')
    execute.add_argument('--coordinates', type=Path, help='Amber rst7/inpcrd, required with --topology')
    execute.add_argument('--parameters', type=Path, default=default_parameters(), help='Residue libraries; default: bundled validated library')
    execute.add_argument('--output', type=Path, required=True, help='New output directory or the same verified run to resume')
    execute.add_argument('--probe-radius', type=float, default=DEFAULT_PROBE_RADIUS_A,
                         help='Na probe Lennard-Jones half-Rmin in Å (default: 1.369, HSA reference model); distinct from PB ion radius')
    execute.add_argument('--percentages', nargs='+', default=['25', '50', '75', '100'],
                         help='Percentages of all carboxyl groups, rounded half up; baseline 0 is always exported')
    execute.add_argument('--force-field',choices=_CHOICES['force_field'],help='Optional: choose interactively for a new terminal run; noninteractive default ff14SB or saved preference')
    execute.add_argument('--engine',choices=_CHOICES['engine'],help='Optional minimizer: choose at runtime; noninteractive default openmm or saved preference')
    execute.add_argument('--engine-executable',help='Installed serial sander/pmemd executable; never downloaded')
    execute.add_argument('--solver',choices=_CHOICES['solver'],help='Optional potential model: choose at runtime; noninteractive default apbs or saved preference. Coulomb is an approximation')
    execute.add_argument('--apbs',help='APBS executable; default: private installation or PATH')
    execute.add_argument('--pbsa',help='Installed Amber PBSA executable when --solver pbsa is selected')
    execute.add_argument('--pbsa-nonlinear',action='store_true',help='PBSA only: full nonlinear PB with damped SOR; default PBSA is linear ICCG. APBS remains nonlinear')
    execute.add_argument('--protonation',choices=['propka','fixed'],default='propka',help='Default: PROPKA at baseline and every chemical step; fixed reproduces supplied native states')
    execute.add_argument('--ph',type=float,default=7.4,help='PROPKA pH (default: 7.4); introduced ammonium groups remain +1')
    execute.add_argument('--formats',nargs='+',choices=['amber','openmm','gromacs'],default=['amber'],help='Physical parameter exports; default: amber. PDB and OpenMM XML/JSON are always retained')
    execute.add_argument('--temperature', type=float, default=298., help='PB temperature in K (default: 298)')
    execute.add_argument('--solute-dielectric', type=float, default=2.)
    execute.add_argument('--solvent-dielectric', type=float, default=78.5)
    execute.add_argument('--salt', type=float, default=.15, help='PB ranking: symmetric monovalent salt concentration in M; minimization uses unscreened GBn2')
    execute.add_argument('--pb-ion-radius', type=float, default=2., help='PB electrolyte ion radius in Å')
    execute.add_argument('--pb-spacing', type=float, default=.75, help='Maximum fine-grid spacing in Å')
    execute.add_argument('--sampling-spacing', type=float, default=1., help='Local 5 Å sphere lattice spacing in Å')
    execute.add_argument('--keep-maps', action='store_true', help='Retain every APBS potential DX map')
    execute.add_argument('--platform', choices=['CPU', 'Reference', 'CUDA', 'OpenCL'], default='CPU')
    execute.add_argument('--beam-width', type=int, default=8)
    execute.add_argument('--max-iterations', type=int, default=20000)
    execute.add_argument('--restraint', type=float, default=10., help='Native backbone restraint in kcal/(mol Å²)')
    execute.add_argument('--tolerance', type=float, default=.1, help='Minimization RMS force tolerance in kcal/(mol Å)')
    execute.add_argument('--export-amber', action='store_true',help=argparse.SUPPRESS)
    execute.add_argument('--verbose', action='store_true')
    return command


def installed_settings():
    path=Path(sys.prefix)/'cationize-protein-settings.json'
    if not path.exists():return {}
    if path.is_symlink():raise ValueError('Installation settings must not be symlinks')
    value=json.loads(path.read_text())
    if not isinstance(value,dict) or value.get('schema_version')!=1:
        raise ValueError('Unsupported installation settings format')
    return value


def _available_choices(args, project):
    """Inspect current capabilities again, including software installed later."""
    from .dependencies import find_executable
    from .forcefields import force_field_status
    openmm = importlib.util.find_spec('openmm') is not None
    serial = {name: bool(find_executable(name, project)) for name in ('sander', 'pmemd')}
    if args.engine_executable and shutil.which(args.engine_executable):
        # The actual backend validates this explicit serial executable later.
        for name in serial:
            serial[name] = True
    solvers = {name: bool(shutil.which(getattr(args, name))) if getattr(args, name) else
               bool(find_executable(name, project)) for name in ('apbs', 'pbsa')}
    return {'engine': {'openmm':openmm, 'amber':any(serial.values()), **serial},
            'solver': {**solvers, 'coulomb':True},
            'force_field': {name:bool(force_field_status(name).get('available')) if openmm else False
                            for name in _CHOICES['force_field']}}


def run_choices(args, settings, project):
    """Select only new-run choices; committed run settings take precedence on resume.

    This reads configuration only. Full input and checkpoint validation remains
    authoritative in workflow.validate_inputs before scientific execution.
    """
    defaults = {'engine':'openmm', 'solver':'apbs', 'force_field':'ff14SB'}
    output = Path(args.output)
    previous = {}
    if output.is_symlink():
        raise ValueError('Output directories must not be symlinks')
    if output.exists():
        marker = output/'RUN.json'
        if not output.is_dir() or marker.is_symlink() or not marker.is_file():
            raise ValueError('Existing output is not a recognized committed run; use a new directory')
        previous = json.loads(marker.read_text())
        if not isinstance(previous, dict) or previous.get('schema_version') != 2:
            raise ValueError('Existing run has an unsupported configuration format')
        if any(key not in previous for key in defaults):
            raise ValueError('Existing run is missing its scientific choices')
    selected = {key:getattr(args,key) or previous.get(key) or settings.get(key,default)
                for key,default in defaults.items()}
    for key,value in selected.items():
        if value not in _CHOICES[key]:
            raise ValueError(f'Unknown saved {key}: {value!r}')
    pending = [key for key in defaults if getattr(args,key) is None]
    if not previous and pending and sys.stdin.isatty():
        available = _available_choices(args, project)
        labels = {'engine':'Minimization engine', 'solver':'Electrostatic model', 'force_field':'Protein force field'}
        print('\nChoose the scientific settings for this run. Enter accepts the displayed default.')
        for key in pending:
            print(f'\n{labels[key]}:')
            options = list(_CHOICES[key])
            for number,value in enumerate(options,1):
                status = 'available' if available[key][value] else 'unavailable'
                suffix = ' (screened approximation)' if value == 'coulomb' else ''
                print(f'  {number}. {value:<12} {status}{suffix}')
            try:
                answer = input(f'Choice [{selected[key]}]: ').strip()
            except (EOFError,KeyboardInterrupt) as exc:
                raise ValueError('Selection cancelled; no run was started') from exc
            if answer.isdecimal() and 1 <= int(answer) <= len(options):
                answer = options[int(answer)-1]
            value = next((v for v in options if v.casefold() == answer.casefold()), answer) if answer else selected[key]
            if value not in options:
                raise ValueError(f'Unknown {labels[key].lower()}: {value!r}')
            if not available[key][value]:
                raise ValueError(f'{value} is unavailable; install it or explicitly choose another available option')
            selected[key] = value
    return selected


def main(argv=None):
    command = parser()
    args = command.parse_args(argv)
    if args.command == 'model':
        from .modeling import build_model
        print(f'Optional native modeling: {args.method}; PDBFixer preparation remains the default.', flush=True)
        try:
            report = build_model(args.fasta,args.output,method=args.method,template=args.pdb,
                esmfold_model=args.esmfold_model,device=args.device,python=args.python,revision=args.model_revision)
        except (ValueError,RuntimeError,OSError,ImportError) as exc:
            print(f'Native modeling: FAIL\nERROR: {exc}',file=sys.stderr)
            return 1
        print('\n'.join(report['messages']))
        print(f'Structural candidate: {args.output}/protein.pdb')
        return 0
    if args.command == 'fetch':
        from .fetch import fetch
        try:
            report = fetch(args.output, pdb_id=args.pdb, uniprot_id=args.uniprot)
        except (ValueError, OSError) as exc:
            print(f'Database retrieval: FAIL\nERROR: {exc}', file=sys.stderr)
            return 1
        print('\n'.join(report['messages']))
        print(f'Downloaded records: {args.output}')
        return 0
    if args.command=='doctor':
        from .dependencies import inspect_dependencies,format_report
        try:
            settings=installed_settings()
            project=Path(settings.get('project_root',Path(__file__).resolve().parents[2]))
            report=inspect_dependencies(project,sys.executable,
                **{key:settings.get(key,default) for key,default in
                   {'engine':'openmm','solver':'apbs','force_field':'ff14SB'}.items()})
        except (ValueError,OSError) as exc:
            command.error(str(exc))
        print(json.dumps(report,indent=2) if args.json else format_report(report))
        return 0 if report['capabilities']['runtime']['available'] else 1
    if args.command=='check':
        from .preflight import check_pdb
        try:report=check_pdb(args.pdb,args.json)
        except (ValueError,OSError) as exc:command.error(str(exc))
        print(f"PDB check: {'PASS' if report['ok'] else 'FAIL'}; {report['residue_count']} residues, {report['atom_count']} atoms")
        for issue in report['issues']:print(f"{issue['severity'].upper()}: {issue['message']} {issue['action']}")
        return 0 if report['ok'] else 1
    if args.command=='prepare':
        from .preparation import prepare_pdb
        try:
            report=prepare_pdb(args.pdb,args.output,pH=args.ph,protonation=args.protonation,
                               uniprot_record=args.uniprot_record)
        except (ValueError,RuntimeError,OSError,ImportError) as exc:
            print(f'PDB preparation: FAIL\nERROR: {exc}',file=sys.stderr)
            return 1
        print('\n'.join(report['messages']))
        print(f'Prepared PDB: {report["pdb"]}')
        return 0
    if args.command == 'parameterize':
        from .parameter_serial import run_serial
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
        try:
            result = run_serial(args.project, args.names)
        except (ValueError, RuntimeError, OSError, ImportError) as exc:
            logging.error('%s', exc)
            return 1
        logging.info('Parameterization complete: %s', ', '.join(result['completed']))
        return 0
    if args.topology is not None and args.coordinates is None:
        command.error('--coordinates is required with --topology')
    if args.pdb is not None and args.coordinates is not None:
        command.error('--coordinates is only valid with --topology')
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING,
                        format='%(levelname)s: %(message)s')
    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.WARNING)
    logging.getLogger('cationize_protein').setLevel(logging.DEBUG if args.verbose else logging.INFO)
    try:
        from .dependencies import find_executable
        settings=installed_settings()
        project=Path(settings.get('project_root',Path(__file__).resolve().parents[2]))
        choices=run_choices(args,settings,project)
        field,engine,solver=choices['force_field'],choices['engine'],choices['solver']
        LOG.info('Run settings: engine=%s; solver=%s; force field=%s; Na probe radius=%.4g A',
                 engine,solver,field,args.probe_radius)
        binary=(args.apbs if solver=='apbs' else args.pbsa) if solver!='coulomb' else 'not-used'
        if solver!='coulomb':
            binary=binary or find_executable(solver,project)
            saved=settings.get('executables',{}).get(solver)
            if not binary and saved and Path(saved).is_file():binary=saved
            if not binary:raise ValueError(f'{solver.upper()} was not found; run install.py or provide --{solver} /path/to/{solver}')
        formats=tuple(dict.fromkeys(args.formats+(['amber'] if args.export_amber else [])))
        apbs = APBSConfig(executable=binary, temperature_K=args.temperature,
                          solute_dielectric=args.solute_dielectric, solvent_dielectric=args.solvent_dielectric,
                          salt_M=args.salt, ion_radius_A=args.pb_ion_radius, spacing_A=args.pb_spacing)
        config = WorkflowConfig(output_dir=args.output, parameters_dir=args.parameters,
                                apbs=apbs, probe_radius_A=args.probe_radius,
                                topology_path=args.topology, coordinates_path=args.coordinates,
                                pdb_path=args.pdb, percentages=tuple(args.percentages),
                                sampling_spacing_A=args.sampling_spacing, keep_maps=args.keep_maps,
                                platform=args.platform, beam_width=args.beam_width,
                                max_iterations=args.max_iterations, restraint_kcal_A2=args.restraint,
                                tolerance_kcal_A=args.tolerance,export_amber=args.export_amber,
                                force_field=field,engine=engine,engine_executable=args.engine_executable,
                                solver=solver,pbsa_nonlinear=args.pbsa_nonlinear,
                                protonation=args.protonation,pH=args.ph,formats=formats)
        result = run(config)
    except (ValueError, RuntimeError, OSError, ImportError) as exc:
        logging.error('%s', exc)
        return 1
    LOG.info('Complete: %d/%d groups modified; stages %s', result['final_modified_groups'],
                 result['total_carboxyl_groups'], ', '.join(result['stages']))
    return 0
