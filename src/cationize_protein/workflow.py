"""Serial progressive amidation with explicit physical and protonation choices.

Each committed step contains the current structure, complete parameters,
PROPKA history and checked exports. Resume never changes a committed artifact.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
import fcntl
import hashlib
import json
import logging
import importlib.metadata
import importlib.util
import numbers
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np

from . import builder, catalog, parameter_assignment
from .electrostatics import APBSConfig, build_candidate_grid, rank_candidates, run_apbs
from .potential_solvers import solve
from .sterics import RotatableBond, choose_rotamer, derive_rotations
from .topology import LEGACY_ALIASES, _REUSED_MAPS, inspect_structure

LOG = logging.getLogger(__name__)
_VALIDATED = {'validated', 'reused_validated_HSA_template',
              'RESP_capped_Amber_and_OpenMM_validated'}


class WorkflowError(RuntimeError):
    """A run cannot proceed without changing or repairing an explicit input."""


@dataclass(frozen=True)
class WorkflowConfig:
    output_dir: Path
    parameters_dir: Path
    apbs: APBSConfig
    probe_radius_A: float
    topology_path: Path | None = None
    coordinates_path: Path | None = None
    pdb_path: Path | None = None
    percentages: tuple = (25, 50, 75, 100)
    sampling_spacing_A: float = 1.0
    keep_maps: bool = False
    platform: str = 'CPU'
    restraint_kcal_A2: float = 10.0
    tolerance_kcal_A: float = 0.1
    max_iterations: int = 20000
    beam_width: int = 8
    export_amber: bool = False
    force_field: str = 'ff14SB'
    engine: str = 'openmm'
    engine_executable: str | None = None
    solver: str = 'apbs'
    pbsa_nonlinear: bool = False
    protonation: str = 'propka'
    pH: float = 7.4
    formats: tuple = ('amber',)


def stage_counts(total_groups: int, percentages=(25, 50, 75, 100)):
    """Sorted output labels/counts, including baseline; decimal ROUND_HALF_UP."""
    if isinstance(total_groups, bool) or not isinstance(total_groups, numbers.Integral) or total_groups < 0:
        raise ValueError('total_groups must be a nonnegative integer')
    parsed = []
    try:
        for value in percentages:
            percentage = Decimal(str(value))
            if not percentage.is_finite() or not 0 <= percentage <= 100 or percentage in parsed:
                raise ValueError('Percentages must be unique finite values between 0 and 100')
            parsed.append(percentage)
    except InvalidOperation as exc:
        raise ValueError('Invalid percentage') from exc
    parsed = sorted(set(parsed) | {Decimal(0)})
    return [(format(p.normalize(), 'f'), int((Decimal(int(total_groups)) * p / 100).quantize(Decimal(1), rounding=ROUND_HALF_UP))) for p in parsed]


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path, value):
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError(path)
    temporary = path.with_name('.' + path.name + '.tmp')
    with temporary.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
    temporary.replace(path)


def verify_artifacts(directory, hashes):
    """Verify committed relative paths without following any symlink component."""
    directory = Path(directory)
    for name, expected in hashes.items():
        path = directory / name
        relative=Path(name)
        unsafe=relative.is_absolute() or not relative.parts or any(p in {'.','..'} for p in relative.parts)
        unsafe=unsafe or directory.is_symlink() or any(directory.joinpath(*relative.parts[:i]).is_symlink()
                                                       for i in range(1,len(relative.parts)+1))
        if unsafe or not path.is_file() or _sha(path) != expected:
            raise WorkflowError(f'Checkpoint artifact is missing, unsafe, or changed: {path}')


def validate_inputs(config: WorkflowConfig):
    """Read-only validation/fingerprint, also guarding automatic run reuse."""
    paired = config.topology_path is not None and config.coordinates_path is not None
    if (config.topology_path is None) != (config.coordinates_path is None) or paired == (config.pdb_path is not None):
        raise ValueError('Supply topology plus coordinates, or one native PDB')
    for name in ('probe_radius_A', 'sampling_spacing_A', 'restraint_kcal_A2', 'tolerance_kcal_A'):
        value = getattr(config, name)
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be finite and positive')
    if not isinstance(config.max_iterations, int) or config.max_iterations < 1 or not isinstance(config.beam_width, int) or config.beam_width < 1:
        raise ValueError('Iteration count and beam width must be positive integers')
    labels = [label for label, _ in stage_counts(0, config.percentages)]
    parameters = Path(config.parameters_dir).resolve()
    if not (parameters / 'manifest.json').is_file():
        raise FileNotFoundError(parameters / 'manifest.json')
    inputs = {}
    for name, value in [('topology', config.topology_path), ('coordinates', config.coordinates_path), ('pdb', config.pdb_path)]:
        if value is not None:
            path = Path(value).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            inputs[name + '_path'], inputs[name + '_sha256'] = str(path), _sha(path)
    from .forcefields import profile
    profile(config.force_field)
    if config.engine not in {'openmm','amber','sander','pmemd'} or config.solver not in {'apbs','pbsa','coulomb'}:
        raise ValueError('Unknown minimization engine or electrostatic solver')
    if config.solver=='pbsa' and config.apbs.spacing_A>1.4:
        raise ValueError('PBSA ipb=2 requires grid spacing no larger than its 1.4 A solvent probe')
    if config.protonation not in {'propka','fixed'} or not np.isfinite(config.pH) or not 0<=config.pH<=14:
        raise ValueError('Choose propka/fixed protonation and a pH between 0 and 14')
    if not config.formats or set(config.formats)-{'amber','openmm','gromacs'}:
        raise ValueError('Choose amber, openmm and/or gromacs export formats')
    if config.protonation=='propka' and importlib.util.find_spec('propka') is None:
        raise WorkflowError('PROPKA is required; run the private installer or explicitly choose --protonation fixed')
    executable = shutil.which(str(config.apbs.executable)) if config.solver!='coulomb' else None
    if config.solver!='coulomb' and executable is None:
        raise WorkflowError(f'{config.solver.upper()} executable is unavailable: {config.apbs.executable}')
    settings = asdict(config.apbs)
    settings['executable'] = str(Path(executable).resolve()) if executable else None
    engine_binary=None
    if config.engine!='openmm':
        from .minimization import amber_executable
        engine_binary=amber_executable(config.engine,config.engine_executable)
    signature = dict(schema_version=2, inputs=inputs, parameters_dir=str(parameters),
                     parameter_bundle_sha256={p.name:_sha(p) for p in sorted(parameters.iterdir())
                         if p.is_file() and p.suffix in {'.lib','.json','.xml','.csv'}},
                     native_forcefield_sha256=_native_field_hashes(config.force_field),
                     apbs=settings, apbs_sha256=_sha(executable) if executable else None, percentages=labels,
                     probe_radius_A=float(config.probe_radius_A), sampling_spacing_A=float(config.sampling_spacing_A),
                     keep_maps=config.keep_maps, platform=config.platform,
                     restraint_kcal_A2=config.restraint_kcal_A2,
                     tolerance_kcal_A=config.tolerance_kcal_A, max_iterations=config.max_iterations,
                     beam_width=config.beam_width, export_amber=config.export_amber,
                     force_field=config.force_field,engine=config.engine,solver=config.solver,pbsa_nonlinear=config.pbsa_nonlinear,
                     engine_executable=engine_binary,engine_sha256=_sha(engine_binary) if engine_binary else None,
                     protonation=config.protonation,pH=config.pH,formats=list(config.formats),
                     dependency_versions={name: _version(name) for name in ('openmm','parmed','propka','numpy','scipy')},
                     implementation_sha256={name: _sha(Path(__file__).with_name(name)) for name in
                                            ('workflow.py', 'topology.py', 'electrostatics.py', 'sterics.py',
                                             'builder.py', 'openmm_backend.py', 'catalog.py', 'off.py',
                                             'forcefields.py','parameter_assignment.py','minimization.py','exports.py','amber_bundle.py',
                                             'protonation.py','native_states.py','preflight.py','potential_solvers.py','gromacs_export.py')})
    output = Path(config.output_dir)
    if output.is_symlink():
        raise WorkflowError('Output directories must not be symlinks')
    if output.exists():
        marker = output / 'RUN.json'
        if not output.is_dir() or not marker.is_file() or marker.is_symlink():
            raise WorkflowError('Existing output is not a recognized committed run; use a new directory')
        if json.loads(marker.read_text()) != signature:
            raise WorkflowError('Existing run inputs or settings differ; use a new output directory')
    return signature


def _version(name):
    try:return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:return 'not installed'


def _native_field_hashes(force_field):
    from .forcefields import data_directory,profile
    try:directory=data_directory()
    except ImportError:return {'runtime':'OpenMM not installed'}
    return {name:_sha(directory/name) for name in (profile(force_field)['xml'],'implicit/gbn2.xml')}


def require_validated_template(parameters_dir, candidate):
    """Validate exactly the selected target; pending candidates are never skipped."""
    directory = Path(parameters_dir)
    manifest = json.loads((directory / 'manifest.json').read_text())
    name = candidate.target_spec['name']
    entries = [entry for entry in manifest.get('templates', []) if entry.get('name') == name]
    entry = entries[0] if len(entries) == 1 else {}
    status = entry.get('status', 'missing or ambiguous manifest entry')
    if status not in _VALIDATED:
        raise WorkflowError(f'Selected site {candidate.candidate_id} requires {name}; template status: {status}')
    for field in ('parent', 'context', 'modified_groups'):
        if entry.get(field) != candidate.target_spec.get(field):
            raise WorkflowError(f'Selected site {candidate.candidate_id}: {name} manifest chemistry differs')
    library = directory / (name + '.lib')
    if not library.is_file() or library.is_symlink() or not entry.get('library_sha256') or _sha(library) != entry['library_sha256']:
        raise WorkflowError(f'Selected site {candidate.candidate_id}: {name} library missing or hash mismatch')
    return library


def _used_templates(info, parameters_dir):
    names = sorted({residue.source_spec['name'] for residue in info.residues if residue.source_spec})
    for name in names:
        require_validated_template(parameters_dir, SimpleNamespace(candidate_id='existing:' + name, target_spec=catalog.by_name(name)))
    return names


def _normalize_labels(structure):
    disulfide = {id(atom.residue) for bond in structure.bonds
                 if bond.atom1.name == bond.atom2.name == 'SG'
                 for atom in (bond.atom1, bond.atom2)}
    for residue in structure.residues:
        residue.name = LEGACY_ALIASES.get(residue.name, residue.name)
        if residue.name == 'HIS':
            hydrogens = {atom.name for atom in residue.atoms} & {'HD1', 'HE2'}
            if not hydrogens:
                raise WorkflowError('HIS needs explicit HD1/HE2 hydrogens or a declared HID/HIE/HIP state')
            residue.name = 'HIP' if len(hydrogens) == 2 else 'HID' if 'HD1' in hydrogens else 'HIE'
        if residue.name == 'CYS' and id(residue) in disulfide:
            residue.name = 'CYX'
    return structure


def _pdb_structure(path, parameters, force_field='ff14SB', pH=7.4):
    raw = []
    last = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            if line[16:17].strip():
                raise WorkflowError('PDB alternate locations must be resolved before use')
            key = line[21:27]
            if key != last:
                original = line[17:20].strip()
                supported = {spec['parent'] for spec in catalog.templates()} | {'HIS'}
                if original not in supported or original in {'CYM'}:
                    raise WorkflowError(f'Unsupported native PDB state {original}; use a prepared parameterized topology')
                raw.append(original)
                last = key
    import parmed
    from openmm.app import Modeller, PDBFile
    from .forcefields import load_forcefield
    pdb = PDBFile(str(path))
    from .preflight import check_pdb
    _verify_declared_disulfides(pdb,check_pdb(path))
    residues = list(pdb.topology.residues())
    if len(raw) != len(residues):
        raise WorkflowError('PDB residue identities are ambiguous; provide one resolved model')
    variants = []
    for residue, original in zip(residues, raw):
        names = {atom.name for atom in residue.atoms()}
        if original in {'HID', 'HIE', 'HIP'}:
            if (original == 'HID' and 'HE2' in names) or (original == 'HIE' and 'HD1' in names):
                raise WorkflowError('Declared histidine state contradicts its explicit ring hydrogens')
            variants.append(original)
        elif residue.name == 'HIS':
            present = names & {'HD1', 'HE2'}
            if not present:
                raise WorkflowError('Ambiguous HIS protonation; provide HID/HIE/HIP labels or explicit ring hydrogens')
            variants.append('HIP' if len(present) == 2 else 'HID' if 'HD1' in present else 'HIE')
        else:
            if original == 'CYX':
                sulfur = next((atom for atom in residue.atoms() if atom.name == 'SG'), None)
                if sulfur is None or not any(sulfur in (a, b) and a.name == b.name == 'SG' for a, b in pdb.topology.bonds()):
                    raise WorkflowError('A declared CYX PDB residue has no SG-SG bond')
            variants.append(None)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(load_forcefield(parameters,force_field,names=[]),pH=pH,variants=variants)
    structure = parmed.openmm.load_topology(modeller.topology, xyz=modeller.positions)
    for residue, variant in zip(structure.residues, variants):
        if variant:
            residue.name = variant
    structure = _normalize_labels(structure)
    typed, system = parameter_assignment.parameterize(structure, parameters,force_field=force_field,names=[])
    return typed, system


def _verify_declared_disulfides(pdb,report):
    observed={tuple(sorted((a.residue.index+1,b.residue.index+1))) for a,b in pdb.topology.bonds()
              if a.name==b.name=='SG'}
    for declaration in report.get('ssbond_declarations',[]):
        pair=tuple(sorted(declaration.get('residue_ordinals',[])))
        if len(pair)!=2 or pair not in observed:
            raise WorkflowError(f'Declared disulfide {pair} is absent from the prepared graph; review its SG coordinates/connectivity')


def _initial_structure(config):
    import parmed
    from openmm.app import NoCutoff
    if config.pdb_path is not None:
        return _pdb_structure(config.pdb_path, config.parameters_dir,config.force_field,config.pH)
    structure = _normalize_labels(parmed.load_file(str(config.topology_path), xyz=str(config.coordinates_path)))
    info = inspect_structure(structure)
    _used_templates(info, config.parameters_dir)
    return parameter_assignment.parameterize(structure,config.parameters_dir,force_field=config.force_field)


def _native_reference(info):
    return {f'{residue.ordinal}|{name}': info.coordinates[index].tolist()
            for residue in info.residues for name, index in residue.native_atom_indices.items()
            if name in {'N', 'CA', 'C', 'O', 'SG'} and index in info.protected_atoms}


def _reference_positions(info, reference):
    positions = info.coordinates.copy()
    for residue in info.residues:
        for name, index in residue.native_atom_indices.items():
            key = f'{residue.ordinal}|{name}'
            if key in reference:
                positions[index] = reference[key]
    return positions


def inherited_atom_map(structure, residue_state, template, target_spec):
    """Map retained heavy atoms and every unchanged attached hydrogen explicitly."""
    old = {atom.name: atom for atom in structure.residues[residue_state.ordinal - 1].atoms}
    product = {atom.name: atom for atom in template.atoms}
    roles = _REUSED_MAPS.get(target_spec.get('reused_from'), {})
    mapping = {}
    for native, index in residue_state.native_atom_indices.items():
        atom = structure.atoms[index]
        target = roles.get(native, native)
        if atom.atomic_number != 1 and target in product:
            mapping[atom.name] = target
    old_graph = {name: set() for name in old}
    new_graph = {name: set() for name in product}
    for bond in structure.bonds:
        if bond.atom1.name in old and bond.atom2.name in old and bond.atom1.residue is structure.residues[residue_state.ordinal - 1] and bond.atom2.residue is bond.atom1.residue:
            old_graph[bond.atom1.name].add(bond.atom2.name)
            old_graph[bond.atom2.name].add(bond.atom1.name)
    for bond in template.bonds:
        new_graph[bond.atom1.name].add(bond.atom2.name)
        new_graph[bond.atom2.name].add(bond.atom1.name)
    for source, target in list(mapping.items()):
        before = sorted(name for name in old_graph[source] if old[name].atomic_number == 1)
        after = sorted(name for name in new_graph[target] if product[name].atomic_number == 1)
        if len(before) != len(after):
            raise WorkflowError(f'Inherited atom {source} changes attached hydrogen count')
        same = sorted(set(before) & set(after))
        mapping.update((name, name) for name in same)
        mapping.update(zip([name for name in before if name not in same], [name for name in after if name not in same]))
    return mapping


def _rotatable_bonds(structure, info, ordinal):
    members = set(info.residues[ordinal - 1].atom_indices)
    graph = [set() for _ in info.atoms]
    for i, j in info.bonds:
        graph[i].add(j)
        graph[j].add(i)
    residue_ids = [str(atom.residue_ordinal) for atom in info.atoms]
    axes = []
    for bond in structure.bonds:
        i, j = bond.atom1.idx, bond.atom2.idx
        if i not in members or j not in members or min(bond.atom1.atomic_number, bond.atom2.atomic_number) == 1 or float(bond.order) > 1.1:
            continue
        # Preserve amide and guanidinium/amidinium resonance when OFF has order 1.
        if {bond.atom1.atomic_number, bond.atom2.atomic_number} == {6, 7}:
            carbon = i if bond.atom1.atomic_number == 6 else j
            if (any(info.atoms[k].atomic_number == 8 for k in graph[carbon])
                    or sum(info.atoms[k].atomic_number == 7 for k in graph[carbon]) >= 2):
                continue
        for base, tip in ((i, j), (j, i)):
            axis = RotatableBond(base, tip)
            try:
                derive_rotations(info.bonds, [axis], residue_ids,
                                 protected_atoms=info.protected_atoms,
                                 disulfide_bonds=info.disulfide_bonds)
            except ValueError:
                continue
            axes.append(axis)
            break
    return sorted(axes, key=lambda axis: (axis.base, axis.tip))


def _geometry(info):
    distances = [np.linalg.norm(info.coordinates[i] - info.coordinates[j]) for i, j in info.bonds]
    if distances and (min(distances) < 0.6 or max(distances) > 2.5):
        raise WorkflowError('Invalid covalent geometry after placement/minimization')
    return dict(min_bond_A=float(min(distances)) if distances else None,
                max_bond_A=float(max(distances)) if distances else None)


def _export_amber(structure, directory):
    """Optional checked interchange files; XML/JSON remains authoritative."""
    import parmed
    top, rst = directory / 'protein.parm7', directory / 'protein.rst7'
    try:
        amber = parmed.amber.AmberParm.from_structure(structure)
        amber.save(str(top), overwrite=False)
        amber.save(str(rst), format='rst7', overwrite=False)
        restored = parmed.load_file(str(top), xyz=str(rst))
        before, after = inspect_structure(structure), inspect_structure(restored)
        if before.bonds != after.bonds or not np.allclose(before.coordinates, after.coordinates, atol=1e-6, rtol=0) or not np.allclose(before.charges_e, after.charges_e, atol=1e-7, rtol=0) or not np.allclose(before.lj_radii_A, after.lj_radii_A, atol=1e-7, rtol=0):
            raise WorkflowError('ParmEd Amber round trip changed physical data')
        return dict(validated=True, files=[top.name, rst.name])
    except Exception as exc:
        for path in (top, rst):
            if path.exists():
                path.rename(path.with_name('UNVALIDATED_' + path.name))
        return dict(validated=False, error=str(exc), authority='Use system.xml and state.json')


def _save_checkpoint(structure, system, directory, metadata, *, export_amber=False,
                     formats=None,force_field='ff14SB',audit_directory=None,preflight=None):
    from openmm import XmlSerializer
    directory.mkdir(parents=True, exist_ok=False)
    info = inspect_structure(structure)
    pdb, graph = builder.export_mutation_pdb(structure)
    (directory / 'protein.pdb').write_text(pdb)
    (directory / 'system.xml').write_text(XmlSerializer.serialize(system))
    graph['requires_parameterization'] = False
    _json(directory / 'state.json', dict(graph=graph, masses=[float(a.mass) for a in structure.atoms],
                                         pb_radii_A=info.pb_radii_A.tolist(), lj_radii_A=info.lj_radii_A.tolist()))
    files = ['protein.pdb', 'system.xml', 'state.json']
    record = dict(metadata)
    if formats is not None or export_amber:
        from .exports import export_system
        physical,_=parameter_assignment.parameterize(structure,metadata['parameters_dir'],force_field=force_field)
        record['physical_exports']=export_system(physical,directory/'exports',
                    formats=formats or ('amber',),force_field=force_field)
        files.extend(str(p.relative_to(directory)) for p in sorted((directory/'exports').rglob('*')) if p.is_file())
    if audit_directory is not None and Path(audit_directory).is_dir():
        shutil.copytree(audit_directory,directory/'protonation')
        files.extend(str(p.relative_to(directory)) for p in sorted((directory/'protonation').rglob('*')) if p.is_file())
    if preflight is not None:
        _json(directory/'PDB_CHECK.json',preflight);files.append('PDB_CHECK.json')
    record['artifacts'] = {name: _sha(directory / name) for name in files}
    _json(directory / 'STEP_COMPLETE.json', record)
    return record


def _read_checkpoint(directory):
    import parmed
    from openmm import XmlSerializer, unit
    directory = Path(directory)
    marker = json.loads((directory / 'STEP_COMPLETE.json').read_text())
    verify_artifacts(directory, marker['artifacts'])
    state = json.loads((directory / 'state.json').read_text())
    graph = state['graph']
    untyped = parmed.Structure()
    for source in graph['residues']:
        residue = parmed.Residue(source['name'], source['source_number'], chain=source['chain'],
                                 insertion_code=source['insertion_code'], segid=source['segid'],
                                 list=untyped.residues)
        untyped.residues.append(residue)
    for atom, mass in zip(graph['atoms'], state['masses']):
        created = parmed.Atom(name=atom['name'], atomic_number=atom['atomic_number'], mass=mass)
        untyped.residues[atom['residue']].add_atom(created)
        untyped.atoms.append(created)
    for i, j, order in graph['bonds']:
        untyped.bonds.append(parmed.Bond(untyped.atoms[i], untyped.atoms[j], order=order))
    system = XmlSerializer.deserialize((directory / 'system.xml').read_text())
    structure = parmed.openmm.load_topology(untyped.topology, system,
                                           xyz=np.array(graph['coordinates_A']) * unit.angstrom)
    for source, residue in zip(graph['residues'], structure.residues):
        residue.name, residue.number, residue.chain = source['name'], source['source_number'], source['chain']
        residue.insertion_code, residue.segid, residue.ter = source['insertion_code'], source['segid'], source['declared_ter']
    for atom, radius in zip(structure.atoms, state['pb_radii_A']):
        atom.solvent_radius = radius
    info = inspect_structure(structure)
    if not np.allclose(info.charges_e, [a['charge_e'] for a in graph['atoms']], atol=1e-8, rtol=0) or not np.allclose(info.lj_radii_A, state['lj_radii_A'], atol=1e-8, rtol=0):
        raise WorkflowError('XML/JSON checkpoint physical parameters disagree')
    return structure, system, marker


def _export_stage(root, label, checkpoint, metadata):
    destination = root / label
    if destination.exists():
        marker = destination / 'STAGE_COMPLETE.json'
        if not marker.is_file():
            if destination.is_symlink() or not destination.is_dir():raise WorkflowError('Unsafe incomplete stage export')
            failed=_attempt_directory(root/'internal'/'partial_exports',0)
            destination.rename(failed/str(label))
            return _export_stage(root,label,checkpoint,metadata)
        saved = json.loads(marker.read_text())
        verify_artifacts(destination, saved['artifacts'])
        if saved['checkpoint_sha256'] != _sha(checkpoint / 'STEP_COMPLETE.json'):
            raise WorkflowError('Stage export points to a different checkpoint')
        return saved
    staging=_attempt_directory(root/'internal'/'stage_exports',0)
    source = json.loads((checkpoint / 'STEP_COMPLETE.json').read_text())
    verify_artifacts(checkpoint, source['artifacts'])
    for name in source['artifacts']:
        (staging/name).parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(checkpoint / name, staging / name)
    result = dict(metadata, checkpoint_sha256=_sha(checkpoint / 'STEP_COMPLETE.json'),
                  source_checkpoint=str(checkpoint.relative_to(root)),
                  protein_charge_e=source['protein_charge_e'], checkpoint_metadata=source,
                  artifacts={name: _sha(staging / name) for name in source['artifacts']})
    _json(staging / 'STAGE_COMPLETE.json', result)
    staging.rename(destination)
    LOG.info('Saved %s%% checkpoint: %s',label,destination)
    return result


def _attempt_directory(root, step):
    root.mkdir(parents=True,exist_ok=True)
    index = 1
    while True:
        directory = root / f'{step:06d}_attempt{index:03d}'
        try:
            directory.mkdir()
            return directory
        except FileExistsError:
            index += 1


def _minimize(structure,config,reference_positions,directory=None):
    from .minimization import minimize
    return minimize(structure,config.parameters_dir,force_field=config.force_field,engine=config.engine,
        executable=config.engine_executable,directory=directory,reference_positions=reference_positions,
        restraint_kcal_A2=config.restraint_kcal_A2,tolerance_kcal_A=config.tolerance_kcal_A,
        max_iterations=config.max_iterations,platform=config.platform)


def _protonate(structure,directory,config,reference,previous=None):
    """Audit one chemical step; rebuild only supported persistent native changes."""
    if config.protonation=='fixed':return structure,dict(mode='fixed',history={})
    LOG.info('Auditing native protonation at pH %.2f; introduced ammonium groups remain +1',config.pH)
    from .protonation import audit
    from .native_states import apply_native_states,supported_native_contexts
    directory=Path(directory)
    contexts=supported_native_contexts(config.force_field)
    record=audit(structure,directory/'before',previous_proposals=previous,pH=config.pH,
                 supported_native_states=set(contexts))
    info=inspect_structure(structure)
    changes=dict(record['applicable_native_states'])
    for ordinal,state in list(changes.items()):
        if info.residues[int(ordinal)-1].context not in contexts[state]:
            del changes[ordinal]
            for decision in record['decisions']:
                if str(decision['residue_ordinal'])==ordinal:
                    decision.update(status='held_unsupported_context',reason=f'{config.force_field} has no matching native terminal/context template for {state}')
            record['has_unsupported_proposals']=True
    record['applicable_native_states']=changes
    result=dict(mode='propka',pH=config.pH,before=record,history=record)
    if changes:
        LOG.info('PROPKA: rebuilding native states %s, then minimizing again',changes)
        structure,_,applied=apply_native_states(structure,changes,config.parameters_dir,force_field=config.force_field)
        current=inspect_structure(structure)
        structure,minimum=_minimize(structure,config,_reference_positions(current,reference),directory/'rebuild_minimum')
        # A post-minimum audit checks the new structure but cannot constitute
        # a second chemical step for the persistence criterion.
        after=audit(structure,directory/'after',previous_proposals=None,pH=config.pH,supported_native_states=set(contexts))
        result.update(state_changes=applied,rebuild_minimization=minimum,after=after,history=after)
    held=[d for d in result['history']['decisions'] if d['status'].startswith('held_unsupported') or d['status']=='held_modified_template']
    if held:LOG.warning('PROPKA: %d proposals retained because their fitted/native template is unavailable; see protonation/POLICY.json',len(held))
    directory.mkdir(parents=True,exist_ok=True)
    _json(directory/'POLICY.json',result)
    return structure,result


def _mutate(structure, info, candidate, library, config, reference,directory=None):
    import parmed
    templates = parmed.amber.AmberOFFLibrary.parse(str(library))
    template = templates[candidate.target_spec['name']]
    state = info.residues[candidate.residue_ordinal - 1]
    mapping = inherited_atom_map(structure, state, template, candidate.target_spec)
    replaced = builder.replace_residue(structure, state.ordinal - 1, template, mapping)
    names = sorted(set(_used_templates(info, config.parameters_dir)) | {candidate.target_spec['name']})
    # A partial state can disappear when its remaining terminal/side group is changed.
    present = {residue.name for residue in replaced.residues}
    names = [name for name in names if name in present]
    typed, _ = parameter_assignment.parameterize(replaced, config.parameters_dir,force_field=config.force_field,names=names)
    current = inspect_structure(typed)
    placement = _geometry(current)
    axes = _rotatable_bonds(typed, current, state.ordinal)
    steric = dict(rotations_degrees=[], reason='No permitted nonamide heavy-atom rotation')
    if axes:
        typed.coordinates, steric = choose_rotamer(
            current.coordinates, current.bonds, current.lj_radii_A, axes,
            residue_ids=[str(atom.residue_ordinal) for atom in current.atoms],
            protected_atoms=current.protected_atoms, disulfide_bonds=current.disulfide_bonds,
            beam_width=config.beam_width)
    current = inspect_structure(typed)
    minimized, minimization = _minimize(typed,config,_reference_positions(current,reference),directory)
    final, system = parameter_assignment.parameterize(minimized,config.parameters_dir,force_field=config.force_field,names=names)
    final_info = inspect_structure(final)
    if len(final_info.residues) != len(info.residues) or sum(len(r.modified_groups) for r in final_info.residues) != sum(len(r.modified_groups) for r in info.residues) + 1:
        raise WorkflowError('Mutation did not add exactly one independent chemical group')
    old_ss = {tuple(sorted((info.atoms[i].residue_ordinal, info.atoms[j].residue_ordinal))) for i, j in info.disulfide_bonds}
    new_ss = {tuple(sorted((final_info.atoms[i].residue_ordinal, final_info.atoms[j].residue_ordinal))) for i, j in final_info.disulfide_bonds}
    if old_ss != new_ss:
        raise WorkflowError('Mutation changed the native disulfide set')
    for before, after in zip(info.residues, final_info.residues):
        roles = ('N', 'C', 'CB')
        if all(name in before.native_atom_indices and name in after.native_atom_indices for name in (*roles, 'CA')):
            first = np.array([info.coordinates[before.native_atom_indices[name]] - info.coordinates[before.native_atom_indices['CA']] for name in roles])
            second = np.array([final_info.coordinates[after.native_atom_indices[name]] - final_info.coordinates[after.native_atom_indices['CA']] for name in roles])
            if np.linalg.det(first) * np.linalg.det(second) <= 0 or abs(np.linalg.det(second)) < .1:
                raise WorkflowError(f'Mutation inverted or flattened alpha stereochemistry at residue {before.ordinal}')
    return final, system, dict(atom_map=mapping, placement_geometry=placement,
                               final_geometry=_geometry(final_info), sterics=steric, minimization=minimization)


def _one_step(structure, info, step, total_groups, steps_root, config, reference,previous_protonation=None):
    attempt = _attempt_directory(steps_root, step)
    LOG.info('Step %d/%d: evaluating %s potential for %d remaining carboxyl groups',
             step,total_groups,config.solver.upper(),len(info.candidates))
    centers = {candidate.candidate_id: info.coordinates[list(candidate.oxygen_indices)] for candidate in info.candidates}
    grid = build_candidate_grid(centers, info.coordinates, info.lj_radii_A,
                                config.probe_radius_A, spacing_A=config.sampling_spacing_A)
    scores,pb=solve(structure,info,grid,attempt/'pb',config.apbs,solver=config.solver,
        provenance=dict(step=step,total_carboxyl_groups=total_groups),keep_map=config.keep_maps,
        pbsa_nonlinear=config.pbsa_nonlinear)
    ranking, unavailable = rank_candidates(scores, grid.indices, grid.blocked)
    _json(attempt / 'ranking.json', dict(rows=[asdict(row) for row in ranking], unavailable=unavailable))
    if not ranking:
        raise WorkflowError('No accessible unmodified carboxyl group remains')
    by_id = {candidate.candidate_id: candidate for candidate in info.candidates}
    candidate = by_id[ranking[0].candidate_id]
    library = require_validated_template(config.parameters_dir, candidate)
    LOG.info('Step %d: %s -> %s (%.6g kBT/e)', step, candidate.candidate_id,
             candidate.target_spec['name'], ranking[0].score_kBT_per_e)
    final, system, mutation = _mutate(structure,info,candidate,library,config,reference,attempt/'minimum')
    final,protonation=_protonate(final,attempt/'protonation',config,reference,previous_protonation)
    final,system=parameter_assignment.parameterize(final,config.parameters_dir,force_field=config.force_field)
    updated = inspect_structure(final)
    modified = sum(len(residue.modified_groups) for residue in updated.residues)
    if modified + len(updated.candidates) != total_groups:
        raise WorkflowError('Total independently targetable group count changed')
    checkpoint = attempt / 'checkpoint'
    marker = _save_checkpoint(final, system, checkpoint,
                              dict(step=step, modified_groups=modified, selected=asdict(candidate),
                                   score=asdict(ranking[0]), library_path=str(library.resolve()),
                                   library_sha256=_sha(library), mutation=mutation, PB=pb,
                                   field_artifacts={path.name: _sha(path) for path in (attempt / 'pb').iterdir() if path.is_file()},
                                   ranking_sha256=_sha(attempt / 'ranking.json'),
                                   protein_charge_e=float(updated.charges_e.sum()),protonation=protonation,
                                   parameters_dir=str(Path(config.parameters_dir).resolve()),force_field=config.force_field),
                              formats=config.formats,force_field=config.force_field,audit_directory=attempt/'protonation')
    return final, system, checkpoint, marker


def run(config: WorkflowConfig):
    """Create or automatically resume one hash-verified progressive run."""
    signature = validate_inputs(config)
    pdb_report=None
    if config.pdb_path is not None:
        from .preflight import check_pdb
        pdb_report=check_pdb(config.pdb_path)
        if not pdb_report['ok']:
            raise WorkflowError('PDB check failed: '+'; '.join(i['message'] for i in pdb_report['issues'] if i['severity']=='error'))
    root = Path(config.output_dir).resolve()
    new = not root.exists()
    if new:
        structure, system = _initial_structure(config)
        initial = inspect_structure(structure)
        modified = sum(len(residue.modified_groups) for residue in initial.residues)
        total = modified + len(initial.candidates)
        plan = stage_counts(total, config.percentages)
        if any(label != '0' and count < modified for label, count in plan):
            raise WorkflowError('A requested percentage is below the already modified input group count')
        reference = _native_reference(initial)
        root.mkdir(parents=True)
    with (root / '.run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise WorkflowError('This output directory already has an active workflow') from exc
        internal = root / 'internal'
        baseline = internal / 'baseline'
        if new:
            _json(root / 'RUN.json', signature)
            _json(root / 'REFERENCE.json', dict(positions=reference, total_groups=total,
                                               initial_modified_groups=modified, stages=plan))
        if not (baseline/'STEP_COMPLETE.json').is_file():
            if not new:
                structure,system=_initial_structure(config)
                initial=inspect_structure(structure)
                modified=sum(len(r.modified_groups) for r in initial.residues)
                total=modified+len(initial.candidates);plan=stage_counts(total,config.percentages)
                reference=_native_reference(initial)
                expected=dict(positions=reference,total_groups=total,initial_modified_groups=modified,stages=plan)
                if (root/'REFERENCE.json').exists():
                    # JSON turns tuple plan entries into lists.
                    if json.loads((root/'REFERENCE.json').read_text())!=json.loads(json.dumps(expected)):
                        raise WorkflowError('Incomplete baseline source reference changed')
                else:_json(root/'REFERENCE.json',expected)
            if baseline.exists():
                if baseline.is_symlink() or not baseline.is_dir():raise WorkflowError('Unsafe incomplete baseline')
                failed=_attempt_directory(internal/'partial_baselines',0)
                baseline.rename(failed/'checkpoint')
            attempt=_attempt_directory(internal/'baseline_attempts',0)
            baseline_audit=attempt/'protonation'
            baseline_minimum=None
            if config.protonation=='propka':
                LOG.info('Preparing baseline: %s restrained minimization with %s',config.engine,config.force_field)
                structure,baseline_minimum=_minimize(structure,config,_reference_positions(initial,reference),attempt/'minimum')
            structure,protonation=_protonate(structure,baseline_audit,config,reference)
            structure,system=parameter_assignment.parameterize(structure,config.parameters_dir,force_field=config.force_field) if config.protonation=='propka' else (structure,system)
            initial=inspect_structure(structure)
            _save_checkpoint(structure, system, baseline, dict(step=0, modified_groups=modified,
                             reference_sha256=_sha(root / 'REFERENCE.json'),
                             template_hashes={name: _sha(Path(config.parameters_dir) / (name + '.lib')) for name in _used_templates(initial, config.parameters_dir)},
                             protein_charge_e=float(initial.charges_e.sum()),protonation=protonation,
                             baseline_minimization=baseline_minimum,force_field=config.force_field,
                             parameters_dir=str(Path(config.parameters_dir).resolve())),formats=config.formats,
                             force_field=config.force_field,audit_directory=baseline_audit,preflight=pdb_report)
        if not (baseline / 'STEP_COMPLETE.json').is_file():
            raise WorkflowError('Baseline is incomplete; preserved files require a new output directory')
        baseline_record = json.loads((baseline / 'STEP_COMPLETE.json').read_text())
        if _sha(root / 'REFERENCE.json') != baseline_record['reference_sha256']:
            raise WorkflowError('Native reference positions or percentage plan changed')
        reference_record = json.loads((root / 'REFERENCE.json').read_text())
        reference, total = reference_record['positions'], reference_record['total_groups']
        plan = [(label, int(count)) for label, count in reference_record['stages']]
        checkpoints = [baseline]
        steps_root = internal / 'steps'
        steps_root.mkdir(exist_ok=True)
        checkpoints.extend(path.parent for path in sorted(steps_root.glob('*/checkpoint/STEP_COMPLETE.json')))
        committed = {}
        for checkpoint in checkpoints:
            record = json.loads((checkpoint / 'STEP_COMPLETE.json').read_text())
            verify_artifacts(checkpoint, record['artifacts'])
            step = record['step']
            if step in committed:
                raise WorkflowError('Duplicate committed step; output history is ambiguous')
            for name, digest in record.get('template_hashes', {}).items():
                if _sha(Path(config.parameters_dir) / (name + '.lib')) != digest:
                    raise WorkflowError('A baseline residue library changed')
            if step:
                verify_artifacts(checkpoint.parent / 'pb', record['field_artifacts'])
                if _sha(checkpoint.parent / 'ranking.json') != record['ranking_sha256']:
                    raise WorkflowError('Committed APBS ranking changed')
                library = Path(record['library_path'])
                if not library.is_file() or _sha(library) != record['library_sha256']:
                    raise WorkflowError('A previously used residue library changed')
            committed[step] = (checkpoint, record)
        if sorted(committed) != list(range(max(committed) + 1)):
            raise WorkflowError('Committed step history has a gap')
        # Reuse each historical checkpoint for the stage with its exact count.
        initial_count = reference_record['initial_modified_groups']
        maximum = max(count for _, count in plan)
        step = max(committed)
        latest, latest_record = committed[step]
        structure, system, _ = _read_checkpoint(latest)
        _used_templates(inspect_structure(structure), config.parameters_dir)
        for label, target in plan:
            target_step = 0 if label == '0' else target - initial_count
            if target_step in committed:
                checkpoint, record = committed[target_step]
                _export_stage(root, label, checkpoint, dict(requested_percent=label,
                              modified_groups=record['modified_groups'], total_carboxyl_groups=total, step=target_step))
        while initial_count + step < maximum:
            info = inspect_structure(structure)
            structure, system, latest, latest_record = _one_step(
                structure,info,step+1,total,steps_root,config,reference,
                previous_protonation=latest_record.get('protonation',{}).get('history',{}))
            step += 1
            committed[step] = (latest, latest_record)
            # Publish reached stages immediately, even if a later template is pending.
            for label, target in plan:
                if label != '0' and target == initial_count + step:
                    _export_stage(root, label, latest, dict(requested_percent=label, modified_groups=target,
                                  total_carboxyl_groups=total, step=step))
        for label, target in plan:
            target_step = 0 if label == '0' else target - initial_count
            checkpoint, record = committed[target_step]
            _export_stage(root, label, checkpoint, dict(requested_percent=label,
                          modified_groups=record['modified_groups'], total_carboxyl_groups=total, step=target_step))
        completed = dict(total_carboxyl_groups=total, final_modified_groups=initial_count + step,
                         stages=[label for label, _ in plan], last_checkpoint=str(latest.relative_to(root)))
        if not (root / 'RUN_COMPLETE.json').exists():
            _json(root / 'RUN_COMPLETE.json', completed)
        return completed
