"""Validated physical Amber, OpenMM and GROMACS exports."""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
from .forcefields import profile

_FORMATS = {'amber','openmm','gromacs'}
_PHYSICAL_FORCES = {'HarmonicBondForce','HarmonicAngleForce','PeriodicTorsionForce',
                    'RBTorsionForce','CMAPTorsionForce','NonbondedForce'}


def _physical_system(structure):
    from openmm.app import NoCutoff
    system = getattr(structure, '_cationize_system', None)
    if system is None:
        if getattr(structure, 'unknown_functional', False):
            raise ValueError('Unknown physical force terms cannot be exported')
        system = structure.createSystem(nonbondedMethod=NoCutoff, constraints=None,
                                        rigidWater=False, removeCMMotion=False)
    rejected = [type(f).__name__ for f in system.getForces() if type(f).__name__ not in _PHYSICAL_FORCES]
    if rejected:
        raise ValueError('Export requires physical vacuum parameters; unsupported forces: '+', '.join(rejected))
    if system.getNumConstraints():
        raise ValueError('Constrained source Systems require explicit physical reparameterization')
    return system


def _nonbonded(system):
    from openmm import unit
    forces = [f for f in system.getForces() if type(f).__name__ == 'NonbondedForce']
    if len(forces) != 1: raise ValueError('Expected exactly one NonbondedForce')
    force = forces[0]
    particles = np.array([[q.value_in_unit(unit.elementary_charge),s.value_in_unit(unit.nanometer),
                            e.value_in_unit(unit.kilojoules_per_mole)] for q,s,e in
                           (force.getParticleParameters(i) for i in range(force.getNumParticles()))])
    exceptions = {}
    for i in range(force.getNumExceptions()):
        a,b,q,s,e = force.getExceptionParameters(i)
        exceptions[tuple(sorted((a,b)))] = np.array([q.value_in_unit(unit.elementary_charge**2),
            s.value_in_unit(unit.nanometer),e.value_in_unit(unit.kilojoules_per_mole)])
    return particles, exceptions


def _single_point(system, xyz):
    from openmm import Context, VerletIntegrator, Platform, unit
    integrator = VerletIntegrator(.001*unit.picoseconds)
    context = Context(system, integrator, Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(np.asarray(xyz)*unit.angstrom)
        state = context.getState(getEnergy=True,getForces=True)
        energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        forces = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole/unit.nanometer))
        if not np.isfinite(energy) or not np.isfinite(forces).all():
            raise ValueError('Non-finite vacuum single-point energy or forces')
        return energy, forces
    finally:
        del context, integrator


def _one_four_scaling(system):
    particles, exceptions = _nonbonded(system)
    coulomb, lj = [], []
    for (a,b), values in exceptions.items():
        qprod = particles[a,0]*particles[b,0]
        epsilon = float(np.sqrt(particles[a,2]*particles[b,2]))
        if abs(values[0]) > 1e-14 and abs(qprod) > 1e-14: coulomb.append(values[0]/qprod)
        if values[2] > 1e-14 and epsilon > 1e-14: lj.append(values[2]/epsilon)
    if not coulomb or not lj:
        raise ValueError('Uniform 1-4 settings cannot be inferred from this source System')
    if np.ptp(coulomb) > 1e-7 or np.ptp(lj) > 1e-7:
        raise ValueError('Uniform 1-4 settings cannot represent nonuniform source scaling')
    return float(np.mean(coulomb)), float(np.mean(lj))


def _validate(source, source_system, reloaded):
    from openmm.app import NoCutoff
    if len(source.atoms) != len(reloaded.atoms): raise ValueError('Export changed atom count')
    graph = lambda s: {tuple(sorted((b.atom1.idx,b.atom2.idx))) for b in s.bonds}
    if graph(source) != graph(reloaded): raise ValueError('Export changed covalent graph')
    if [(a.residue.name,a.name,a.atomic_number) for a in source.atoms] != [(a.residue.name,a.name,a.atomic_number) for a in reloaded.atoms]:
        raise ValueError('Export changed atom identity or order')
    rebuilt = reloaded.createSystem(nonbondedMethod=NoCutoff, constraints=None, rigidWater=False,removeCMMotion=False)
    qlj, exceptions = _nonbonded(source_system); actual, actual_exceptions = _nonbonded(rebuilt)
    particle_delta=np.abs(qlj-actual)
    # OpenMM may use sigma=1 nm for zero-epsilon OH hydrogens. Amber
    # legitimately serializes their A/B coefficients as zero. Sigma then has
    # no LJ contribution under either model; q and epsilon still must match.
    no_lj=(qlj[:,2]==0)&(actual[:,2]==0)
    particle_delta[no_lj,1]=0.
    errors = np.max(particle_delta,axis=0)
    if np.any(errors > [1e-6,1e-6,1e-6]): raise ValueError('Export changed charge or LJ sigma/epsilon: '+str(errors))
    if set(exceptions) != set(actual_exceptions): raise ValueError('Export changed 1-4/exclusion graph')
    exception_error = 0.
    for pair, expected in exceptions.items():
        observed = actual_exceptions[pair]
        # Sigma has no physical effect for a zero-epsilon excluded pair.
        selected = [0,2] if expected[2] == observed[2] == 0 else [0,1,2]
        delta = float(np.max(np.abs(expected[selected]-observed[selected])))
        exception_error = max(exception_error, delta)
    if exception_error > 1e-6: raise ValueError('Export changed 1-4 electrostatic/LJ settings')
    cmap = lambda s: sum(f.getNumTorsions() for f in s.getForces() if type(f).__name__ == 'CMAPTorsionForce')
    if cmap(source_system) != cmap(rebuilt): raise ValueError('Export lost CMAP terms')
    xyz = np.asarray(source.coordinates)
    before, force_before = _single_point(source_system, xyz)
    after, force_after = _single_point(rebuilt, xyz)
    energy_error = abs(after-before); force_error = float(np.max(np.abs(force_after-force_before)))
    energy_tolerance = .005 + abs(before)*1e-7
    force_tolerance = .005 + float(np.max(np.abs(force_before)))*1e-7
    if energy_error > energy_tolerance or force_error > force_tolerance:
        raise ValueError(f'Export changed physical energy/forces: dE={energy_error}, max dF={force_error}')
    return dict(passed=True, atoms=len(source.atoms), bonds=len(source.bonds), cmap_terms=cmap(rebuilt),
        max_charge_error_e=float(errors[0]),max_sigma_error_nm=float(errors[1]),max_epsilon_error_kJ_mol=float(errors[2]),
        zero_LJ_particles=int(no_lj.sum()),max_exception_error=exception_error,vacuum_energy_error_kJ_mol=energy_error,
        max_force_error_kJ_mol_nm=force_error,energy_tolerance_kJ_mol=energy_tolerance,
        force_tolerance_kJ_mol_nm=force_tolerance,coordinates_for_comparison='identical original coordinates')


def export_system(structure, directory, formats=('amber',), force_field='ff14SB'):
    """Export a typed physical structure; reject any unverified parameter loss."""
    formats = tuple(dict.fromkeys(formats))
    if not formats or set(formats)-_FORMATS:
        raise ValueError('Unsupported export format: '+', '.join(str(x) for x in formats))
    profile(force_field)
    import parmed as pmd
    from openmm import XmlSerializer
    from openmm.app import NoCutoff
    system = _physical_system(structure)
    # Capture the physical source independently of all writers.
    typed = pmd.openmm.load_topology(structure.topology,system,xyz=np.asarray(structure.coordinates))
    for a,b in zip(structure.atoms,typed.atoms):
        b.solvent_radius=float(a.solvent_radius);b.screen=float(a.screen)
    source_validation = _validate(typed,system,typed)
    directory = Path(directory);directory.mkdir(parents=True,exist_ok=True)
    record = dict(status='validated_physical_exports',force_field=force_field,
        custom_charge_reference='ff14SB RESP retained across profiles',exports={})
    with tempfile.TemporaryDirectory(prefix='.physical-export-',dir=directory) as temporary:
        dest = Path(temporary)
        (dest/'system.xml').write_text(XmlSerializer.serialize(system))
        typed.save(str(dest/'protein.pdb'),format='pdb',overwrite=False)
        (dest/'state.json').write_text(json.dumps(dict(coordinates_A=np.asarray(typed.coordinates).tolist(),
            atoms=[dict(name=a.name,residue=a.residue.name,residue_index=a.residue.idx,atomic_number=a.atomic_number) for a in typed.atoms],
            bonds=[[b.atom1.idx,b.atom2.idx] for b in typed.bonds],assignment=getattr(structure,'_cationize_assignment',{})),indent=2)+'\n')
        amber_validation = None
        if 'openmm' in formats:
            record['exports']['openmm'] = dict(files=['system.xml','state.json','protein.pdb'],validation=source_validation)
        if 'amber' in formats:
            parm = pmd.amber.AmberParm.from_structure(typed)
            parm.save(str(dest/'system.parm7'),overwrite=False)
            parm.save(str(dest/'system.rst7'),overwrite=False)
            reloaded = pmd.load_file(str(dest/'system.parm7'),xyz=str(dest/'system.rst7'))
            amber_validation = _validate(typed,system,reloaded)
            if 'amber' in formats:
                record['exports']['amber'] = dict(files=['system.parm7','system.rst7'],validation=amber_validation)
        if 'gromacs' in formats:
            from .gromacs_export import write_gromacs
            files=write_gromacs(typed,dest)
            reloaded = pmd.load_file(str(dest/'system.top'),xyz=str(dest/'system.gro'))
            record['exports']['gromacs'] = dict(files=files,validation=_validate(typed,system,reloaded),
                engine_runtime_tested=False)
        from .amber_bundle import write_bundle
        record['complete_parameter_bundle']=write_bundle(typed,system,dest,force_field,
                                     source_types=getattr(structure,'_cationize_atom_types',None))
        record['files_sha256'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(dest.iterdir())}
        (dest/'EXPORTS.json').write_text(json.dumps(record,indent=2)+'\n')
        for p in dest.iterdir():
            if (directory/p.name).exists(): raise FileExistsError('Export refuses to overwrite '+str(directory/p.name))
        for p in dest.iterdir(): os.replace(p,directory/p.name)
    return record
