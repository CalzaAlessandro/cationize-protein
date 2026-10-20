"""Selected-engine restrained minimization, always serial and never MD."""
import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np


def amber_executable(engine='amber', executable=None):
    choices=[executable] if executable else ['sander','pmemd'] if engine=='amber' else [engine]
    for name in choices:
        found=shutil.which(str(name))
        if not found:
            found=str(Path(os.environ.get('AMBERHOME',''))/'bin'/str(name))
            if not Path(found).is_file() or not os.access(found,os.X_OK):continue
        basename=Path(found).name.lower()
        if 'mpi' in basename or 'cuda' in basename or basename not in {'sander','pmemd'}:
            raise ValueError('The portable Amber minimizer requires serial sander or pmemd')
        return str(Path(found).resolve())
    raise RuntimeError('Serial Amber sander/pmemd was not found; select --engine openmm or provide --engine-executable')


def group_restraints(indices, weight):
    indices=sorted(set(int(i) for i in indices))
    if not indices or indices[0]<0 or not math.isfinite(weight) or weight<=0:raise ValueError('Invalid restraint group')
    ranges=[]
    for index in indices:
        if ranges and index+1==ranges[-1][1]+1:ranges[-1][1]=index+1
        else:ranges.append([index+1,index+1])
    lines=['Native backbone and disulfide atoms',f'{weight:.12g}']
    for offset in range(0,len(ranges),7):
        lines.append('ATOM '+' '.join(str(i) for pair in ranges[offset:offset+7] for i in pair))
    return '\n'.join(lines+['END','END'])+'\n'


def amber_input(indices, *,restraint,tolerance,max_iterations):
    # Amber DRMS is the component RMS; OpenMM tolerance is the RMS vector norm.
    drms=tolerance/math.sqrt(3)
    return (f'Protein restrained minimization only\n&cntrl\n imin=1,ntx=1,irest=0,ntb=0,igb=8,saltcon=0.0,\n'
            f' cut=9999.,maxcyc={max_iterations},ncyc={min(500,max_iterations//2)},ntmin=1,drms={drms:.12g},\n'
            f" ntr=1,restraintmask='',ntpr=100,ntxo=1,\n/\n"+group_restraints(indices,restraint))


def parse_amber_minimum(text,tolerance):
    complete='Total wall time:' in text or ('TIMINGS' in text and re.search(r'Run\s+done\s+at',text))
    if 'FINAL RESULTS' not in text or not complete:
        raise RuntimeError('Amber minimization did not finish; inspect min.out')
    tail=text.rsplit('FINAL RESULTS',1)[1]
    match=re.search(r'(?m)^\s*(\d+)\s+([-+\w.]+)\s+([-+\w.]+)\s+([-+\w.]+)\s+\S+\s+\d+',tail)
    if not match:raise RuntimeError('Missing final Amber energy/gradient')
    values=[float(v.replace('D','E')) for v in match.groups()[1:]]
    if not np.isfinite(values).all() or values[1]*math.sqrt(3)>tolerance*1.05:
        raise RuntimeError('Amber minimization has a nonfinite or unconverged final gradient')
    return dict(steps=int(match.group(1)),energy_kcal_mol=values[0],
                component_RMS_force_kcal_mol_A=values[1],RMS_force_kcal_mol_A=values[1]*math.sqrt(3),
                target_RMS_force=tolerance)


def minimize(structure,folder,*,force_field='ff14SB',engine='openmm',executable=None,
             directory=None,reference_positions=None,restraint_kcal_A2=10.,
             tolerance_kcal_A=.1,max_iterations=20000,platform='CPU',names=None):
    from .parameter_assignment import parameterize
    from .topology import inspect_structure
    if engine=='openmm':
        from openmm import Context,VerletIntegrator,CustomExternalForce,LocalEnergyMinimizer,Platform,unit
        typed,system=parameterize(structure,folder,force_field=force_field,implicit=True,names=names)
        info=inspect_structure(typed)
        ref=np.asarray(reference_positions if reference_positions is not None else typed.coordinates,float)
        if ref.shape!=info.coordinates.shape or not np.isfinite(ref).all():raise ValueError('Invalid reference coordinates')
        restraint=CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint.addGlobalParameter('k',restraint_kcal_A2*418.4)
        for name in ('x0','y0','z0'):restraint.addPerParticleParameter(name)
        for i in info.protected_atoms:restraint.addParticle(int(i),(ref[int(i)]/10).tolist())
        system.addForce(restraint)
        integrator=VerletIntegrator(.001*unit.picoseconds)
        properties={'Threads':'1'} if platform=='CPU' else {}
        context=Context(system,integrator,Platform.getPlatformByName(platform),properties)
        try:
            context.setPositions(info.coordinates*unit.angstrom)
            LocalEnergyMinimizer.minimize(context,tolerance_kcal_A*41.84*unit.kilojoules_per_mole/unit.nanometer,max_iterations)
            state=context.getState(getPositions=True,getEnergy=True,getForces=True)
            xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
            forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole/unit.nanometer))/41.84
            rms=float(np.sqrt(np.mean(np.sum(forces*forces,axis=1))))
            energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
            if not np.isfinite(xyz).all() or not np.isfinite([rms,energy]).all() or rms>tolerance_kcal_A*1.05:
                raise RuntimeError(f'OpenMM minimization did not converge: RMS force {rms:.6g}')
            typed.coordinates=xyz
        finally:del context,integrator
        # Never attach GB/restraint forces to a structure used by export code.
        typed,_=parameterize(typed,folder,force_field=force_field,names=names)
        return typed,dict(backend='OpenMM GBn2',force_field=force_field,implicit_salt_M=0.,platform=platform,threads=1 if platform=='CPU' else None,
            restraint_kcal_mol_A2=restraint_kcal_A2,RMS_force_kcal_mol_A=rms,
            target_RMS_force=tolerance_kcal_A,potential_energy_kJ_mol=energy)
    import parmed
    from parmed.tools.actions import changeRadii
    from .catalog import by_name
    binary=amber_executable(engine,executable)
    if directory is None:raise ValueError('Amber minimization needs its own output directory')
    out=Path(directory);out.mkdir(parents=True,exist_ok=False)
    typed,_=parameterize(structure,folder,force_field=force_field,names=names)
    info=inspect_structure(typed)
    ref=np.asarray(reference_positions if reference_positions is not None else info.coordinates,float)
    if ref.shape!=info.coordinates.shape or not np.isfinite(ref).all():raise ValueError('Invalid reference coordinates')
    parm=parmed.amber.AmberParm.from_structure(typed)
    original_names=[r.name for r in parm.residues]
    # Amber GBn2 uses mbondi3; APBS continues to use the separately assigned mbondi2 radii.
    for residue in parm.residues:
        try:residue.name=by_name(residue.name)['parent']
        except StopIteration:pass
    changeRadii(parm,'mbondi3').execute()
    for residue,name in zip(parm.residues,original_names):residue.name=name
    parm.save(str(out/'protein.parm7'),overwrite=False)
    parm.save(str(out/'input.rst7'),format='rst7',overwrite=False)
    parm.coordinates=ref;parm.save(str(out/'reference.rst7'),format='rst7',overwrite=False)
    (out/'min.in').write_text(amber_input(info.protected_atoms,restraint=restraint_kcal_A2,
        tolerance=tolerance_kcal_A,max_iterations=max_iterations))
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (out/'engine.log').open('x') as log:
        result=subprocess.run([binary,'-O','-i','min.in','-o','min.out','-p','protein.parm7',
            '-c','input.rst7','-ref','reference.rst7','-r','min.rst7','-inf','min.mdinfo'],
            cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
    if result.returncode:raise RuntimeError('Amber minimization failed; inspect '+str(out/'engine.log'))
    record=parse_amber_minimum((out/'min.out').read_text(),tolerance_kcal_A)
    xyz=np.asarray(parmed.load_file(str(out/'protein.parm7'),xyz=str(out/'min.rst7')).coordinates)
    if xyz.shape!=info.coordinates.shape or not np.isfinite(xyz).all():raise RuntimeError('Invalid Amber minimum coordinates')
    typed.coordinates=xyz
    typed,_=parameterize(typed,folder,force_field=force_field,names=names)
    record.update(backend='Amber GBn2',force_field=force_field,implicit_salt_M=0.,executable=binary,
        executable_sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest(),threads=1,
        restraint_atoms=list(info.protected_atoms),restraint_kcal_mol_A2=restraint_kcal_A2,
        radii='Amber GBn2 uses mbondi3; PB ranking uses separately assigned mbondi2')
    return typed,record
