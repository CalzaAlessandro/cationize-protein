"""Explicit alternative electrostatic models; no silent fallback from PB.

APBS remains the nonlinear PB reference. Amber PBSA is a distinct PB
implementation, while screened Coulomb is a uniform-dielectric approximation
without a protein/solvent interface or reaction field.
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import numpy as np
from .electrostatics import run_apbs, read_dx, APBSError


def screened_coulomb(coordinates, charges, points, *, temperature_K=298., dielectric=78.5, salt_M=.15):
    xyz, q, points = np.asarray(coordinates,float), np.asarray(charges,float), np.asarray(points,float)
    if xyz.ndim != 2 or xyz.shape[1:] != (3,) or points.ndim != 2 or points.shape[1:] != (3,) or q.shape != (len(xyz),):
        raise ValueError('Invalid charge/sample coordinate dimensions')
    if not all(np.isfinite(x).all() for x in (xyz,q,points)) or not np.isfinite([temperature_K,dielectric,salt_M]).all() or min(temperature_K,dielectric)<=0 or salt_M<0:
        raise ValueError('Invalid screened-Coulomb parameters')
    e, eps0, kb, na = 1.602176634e-19, 8.8541878128e-12, 1.380649e-23, 6.02214076e23
    bjerrum_A = e*e/(4*np.pi*eps0*dielectric*kb*temperature_K)*1e10
    kappa_A = np.sqrt(2*e*e*na*1000*salt_M/(eps0*dielectric*kb*temperature_K))*1e-10
    result = np.empty(len(points))
    for offset in range(0,len(points),256):
        r=np.linalg.norm(points[offset:offset+256,None,:]-xyz[None,:,:],axis=2)
        if np.any(r<1e-10):raise ValueError('A sample point coincides with an atomic charge')
        result[offset:offset+256]=bjerrum_A*((np.exp(-kappa_A*r)/r)@q)
    return result


def pbsa_to_thermal(values, temperature_K):
    return np.asarray(values)/(0.00198720425864083*temperature_K)


def pbsa_input(*,temperature_K,salt_M,spacing_A,solute_dielectric,solvent_dielectric,ion_radius_A,fillratio,nonlinear=False):
    return f'''Protein potential only; no dynamics
&cntrl
 imin=1, ntx=1, ipb=2, inp=0,
/
&pb
 npbverb=1, radiopt=0, npbopt={1 if nonlinear else 0}, solvopt={6 if nonlinear else 1}, accept=0.0001, maxitn=10000,
 epsin={solute_dielectric:.12g}, epsout={solvent_dielectric:.12g}, istrng={salt_M*1000:.12g}, pbtemp={temperature_K:.12g},
 dprob=1.4, iprob={ion_radius_A:.12g}, space={spacing_A:.12g}, fillratio={fillratio:.12g}, nfocus=1,
 sasopt=2, phiout=1, phiform=2,
/
'''


def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_pbsa(structure, info, sample_points, out, config, *, keep_map=False,nonlinear=False):
    import parmed as pmd
    from .dependencies import external_program_environment
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    binary=shutil.which(str(config.executable))
    if binary is None:raise APBSError('PBSA executable is unavailable')
    # PBSA consumes charges, coordinates and intrinsic radii from this dedicated
    # input. Round-trip them rather than letting PBSA infer standard radii.
    amber=pmd.amber.AmberParm.from_structure(structure)
    top,rst=out/'protein.parm7',out/'protein.rst7'
    amber.save(str(top),overwrite=False);amber.save(str(rst),format='rst7',overwrite=False)
    restored=pmd.load_file(str(top),xyz=str(rst))
    np.testing.assert_allclose([a.charge for a in restored.atoms],info.charges_e,atol=2e-7,rtol=0)
    np.testing.assert_allclose([a.solvent_radius for a in restored.atoms],info.pb_radii_A,atol=1e-7,rtol=0)
    span=float(np.max(np.ptp(info.coordinates,axis=0)))
    fillratio=max(2.,1+2*config.padding_A/max(span,1.))
    (out/'pbsa.in').write_text(pbsa_input(temperature_K=config.temperature_K,salt_M=config.salt_M,
        spacing_A=config.spacing_A,solute_dielectric=config.solute_dielectric,solvent_dielectric=config.solvent_dielectric,
        ion_radius_A=config.ion_radius_A,fillratio=fillratio,nonlinear=nonlinear))
    started=time.monotonic()
    with (out/'pbsa.log').open('x') as log:
        process=subprocess.run([binary,'-O','-i','pbsa.in','-o','pbsa.out','-p','protein.parm7','-c','protein.rst7'],
            cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=config.timeout_seconds,
            env=external_program_environment(binary))
    text=(out/'pbsa.log').read_text(errors='replace')+((out/'pbsa.out').read_text(errors='replace') if (out/'pbsa.out').exists() else '')
    if process.returncode or re.search(r'(?im)^\s*(?:fatal|error)(?!\s+tolerance)|failed to converg|not converg|maximum.*iterations.*(?:reached|exceed)',text):
        raise APBSError('PBSA failed or did not converge; inspect pbsa.log/pbsa.out')
    path=out/'pbsa_phi.dx'
    if not path.exists():raise APBSError('PBSA did not generate its potential DX map')
    field=read_dx(path,max_items=config.max_grid_points)
    values=pbsa_to_thermal(field.sample(sample_points),config.temperature_K)
    if not np.isfinite(values).all():raise APBSError('Non-finite PBSA potential')
    record=dict(solver='Amber PBSA',physical_model=('nonlinear' if nonlinear else 'linear')+' Poisson-Boltzmann; smooth density surface',
        config={**asdict(config),'executable':str(binary)},potential_units='kBT/e',native_DX_units='kcal/(mol e)',
        conversion='DX/(R_kcal*T)',solver_sha256=_sha(binary),topology_sha256=_sha(top),DX_sha256=_sha(path),
        wall_seconds=time.monotonic()-started,map_retained=keep_map,
        radius_policy='RADIOPT=0 requests topology radii; PHIOUT visualization surface follows PBSA conventions and is not identical to APBS smol.',
        reference='https://ambermd.org/doc12/Amber25.pdf, PBSA sections6.2 and6.5')
    if not keep_map:path.unlink()
    (out/'PB_COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    return values,record


def solve(structure, info, grid, out, config, *, solver='apbs', keep_map=False, provenance=None,pbsa_nonlinear=False):
    out=Path(out)
    if solver=='apbs':
        from .dependencies import external_program_environment
        return run_apbs(info.coordinates,info.charges_e,info.pb_radii_A,grid.points,out,config,
            atom_names=[a.name for a in info.atoms],residue_names=[a.residue_name for a in info.atoms],
            residue_numbers=[a.residue_ordinal for a in info.atoms],provenance=provenance,keep_map=keep_map)
    if solver=='pbsa':return run_pbsa(structure,info,grid.points,out,config,keep_map=keep_map,nonlinear=pbsa_nonlinear)
    if solver!='coulomb':raise ValueError('Unknown electrostatic model: '+str(solver))
    out.mkdir(parents=True,exist_ok=False)
    scores=np.zeros(len(grid.points))
    scores[~grid.blocked]=screened_coulomb(info.coordinates,info.charges_e,grid.points[~grid.blocked],
        temperature_K=config.temperature_K,dielectric=config.solvent_dielectric,salt_M=config.salt_M)
    record=dict(solver='screened Coulomb',physical_model='Uniform-dielectric Debye-Huckel approximation; no dielectric boundary or protein reaction field',
        potential_units='kBT/e',dielectric=config.solvent_dielectric,salt_M=config.salt_M,temperature_K=config.temperature_K,
        accessible_points=int((~grid.blocked).sum()),blocked_points=int(grid.blocked.sum()),
        comparison='Explicit approximate ranking method; not APBS or PBSA and never selected as an automatic fallback')
    np.save(out/'potential_samples.npy',scores)
    (out/'PB_COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    return scores,record
