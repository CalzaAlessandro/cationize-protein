"""Parameter assignment and restrained minimization without Amber executables."""
import io,json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from .catalog import by_name
from .off import read_units

LEGACY_ATOMS={'ASC':{'N':'N1','CA':'C1','C':'C2','O':'O1'},'GLC':{'N':'N1','CA':'C1','C':'C2','O':'O1'},'ASX':{'N':'N1','CA':'C1','C':'C2','O':'O1'},'LET':{'N':'N1','CA':'C1','C':'C2','O':'O1'}}

def template_xml(folder,names=None):
 root=ET.Element('ForceField');residues=ET.SubElement(root,'Residues');count=0
 for p in sorted(Path(folder).glob('*.lib')):
  if names is not None and p.stem not in names:continue
  name=p.stem;spec=by_name(name);unit=read_units(p)[name];r=ET.SubElement(residues,'Residue',name=name)
  for a in unit['atoms']:ET.SubElement(r,'Atom',name=a['name'],type='protein-'+a['amber_type'],charge=str(a['initial_template_charge_e']))
  for a,b in unit['bonds']:ET.SubElement(r,'Bond',atomName1=unit['atoms'][a]['name'],atomName2=unit['atoms'][b]['name'])
  roles=LEGACY_ATOMS.get(name,{})
  if spec['context']!='n_terminal':ET.SubElement(r,'ExternalBond',atomName=roles.get('N','N'))
  if spec['context']!='c_terminal':ET.SubElement(r,'ExternalBond',atomName=roles.get('C','C'))
  if spec['parent']=='CYX':ET.SubElement(r,'ExternalBond',atomName='SG')
  count+=1
 if names is not None and count!=len(set(names)):raise ValueError('One or more requested templates are unavailable')
 return ET.tostring(root,encoding='unicode')

def forcefield(folder,implicit=False,names=None):
 from openmm.app import ForceField
 files=['amber14/protein.ff14SB.xml']
 if implicit:files.append('implicit/gbn2.xml')
 return ForceField(*files,io.StringIO(template_xml(folder,names)))

def parameterize(structure,folder,implicit=False,names=None):
 import parmed
 from openmm.app import NoCutoff
 from openmm import unit
 ff=forcefield(folder,implicit,names);topology=structure.topology
 # GBn2 uses the parent ARG name for guanidinium H radii. Keep explicit
 # modified templates for matching while exposing parent names to GB rules.
 renamed={};residue_templates={}
 if implicit:
  for residue in topology.residues():
   try:spec=by_name(residue.name)
   except StopIteration:continue
   renamed[residue]=residue.name;residue_templates[residue]=residue.name;residue.name=spec['parent']
 try:system=ff.createSystem(topology,nonbondedMethod=NoCutoff,constraints=None,rigidWater=False,removeCMMotion=False,residueTemplates=residue_templates)
 finally:
  for residue,name in renamed.items():residue.name=name
 typed=parmed.openmm.load_topology(topology,system,xyz=np.asarray(structure.coordinates)*unit.angstrom)
 # Bondi radii with the same H(N) adjustment used by Amber mbondi2, in Python.
 from parmed.tools.actions import changeRadii
 changeRadii(typed,'mbondi2').execute()
 return typed,system

def minimize(structure,folder,reference_positions=None,restraint_kcal_A2=10.,tolerance_kcal_A=.1,max_iterations=20000,platform='CPU',names=None):
 from openmm import Context,VerletIntegrator,CustomExternalForce,LocalEnergyMinimizer,Platform,unit
 from .topology import inspect_structure
 typed,system=parameterize(structure,folder,implicit=True,names=names)
 info=inspect_structure(typed)
 ref=np.asarray(reference_positions if reference_positions is not None else typed.coordinates,float)
 if ref.shape!=np.asarray(typed.coordinates).shape:raise ValueError('Reference coordinate mismatch')
 restraint=CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)');restraint.addGlobalParameter('k',restraint_kcal_A2*418.4)
 for key in ('x0','y0','z0'):restraint.addPerParticleParameter(key)
 for i in info.protected_atoms:restraint.addParticle(int(i),(ref[int(i)]/10).tolist())
 system.addForce(restraint);integrator=VerletIntegrator(.001*unit.picoseconds)
 plat=Platform.getPlatformByName(platform);properties={'Threads':'1'} if platform=='CPU' else {}
 context=Context(system,integrator,plat,properties);context.setPositions(np.asarray(typed.coordinates)*unit.angstrom)
 LocalEnergyMinimizer.minimize(context,tolerance_kcal_A*41.84*unit.kilojoules_per_mole/unit.nanometer,max_iterations)
 state=context.getState(getPositions=True,getEnergy=True,getForces=True);xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom));forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole/unit.nanometer))/41.84
 rms=float(np.sqrt(np.mean(np.sum(forces*forces,axis=1))))
 if not np.isfinite(xyz).all() or not np.isfinite(rms):raise ValueError('Non-finite minimized state')
 typed.coordinates=xyz;record=dict(backend='OpenMM GBn2',platform=platform,restraint_kcal_mol_A2=restraint_kcal_A2,RMS_force_kcal_mol_A=rms,target_RMS_force=tolerance_kcal_A,potential_energy_kJ_mol=float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)),forcefield='ff14SB with custom RESP templates',radii='mbondi2 for subsequent APBS; GBn2 uses OpenMM implicit-solvent parameters')
 del context,integrator
 if rms>tolerance_kcal_A*1.05:raise RuntimeError('Minimization did not meet the requested RMS force tolerance: '+str(record))
 return typed,record
