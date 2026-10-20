"""Build exact capped chemical models from the user's installed Amber templates."""
import json,os,hashlib
from pathlib import Path
import numpy as np
from .catalog import templates,by_name
from .off import read_units

def amidation_atom_types(atoms,carbon,oxygen):
 # The reacted COO- carbon and retained oxygen now form a neutral amide,
 # with the ff14SB C/O types used by the corresponding ASN/GLN carbonyl.
 if atoms[carbon]['atomic_number']!=6 or atoms[oxygen]['atomic_number']!=8:
  raise ValueError('Amidation requires a carbonyl carbon and oxygen')
 atoms[carbon]['amber_type']='C';atoms[oxygen]['amber_type']='O'

def build(spec,units):
 from rdkit import Chem
 from rdkit.Chem import rdDetermineBonds,AllChem
 atoms=[];edges=[];positions=[]
 def fragment(unit,residue,label,central=False):
  source=units[unit];offset=len(atoms);lut={}
  for i,row in enumerate(source['atoms']):
   lut[row['name']]=offset+i
   atoms.append(dict(row,residue=residue,residue_index=label,central=central));positions.append(source['positions'][i])
  edges.extend((a+offset,b+offset) for a,b in source['bonds'])
  return lut
 name=spec['name'];context=spec['context'];parent=spec['parent']
 # A CYX model includes an explicitly disulfide-linked native capped Cys partner.
 unit=spec['source_unit'];special_thiolate=unit not in units and parent=='CYM'
 if special_thiolate:unit='CCYS'
 if unit not in units:raise ValueError('Missing Amber unit '+unit)
 central=fragment(unit,name,2,True);fixed=[]
 if context!='n_terminal':
  cap=fragment('ACE','ACE',1);edges.append((cap['C'],central['N']));fixed.extend(cap.values())
 if context!='c_terminal':
  cap=fragment('NME','NME',3);edges.append((central['C'],cap['N']));fixed.extend(cap.values())
 if parent=='CYX':
  partner=fragment('CYX','CYX',5);a=fragment('ACE','ACE',4);b=fragment('NME','NME',6)
  edges.extend([(a['C'],partner['N']),(partner['C'],b['N']),(central['SG'],partner['SG'])]);fixed.extend(partner.values());fixed.extend(a.values());fixed.extend(b.values())
 remove=set()
 if special_thiolate:
  remove.add(central['HG'])
  native_s=next(a for a in units['CYM']['atoms'] if a['name']=='SG')
  atoms[central['SG']]['amber_type']=native_s['amber_type']
  atoms[central['SG']]['initial_template_charge_e']=native_s['initial_template_charge_e']
 def modify(carbon,oxygen,prefix):
  remove.add(central[oxygen]);last=central[carbon]
  retained={'OD2':'OD1','OE2':'OE1','OXT':'O'}[oxygen]
  amidation_atom_types(atoms,central[carbon],central[retained])
  # Preserve all unchanged atom types; only the explicitly introduced group is new.
  pattern=[('N','N',1),('C','CT',2),('C','CT',2),('N','N3',3)]
  for j,(element,atype,nh) in enumerate(pattern,1):
   idx=len(atoms);an=f'{prefix}{j}'
   atoms.append(dict(name=an,amber_type=atype,atomic_number={'N':7,'C':6}[element],initial_template_charge_e=0.,residue=name,residue_index=2,central=True));positions.append([0.,0.,0.]);edges.append((last,idx));last=idx
   for h in range(nh):
    hi=len(atoms);ht='H' if element=='N' else ('H1' if j==2 else 'HP')
    atoms.append(dict(name=f'H{prefix}{j}{h+1}',amber_type=ht,atomic_number=1,initial_template_charge_e=0.,residue=name,residue_index=2,central=True));positions.append([0.,0.,0.]);edges.append((idx,hi))
 for group in spec['modified_groups']:
  if group=='side':modify('CG' if parent=='ASP' else 'CD','OD2' if parent=='ASP' else 'OE2','S')
  else:modify('C','OXT','T')
 keep=[i for i in range(len(atoms)) if i not in remove];renumber={old:new for new,old in enumerate(keep)}
 atoms=[atoms[i] for i in keep];edges=[(renumber[a],renumber[b]) for a,b in edges if a in renumber and b in renumber];positions=np.array([positions[i] for i in keep]);central={k:renumber[v] for k,v in central.items() if v in renumber};fixed=[renumber[i] for i in fixed]
 rw=Chem.RWMol()
 for a in atoms:rw.AddAtom(Chem.Atom(a['atomic_number']))
 for a,b in edges:rw.AddBond(a,b,Chem.BondType.SINGLE)
 mol=rw.GetMol();conf=Chem.Conformer(len(atoms))
 for i,x in enumerate(positions):conf.SetAtomPosition(i,x)
 mol.AddConformer(conf)
 rdDetermineBonds.DetermineBondOrders(mol,charge=spec['formal_charge'],allowChargedFragments=True,embedChiral=False)
 mol.UpdatePropertyCache(strict=False)
 Chem.SanitizeMol(mol)
 # Preserve all inherited alpha/beta stereochemistry from Amber's exact geometry.
 # Only centers whose neighbors already have known coordinates are assigned here.
 Chem.AssignAtomChiralTagsFromStructure(mol)
 for a in mol.GetAtoms():
  if not atoms[a.GetIdx()]['name'] in ('CA','CB'):a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
 Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
 expected={i:a.GetProp('_CIPCode') for i,a in enumerate(mol.GetAtoms()) if a.HasProp('_CIPCode')}
 Chem.SanitizeMol(mol)
 if Chem.GetFormalCharge(mol)!=spec['formal_charge']:raise ValueError('Incorrect model formal charge')
 if abs(sum(atoms[i]['initial_template_charge_e'] for i in fixed))>2e-5:raise ValueError('Caps/partner are not charge neutral')
 for i,a in enumerate(atoms):a['index']=i+1
 phi=[]
 if parent=='PRO':
  tip=next(i for i,a in enumerate(atoms) if a['central'] and a['name']=='T1');phi.append([central['N'],central['CA'],central['C'],tip])
 elif context!='n_terminal':
  ace=next(i for i,a in enumerate(atoms) if a['residue_index']==1 and a['name']=='C');phi.append([ace,central['N'],central['CA'],central['C']])
 else:
  nme=next(i for i,a in enumerate(atoms) if a['residue']=='NME' and a['name']=='N');phi.append([central['N'],central['CA'],central['C'],nme])
 mapping=dict(**spec,atoms=atoms,bonds=[[a+1,b+1] for a,b in edges],bond_orders=[mol.GetBondBetweenAtoms(a,b).GetBondTypeAsDouble() for a,b in edges],caps=[[i+1 for i in fixed]],central_atoms=[i+1 for i,a in enumerate(atoms) if a['central']],chirality_CIP={str(i+1):x for i,x in expected.items()},torsions=[[v+1 for v in t] for t in phi],charge=spec['formal_charge'],multiplicity=1,source_note='Amber installed residue types and native atoms; ACE/NME or explicit neutral disulfide partner fixed during RESP')
 return mol,mapping

def prepare(output,amberhome=None,names=None):
 from rdkit import Chem
 from rdkit.Chem import AllChem,rdMolTransforms
 amberhome=Path(amberhome or os.environ['AMBERHOME']);lib=amberhome/'dat/leap/lib';units={}
 sources={}
 for filename in ('amino12.lib','aminont12.lib','aminoct12.lib'):
  path=lib/filename;units.update(read_units(path));sources[filename]=hashlib.sha256(path.read_bytes()).hexdigest()
 output=Path(output);output.mkdir(parents=True,exist_ok=True);records=[]
 generation_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 for spec in templates():
  name=spec['name']
  if spec['reused_from'] or (names and name not in names):continue
  dest=output/name
  if dest.exists():
   if (dest/'atom_mapping.json').exists() and json.loads((dest/'atom_mapping.json').read_text()).get('generation_sha256')==generation_sha and all((dest/c/'opt_esp.gjf').exists() for c in ('alpha','extended')):
    old=json.loads((dest/'atom_mapping.json').read_text());records.append(dict(name=name,atoms=len(old['atoms']),charge=old['formal_charge']));continue
   failed=output.parent/'model_failures';failed.mkdir(exist_ok=True);index=1
   while (failed/f'{name}_{index:03d}').exists():index+=1
   dest.rename(failed/f'{name}_{index:03d}')
  mol,mapping=build(spec,units);mapping['source_library_sha256']=sources;mapping['generation_sha256']=generation_sha;dest.mkdir()
  for j,label in enumerate(('alpha','extended')):
   m=Chem.Mol(mol);m.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=27411+j;params.enforceChirality=True
   if AllChem.EmbedMolecule(m,params)!=0:raise ValueError('Embedding failed '+name)
   if AllChem.MMFFHasAllMoleculeParams(m):AllChem.MMFFOptimizeMolecule(m,mmffVariant='MMFF94s',maxIters=2000);method='MMFF94s'
   else:AllChem.UFFOptimizeMolecule(m,maxIters=2000);method='UFF'
   target=(-45. if spec['context']=='n_terminal' or spec['parent']=='PRO' else -60.) if j==0 else 135.
   for t in mapping['torsions']:rdMolTransforms.SetDihedralDeg(m.GetConformer(),*(i-1 for i in t),target)
   Chem.AssignStereochemistryFrom3D(m)
   for i,expected in mapping['chirality_CIP'].items():
    atom=m.GetAtomWithIdx(int(i)-1)
    if not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode')!=expected:raise ValueError('Stereochemistry changed '+name)
   xyz=np.array(m.GetConformer().GetPositions());folder=dest/label;folder.mkdir();Chem.MolToMolFile(m,str(folder/'initial.mol'))
   for i,(a,b) in enumerate(mapping['bonds']):
    d=np.linalg.norm(xyz[a-1]-xyz[b-1]);limit=2.3 if 16 in (mapping['atoms'][a-1]['atomic_number'],mapping['atoms'][b-1]['atomic_number']) else 1.9
    if not .7<d<limit:raise ValueError(f'Bad starting bond {name} {a} {b}: {d}')
   route='#p HF/6-31G(d) Opt=(ModRedundant,CalcFC,MaxCycles=200) SCF=(Tight,XQC) NoSymm Pop=MK IOp(6/33=2,6/42=6)'
   text=f'%chk=state.chk\n%mem=12GB\n%nprocshared=16\n{route}\n\n{name} {label} capped model\n\n{spec["formal_charge"]} 1\n'
   text+=''.join(f'{a["atomic_number"]} {x[0]:.9f} {x[1]:.9f} {x[2]:.9f}\n' for a,x in zip(mapping['atoms'],xyz))
   constraints=''.join('D '+' '.join(map(str,t))+' F\n' for t in mapping['torsions'])
   # Preserve the specified protonation state during gas-phase geometry optimization.
   for atom in m.GetAtoms():
    if atom.GetAtomicNum() in (7,8,16):
     for neighbor in atom.GetNeighbors():
      if neighbor.GetAtomicNum()==1:constraints+=f'B {atom.GetIdx()+1} {neighbor.GetIdx()+1} F\n'
   text+='\n'+constraints+'\n';(folder/'opt_esp.gjf').write_text(text)
   (folder/'initial_geometry.json').write_text(json.dumps(dict(method=method,formal_charge=spec['formal_charge'],frozen_angle_deg=target,production_forcefield=False),indent=2))
  (dest/'atom_mapping.json').write_text(json.dumps(mapping,indent=2)+'\n');records.append(dict(name=name,atoms=len(mapping['atoms']),charge=spec['formal_charge']))
 (output/'models.json').write_text(json.dumps(records,indent=2)+'\n');return records
