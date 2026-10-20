"""Two-stage multi-conformer RESP with fixed canonical caps and topology export."""
import json,re,subprocess,hashlib,shutil
from pathlib import Path
import numpy as np
from .qm import completed,digest
BOHR=0.529177210903

def command(args,folder,log):
 with (folder/log).open('x') as f:subprocess.run(list(map(str,args)),cwd=folder,stdout=f,stderr=subprocess.STDOUT,check=True)

def mk_coordinates(text,natoms):
 rows=[]
 for m in re.finditer(r'Atomic Center\s+(\d+) is at\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)',text):
  if int(m[1])==1:rows=[]
  if int(m[1])!=len(rows)+1:raise ValueError('Invalid MK atom order')
  rows.append([float(m[i]) for i in (2,3,4)])
 if len(rows)!=natoms:raise ValueError('Missing final MK coordinates')
 return np.array(rows)

def read_esp(path):
 lines=Path(path).read_text().splitlines();header=list(map(int,lines[0].split()))
 if len(header) not in (2,3) or (len(header)==3 and header[2]!=0):raise ValueError('Unsupported ESP header or multipole data')
 nat,npts=header[:2]
 xyz=np.array([[float(x) for x in l.split()] for l in lines[1:1+nat]])
 grid=np.array([[float(x) for x in l.split()] for l in lines[1+nat:1+nat+npts]])
 if nat<1 or npts<1 or xyz.shape!=(nat,3) or grid.shape!=(npts,4) or not np.isfinite(xyz).all() or not np.isfinite(grid).all():raise ValueError('Malformed ESP')
 return xyz,grid[:,1:],grid[:,0]

def safe_ac_names(atoms):
 symbols={1:'H',6:'C',7:'N',8:'O',16:'S'}
 if len(atoms)>999:raise ValueError('Model too large for four-character surrogate names')
 return [symbols[a['atomic_number']]+f'{i:03d}' for i,a in enumerate(atoms,1)]

def validate_resp_elements(text,expected,nmol):
 lines=[x.strip() for x in text.splitlines() if x.strip()]
 try:at=next(i for i,x in enumerate(lines) if x.lower()=='&end')+1
 except StopIteration:raise ValueError('RESP control namelist missing')
 blocks=[]
 while at+2<len(lines):
  # Each molecule starts with weight, a textual title, then charge and atom count.
  try:weight=float(lines[at]);header=[int(v) for v in lines[at+2].split()]
  except ValueError:at+=1;continue
  if weight<=0 or not any(c.isalpha() for c in lines[at+1]) or len(header)!=2 or header[1]!=len(expected):at+=1;continue
  rows=lines[at+3:at+3+len(expected)]
  try:z=[int(row.split()[0]) for row in rows]
  except (ValueError,IndexError):raise ValueError('Malformed RESP atomic-number block')
  if z!=list(expected):raise ValueError(f'RESP elements differ: {z} != {list(expected)}')
  blocks.append(z);at+=3+len(expected)
 if len(blocks)!=nmol:raise ValueError(f'Expected {nmol} RESP molecule blocks, found {len(blocks)}')
 return True

def fit(project,name):
 from rdkit import Chem
 project=Path(project).resolve();models=project/'work/models'/name;m=json.loads((models/'atom_mapping.json').read_text());atoms=m['atoms'];nat=len(atoms)
 root=project/'work/fits'/name;root.mkdir(parents=True,exist_ok=True)
 marker=root/'FIT_COMPLETE.json';library=project/'parameters'/(name+'.lib')
 if marker.exists():
  r=json.loads(marker.read_text())
  if digest(library)!=r['library_sha256'] or r.get('model_sha256')!=digest(models/'atom_mapping.json') or r.get('fit_code_sha256')!=digest(Path(__file__)):
   raise RuntimeError('Existing fit provenance differs; preserve it and use an isolated refit project')
  for label,expected in r.get('qm_input_sha256',{}).items():
   if digest(models/label/'opt_esp.gjf')!=expected:raise RuntimeError('QM input changed after fitting')
  return r
 # Preserve failed partial fit attempts; a new attempt never overwrites them.
 attempt=max([int(p.name.split('_')[-1]) for p in root.glob('attempt_*') if p.is_dir()]+[0])+1
 out=root/f'attempt_{attempt:03d}';out.mkdir();esps=[];positions=[];evidence=[]
 for i,label in enumerate(('alpha','extended'),1):
  qm=completed(models/label,project/'work/qm'/name/label)
  if qm is None:raise RuntimeError('Missing validated QM '+name+'/'+label)
  xyz=mk_coordinates((qm/'gaussian.log').read_text(),nat)
  template=Chem.MolFromMolFile(str(models/label/'initial.mol'),removeHs=False,sanitize=True);conf=template.GetConformer()
  for j,x in enumerate(xyz):conf.SetAtomPosition(j,x)
  Chem.AssignStereochemistryFrom3D(template)
  for j,cip in m['chirality_CIP'].items():
   atom=template.GetAtomWithIdx(int(j)-1)
   if not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode')!=cip:raise ValueError('QM stereochemistry changed')
  for a,b in m['bonds']:
   z=[atoms[a-1]['atomic_number'],atoms[b-1]['atomic_number']];d=np.linalg.norm(xyz[a-1]-xyz[b-1]);limit=1.35 if 1 in z else 2.35 if 16 in z else 1.9
   if not .7<d<limit:raise ValueError(f'QM bond/proton transfer {a}-{b}: {d}')
  esp=out/f'conf_{i}.esp';command(['espgen','-i',qm/'gaussian.log','-o',esp],out,f'espgen_{i}.log');data=read_esp(esp)
  np.testing.assert_allclose(data[0]*BOHR,xyz,atol=2e-5,rtol=0)
  positions.append(xyz);esps.append(data);evidence.append(dict(conformation=label,qm_log_sha256=digest(qm/'gaussian.log'),points=len(data[1])))
 acnames=safe_ac_names(atoms)
 lines=[f'CHARGE {m["formal_charge"]:9.2f} ( {m["formal_charge"]} )\n']
 for atom,x in zip(atoms,positions[0]):
  lines.append(f"ATOM  {atom['index']:5d} {acnames[atom['index']-1]:<4s} {atom['residue']:3s} {atom['residue_index']:5d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{atom['initial_template_charge_e']:10.6f}{atom['amber_type']:>10s}\n")
 for serial,((a,b),order) in enumerate(zip(m['bonds'],m['bond_orders']),1):
  lines.append(f"BOND{serial:5d}{a:5d}{b:5d}{(9 if order==1.5 else int(order)):5d}     {acnames[a-1]} {acnames[b-1]}\n")
 (out/'model.ac').write_text(''.join(lines));(out/'all.esp').write_text(''.join((out/f'conf_{i}.esp').read_text() for i in (1,2)))
 fixed={v for group in m['caps'] for v in group}
 (out/'fixed_caps.dat').write_text(''.join(f"CHARGE {a['initial_template_charge_e']:.8f} {a['index']} {acnames[a['index']-1]}\n" for a in atoms if a['index'] in fixed))
 for stage in (1,2):
  command(['respgen','-i','model.ac','-o',f'stage{stage}.in','-f',f'resp{stage}','-n','2','-a','fixed_caps.dat'],out,f'respgen{stage}.log')
  validate_resp_elements((out/f'stage{stage}.in').read_text(),[a['atomic_number'] for a in atoms],2)
  command(['resp','-i',f'stage{stage}.in','-o',f'stage{stage}.out','-p',f'stage{stage}.punch','-e','all.esp','-q','QIN' if stage==1 else 'stage1.charges','-t',f'stage{stage}.charges','-s',f'stage{stage}.esout'],out,f'resp{stage}.log')
  if re.search(r'NaN|Singular|not converg|Error',(out/f'stage{stage}.out').read_text(),re.I):raise ValueError('RESP fitting did not converge')
 q=np.array([float(x) for x in (out/'stage2.charges').read_text().split()]).reshape(2,nat)
 if not np.isfinite(q).all() or not np.allclose(q[0],q[1],atol=2e-6):raise ValueError('Conformer charge equivalence failed')
 q=q[0];central=np.array(m['central_atoms'])-1
 if abs(q[central].sum()-m['formal_charge'])>3e-5:raise ValueError('Residue net charge mismatch')
 if any(abs(q[i-1]-atoms[i-1]['initial_template_charge_e'])>1e-6 for i in fixed):raise ValueError('Cap charges changed')
 metrics=[]
 for label,(coords,points,pot) in zip(('alpha','extended'),esps):
  pred=(1/np.linalg.norm(points[:,None,:]-coords[None,:,:],axis=2))@q;rrms=float(np.sqrt(np.sum((pred-pot)**2)/np.sum(pot**2)))
  if not np.isfinite(rrms) or rrms>.3:raise ValueError(f'ESP relative RMS too high: {rrms}; inspect chemistry and fitting')
  metrics.append(dict(conformation=label,relative_RMS=rrms))
 ids={int(old)+1:i+1 for i,old in enumerate(central)};bonds=[(ids[a],ids[b],order) for (a,b),order in zip(m['bonds'],m['bond_orders']) if a in ids and b in ids]
 mol2=f'@<TRIPOS>MOLECULE\n{name}\n{len(central)} {len(bonds)} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n'
 for i,j in enumerate(central,1):
  a=atoms[j];x=positions[0][j];mol2+=f'{i} {a["name"]} {x[0]:.8f} {x[1]:.8f} {x[2]:.8f} {a["amber_type"]} 1 {name} {q[j]:.8f}\n'
 mol2+='@<TRIPOS>BOND\n'+''.join(f'{i} {a} {b} {"ar" if o==1.5 else int(o)}\n' for i,(a,b,o) in enumerate(bonds,1))+f'@<TRIPOS>SUBSTRUCTURE\n1 {name} 1 RESIDUE 0 A {name} 0 ROOT\n';(out/(name+'.mol2')).write_text(mol2)
 # Build the complete fitted capped model with explicit atoms/bonds, including disulfide partner.
 allmol=f'@<TRIPOS>MOLECULE\nmodel\n{nat} {len(m["bonds"])} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n'
 for i,(a,x) in enumerate(zip(atoms,positions[0]),1):allmol+=f'{i} A{i} {x[0]:.8f} {x[1]:.8f} {x[2]:.8f} {a["amber_type"]} 1 MOD {q[i-1]:.8f}\n'
 allmol+='@<TRIPOS>BOND\n'+''.join(f'{i} {a} {b} {"ar" if o==1.5 else int(o)}\n' for i,((a,b),o) in enumerate(zip(m['bonds'],m['bond_orders']),1))+'@<TRIPOS>SUBSTRUCTURE\n1 MOD 1 RESIDUE 0 A MOD 0 ROOT\n';(out/'complete.mol2').write_text(allmol)
 text=f'source leaprc.protein.ff14SB\n{name} = loadmol2 {out/name}.mol2\n'
 head='null' if m['context']=='n_terminal' else f'{name}.1.N';tail='null' if m['context']=='c_terminal' else f'{name}.1.C'
 text+=f'set {name} head {head}\nset {name} tail {tail}\nset {name}.1 connect0 {head}\nset {name}.1 connect1 {tail}\n'
 if m['parent']=='CYX':text+=f'set {name}.1 connect2 {name}.1.SG\n'
 text+=f'saveoff {name} {out/name}.lib\np = loadmol2 {out/"complete.mol2"}\ncheck p\nsaveamberparm p {out/"capped.parm7"} {out/"capped.rst7"}\nquit\n';(out/'tleap.in').write_text(text)
 command(['tleap','-f','tleap.in'],out,'tleap.log')
 log=(out/'tleap.log').read_text()
 if not (out/'capped.parm7').exists() or re.search(r'Could not find|does not have a type|FATAL|ATTN',log,re.I):raise ValueError('Missing/unvalidated Amber bonded terms: '+str(out))
 # Publish only the validated template, atom map, and provenance, not source Amber libraries.
 if library.exists():raise FileExistsError(library)
 shutil.copy2(out/(name+'.lib'),library)
 mapping={k:v for k,v in m.items() if k not in ('atoms','bonds','bond_orders','central_atoms','caps','chirality_CIP','torsions')}
 mapping['atoms']=[dict(atoms[j],index=i+1,model_index=int(j)+1,residue_index=1,charge_e=float(q[j])) for i,j in enumerate(central)]
 mapping['bonds']=[(a,b,o) for a,b,o in bonds];mapping['central_atoms']=list(range(1,len(central)+1));mapping['caps']=[]
 mapping['chirality_CIP']={str(ids[int(k)]):v for k,v in m['chirality_CIP'].items() if int(k) in ids}
 mapping['torsions']=[[ids[k] for k in t] for t in m['torsions'] if all(k in ids for k in t)]
 mapping['model_provenance']=dict(atom_mapping_sha256=digest(models/'atom_mapping.json'),residue_to_capped_model_atom_ids=[int(j)+1 for j in central],cap_atom_ids=m['caps'],full_model_torsions=m['torsions'])
 (project/'parameters'/(name+'.json')).write_text(json.dumps(mapping,indent=2)+'\n')
 result=dict(name=name,status='RESP_and_capped_topology_validated',library_sha256=digest(library),charge=float(q[central].sum()),metrics=metrics,QM=evidence,forcefield='ff14SB atom types and bonded analogy; torsions not refitted',geometry_validation='QM bonds and inherited stereochemistry checked; no full-protein transferability claim',model_sha256=digest(models/'atom_mapping.json'),fit_code_sha256=digest(Path(__file__)),qm_input_sha256={label:digest(models/label/'opt_esp.gjf') for label in ('alpha','extended')})
 marker.write_text(json.dumps(result,indent=2)+'\n');return result
