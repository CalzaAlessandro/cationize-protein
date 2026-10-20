"""One template: two validated Gaussian CPU calculations followed by RESP."""
import argparse,json
from pathlib import Path
from .qm import run
from .resp import fit

def worker(project,name):
 project=Path(project).resolve()
 for label in ('alpha','extended'):
  print(f'[QM] {name} / {label}',flush=True);run(project/'work/models'/name/label,project/'work/qm'/name/label)
 print(f'[RESP] {name}: geometry, charge, ESP and capped topology checks',flush=True)
 result=fit(project,name);print(json.dumps(result,indent=2),flush=True)
 from .parameter_validation import validate
 print('[EXPORT] Checking custom template in OpenMM ff14SB and GBn2',flush=True)
 validate(project,name)
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('project',type=Path);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--name');g.add_argument('--index',type=int);a=p.parse_args()
 name=a.name if a.name else json.loads((a.project/'work/array_templates.json').read_text())[a.index]
 worker(a.project,name)
