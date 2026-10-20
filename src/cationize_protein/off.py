"""Read Amber OFF units needed for model construction, without copying the library."""
import re,shlex
from pathlib import Path

def read_units(path):
 sections={}
 for line in Path(path).read_text().splitlines():
  m=re.match(r'^!entry\.([^.]+)\.unit\.([^ ]+)',line)
  if m:current=(m[1],m[2]);sections[current]=[]
  elif line.startswith('!'):current=None
  elif 'current' in locals() and current and line.strip():sections[current].append(shlex.split(line))
 result={}
 for name,section in sections:
  if section!='atoms':continue
  atoms=[]
  for row in sections[name,'atoms']:
   z=int(row[6]);z=z if z>0 else {'H':1,'C':6,'N':7,'O':8,'S':16}[row[0].lstrip('0123456789')[0]]
   atoms.append(dict(name=row[0],amber_type=row[1],atomic_number=z,initial_template_charge_e=float(row[7])))
  result[name]=dict(atoms=atoms,positions=[[float(v) for v in row] for row in sections.get((name,'positions'),[])],bonds=[(int(row[0])-1,int(row[1])-1) for row in sections.get((name,'connectivity'),[])])
 return result
