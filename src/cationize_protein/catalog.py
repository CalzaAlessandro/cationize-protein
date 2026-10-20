"""Complete, explicit chemical states for fixed-ammonium carboxyl amidation."""
TERMINALS={'ALA':'ALT','ARG':'ART','ASN':'ANT','ASP':'AST','CYS':'CYT','GLN':'GNT','GLU':'GLT','GLY':'GYT','HIE':'HET','ILE':'ILT','LEU':'LET','LYS':'LYT','MET':'MTT','PHE':'PHT','PRO':'PRT','SER':'SET','THR':'THT','TRP':'TRT','TYR':'TYT','VAL':'VAT','HID':'HDT','HIP':'HPT','CYX':'CXT','CYM':'CMT'}
SIDE_CHARGE={'ARG':1,'LYS':1,'HIP':1,'ASP':-1,'GLU':-1,'CYM':-1}
REUSED={'ASC':'ASL','GLC':'GLS','ASX':'NSL','LET':'LCL'}

def templates():
 rows=[]
 def add(name,parent,context,groups):
  native=SIDE_CHARGE.get(parent,0)+(1 if context=='n_terminal' else -1 if context=='c_terminal' else 0)
  rows.append(dict(name=name,parent=parent,context=context,modified_groups=sorted(groups),formal_charge=native+2*len(groups),source_unit=('N' if context=='n_terminal' else 'C' if context=='c_terminal' else '')+parent,reused_from=REUSED.get(name),disulfide=parent=='CYX'))
 for parent,name,nterm in [('ASP','ASC','ASX'),('GLU','GLC','GLX')]:
  add(name,parent,'internal',{'side'});add(nterm,parent,'n_terminal',{'side'})
 for parent,name in TERMINALS.items():add(name,parent,'c_terminal',{'terminal'})
 for parent,side,both in [('ASP','ASZ','ASD'),('GLU','GLZ','GLD')]:
  add(side,parent,'c_terminal',{'side'});add(both,parent,'c_terminal',{'side','terminal'})
 return rows

def resolve(parent,context,groups):
 matches=[s for s in templates() if (s['parent'],s['context'],set(s['modified_groups']))==(parent,context,set(groups))]
 if len(matches)!=1:raise ValueError(f'Unsupported chemical state: {parent} / {context} / {sorted(groups)}')
 return matches[0]

def by_name(name):
 return next(s for s in templates() if s['name']==name)
