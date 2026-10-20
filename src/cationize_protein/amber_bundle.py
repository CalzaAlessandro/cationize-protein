"""Complete per-system parameters, independent of an installed Amber data tree.

The physical XML/JSON and parm7 retain every term, including CMAP. The ordinary
Amber frcmod is accompanied by explicit atom types, charges and full-system OFF
connectivity. A Python rebuild preserves CMAP even when a LEaP version cannot
read a residue-based CMAP file or recognize modified backbone atom names.
"""
import copy
import csv
import io
import json
from pathlib import Path
import string
import numpy as np


def _retype(structure,source_types=None):
    """Keep distinct force-field types; LJ equality alone is insufficient."""
    s=copy.deepcopy(structure)
    source_types=source_types or [f'atom:{i}' for i in range(len(s.atoms))]
    if len(source_types)!=len(s.atoms):raise ValueError('Incomplete source type identity')
    codes={};alphabet=string.ascii_letters+string.digits
    for a,source in zip(s.atoms,source_types):
        key=(source,float(a.mass),float(a.rmin),float(a.epsilon))
        if key not in codes:
            i=len(codes)
            if i>=26*len(alphabet):raise ValueError('Two-character Amber export type capacity exceeded')
            codes[key]=string.ascii_lowercase[i//len(alphabet)]+alphabet[i%len(alphabet)]
        a.type=codes[key]
        a.atom_type=copy.copy(a.atom_type);a.atom_type.name=a.type
    return s,source_types


def _complete_torsions(structure,parameters,system):
    """LEaP needs an entry even when the physical torsion amplitude is zero."""
    from parmed import DihedralType,AngleType
    from .exports import _one_four_scaling
    coulomb,lj=_one_four_scaling(system)
    graph={a.idx:set() for a in structure.atoms}
    for bond in structure.bonds:
        i,j=bond.atom1.idx,bond.atom2.idx;graph[i].add(j);graph[j].add(i)
    added=0
    for j,neighbors in graph.items():
        for i in neighbors:
            for k in neighbors-{i}:
                key=tuple(structure.atoms[a].type for a in (i,j,k))
                if key not in parameters.angle_types and key[::-1] not in parameters.angle_types:
                    # An absent OpenMM angle has zero energy; LEaP still needs
                    # an explicit entry before it can rebuild the topology.
                    parameters.angle_types[key]=AngleType(0.,109.5)
    for j in graph:
        for k in graph[j]:
            for i in graph[j]-{k}:
                for l in graph[k]-{j,i}:
                    key=tuple(structure.atoms[a].type for a in (i,j,k,l))
                    if key not in parameters.dihedral_types and key[::-1] not in parameters.dihedral_types:
                        parameters.dihedral_types[key]=DihedralType(0.,1,0.,scee=1/coulomb,scnb=1/lj);added+=1
    return added,graph


def _hybridization(atom,graph,atoms):
    if atom.atomic_number==6:return 'sp2' if len(graph[atom.idx])==3 else 'sp3'
    if atom.atomic_number==7:
        if len(graph[atom.idx])>=4:return 'sp3'
        for i in graph[atom.idx]:
            if atoms[i].atomic_number==6 and len(graph[i])==3:return 'sp2'
        return 'sp3'
    if atom.atomic_number==8:return 'sp2' if len(graph[atom.idx])==1 else 'sp3'
    return 'sp3'


def _frcmod(parameters):
    """Conventional fixed-width frcmod; full precision is retained in JSON/XML."""
    p=parameters;lines=['Complete physical parameters; CMAP is in PARAMETERS.json/system.xml','MASS']
    for name,t in sorted(p.atom_types.items()):lines.append(f'{name:<6}{t.mass:6.3f}')
    for heading,types,fields in [('BOND',p.bond_types,('k','req')),('ANGLE',p.angle_types,('k','theteq'))]:
        lines+=['',heading];done=set()
        for key,t in sorted(types.items()):
            identity=min(key,key[::-1])
            if identity in done:continue
            done.add(identity)
            lines.append('-'.join(f'{n:<2}' for n in key)+f'   {getattr(t,fields[0]):8.3f}  {getattr(t,fields[1]):6.3f}')
    lines+=['','DIHE'];done=set()
    for key,t in sorted(p.dihedral_types.items()):
        identity=min(key,key[::-1])
        if identity in done:continue
        done.add(identity);terms=[t] if hasattr(t,'phi_k') else list(t)
        for i,d in enumerate(terms):
            period=d.per*(-1 if i<len(terms)-1 else 1)
            lines.append('-'.join(f'{n:<2}' for n in key)+f' {1:4d} {d.phi_k:14.8f} {d.phase:8.3f} {period:5.1f} SCEE={d.scee} SCNB={d.scnb}')
    lines+=['','IMPROPER'];done=set()
    for key,t in sorted(p.improper_periodic_types.items()):
        # Amber improper central atom is the third entry; other three commute.
        identity=(key[2],tuple(sorted((key[0],key[1],key[3]))))
        if identity in done:continue
        done.add(identity)
        lines.append('-'.join(f'{n:<2}' for n in key)+f' {t.phi_k:14.8f} {t.phase:8.3f} {t.per:5.1f}')
    lines+=['','NONB']
    for name,t in sorted(p.atom_types.items()):lines.append(f'{name:<2}  {t.rmin:12.8f} {t.epsilon:12.8f}')
    return '\n'.join(lines)+'\n\n'


def _off(structure):
    import parmed
    container=parmed.modeller.ResidueTemplateContainer.from_structure(structure,term_decorate=False)
    container.name='protein'
    # Connectivity is serialized globally below, including chain/disulfide links.
    for r in container:r.head=r.tail=None;r.connections=[]
    stream=io.StringIO();parmed.amber.AmberOFFLibrary.write({'protein':container},stream)
    lines=stream.getvalue().splitlines();start=next(i for i,l in enumerate(lines) if l.startswith('!entry.protein.unit.connectivity '))+1
    end=next(i for i in range(start,len(lines)) if lines[i].startswith('!entry.'))
    lines[start:end]=[f' {b.atom1.idx+1} {b.atom2.idx+1} 1' for b in structure.bonds]
    start=next(i for i,l in enumerate(lines) if l.startswith('!entry.protein.unit.positions '))+1
    end=next(i for i in range(start,len(lines)) if lines[i].startswith('!entry.'))
    lines[start:end]=[' '+' '.join(f'{v:.15g}' for v in xyz) for xyz in structure.coordinates]
    return '\n'.join(lines)+'\n'


def write_bundle(structure,system,directory,force_field,source_types=None):
    """Write complete ordinary frcmod plus all physical parameters and rebuild."""
    import parmed
    from openmm import unit
    directory=Path(directory)
    s,source_types=_retype(structure,source_types)
    from .exports import _one_four_scaling
    coulomb_scale,lj_scale=_one_four_scaling(system)
    # ParmEd's OpenMM reader stores 1-4 interactions as explicit adjustments.
    # Its bare torsion types may retain SCEE/SCNB=1 despite scaled physical
    # exceptions; copy the actual uniform scales before printing frcmod/JSON.
    for dihedral in s.dihedrals:
        for term in ([dihedral.type] if hasattr(dihedral.type,'phi_k') else dihedral.type):
            term.scee=1/coulomb_scale;term.scnb=1/lj_scale
    parameters=parmed.amber.AmberParameterSet.from_structure(s)
    zero_terms,graph=_complete_torsions(s,parameters,system)
    (directory/'parameters.frcmod').write_text(_frcmod(parameters))
    (directory/'protein.lib').write_text(_off(s))
    atoms=[dict(index=i,name=a.name,residue=a.residue.name,residue_index=a.residue.idx,
                residue_number=a.residue.number,chain=a.residue.chain,insertion_code=a.residue.insertion_code,
                atomic_number=a.atomic_number,mass_Da=float(a.mass),charge_e=float(a.charge),
                PB_radius_A=float(a.solvent_radius),GB_screen=float(a.screen),
                amber_export_type=a.type,source_force_field_type=source_types[i],
                LJ_half_Rmin_A=float(a.rmin),LJ_epsilon_kcal_mol=float(a.epsilon)) for i,a in enumerate(s.atoms)]
    rows=lambda objects,atom_names,attributes:[dict(atoms=[getattr(v,n).idx for n in atom_names],
         **{name:float(getattr(v.type,attribute)) for name,attribute in attributes}) for v in objects]
    record=dict(schema_version=1,force_field=force_field,atom_indices='zero based; same order as system.parm7/XML',
        units='Angstrom, degree, elementary charge, dalton, kcal/mol; harmonic Amber terms k*(x-x0)^2',
        atoms=atoms,coordinates_A=np.asarray(s.coordinates).tolist(),
        bonds=rows(s.bonds,('atom1','atom2'),(('k_kcal_mol_A2','k'),('equilibrium_A','req'))),
        angles=rows(s.angles,('atom1','atom2','atom3'),(('k_kcal_mol_rad2','k'),('equilibrium_degrees','theteq'))),
        torsions=[],exceptions=[],cmap_terms=[],cmap_maps=[])
    for d in s.dihedrals:
        for t in ([d.type] if hasattr(d.type,'phi_k') else d.type):
            record['torsions'].append(dict(atoms=[a.idx for a in (d.atom1,d.atom2,d.atom3,d.atom4)],
                improper=d.improper,ignore_14=d.ignore_end,k_kcal_mol=float(t.phi_k),
                periodicity=int(t.per),phase_degrees=float(t.phase),SCEE=float(t.scee),SCNB=float(t.scnb)))
    for f in system.getForces():
        if type(f).__name__=='NonbondedForce':
            for i in range(f.getNumExceptions()):
                a,b,q,sigma,epsilon=f.getExceptionParameters(i)
                record['exceptions'].append(dict(atoms=[a,b],charge_product_e2=q.value_in_unit(unit.elementary_charge**2),
                    sigma_A=sigma.value_in_unit(unit.angstrom),epsilon_kcal_mol=epsilon.value_in_unit(unit.kilocalories_per_mole)))
    for cm in s.cmaps:
        record['cmap_terms'].append(dict(atoms=[a.idx for a in (cm.atom1,cm.atom2,cm.atom3,cm.atom4,cm.atom5)],map_index=cm.type.idx))
    for cm in s.cmap_types:
        record['cmap_maps'].append(dict(index=cm.idx,resolution=cm.resolution,grid_kcal_mol=list(cm.grid),
            convention='ParmEd/Amber grid; angles -180 to 180 degrees; first index times resolution plus second index'))
    (directory/'PARAMETERS.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    with (directory/'ATOM_TYPES.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(atoms[0]));writer.writeheader();writer.writerows(atoms)
    # Rebuild uses complete numeric data; no installed Amber parameter library.
    (directory/'rebuild.py').write_text('''"""Rebuild Amber files with OpenMM+ParmEd; no Amber installation needed."""
import json,sys
from pathlib import Path
import numpy as np
import parmed as pmd
from openmm import XmlSerializer,unit
root=Path(__file__).resolve().parent
out=Path(sys.argv[1]) if len(sys.argv)>1 else root/'rebuilt'
out.mkdir(parents=True,exist_ok=False)
record=json.loads((root/'PARAMETERS.json').read_text())
s=pmd.Structure()
for a in record['atoms']:
    atom=pmd.Atom(name=a['name'],atomic_number=a['atomic_number'],mass=a['mass_Da'])
    s.add_atom(atom,a['residue'],a['residue_index']+1,chain=a['chain'],inscode=a['insertion_code'])
for b in record['bonds']:
    i,j=b['atoms'];s.bonds.append(pmd.Bond(s.atoms[i],s.atoms[j]))
system=XmlSerializer.deserialize((root/'system.xml').read_text())
typed=pmd.openmm.load_topology(s.topology,system,xyz=np.asarray(record['coordinates_A'])*unit.angstrom)
for atom,source in zip(typed.atoms,record['atoms']):
    atom.solvent_radius=source['PB_radius_A'];atom.screen=source['GB_screen']
parm=pmd.amber.AmberParm.from_structure(typed)
parm.save(str(out/'system.parm7'),overwrite=False)
parm.save(str(out/'system.rst7'),format='rst7',overwrite=False)
print('Rebuilt complete physical parameters, including CMAP:',out)
''')
    # A pure LEaP rebuild is appropriate only when there is no CMAP to lose.
    elements={1:'H',6:'C',7:'N',8:'O',16:'S'}
    atom_type_rows={}
    for atom in s.atoms:
        entry=(elements[atom.atomic_number],_hybridization(atom,graph,s.atoms))
        if atom.type in atom_type_rows and atom_type_rows[atom.type]!=entry:
            raise ValueError('Ambiguous hybridization for export type '+atom.type)
        atom_type_rows[atom.type]=entry
    declarations=['addAtomTypes {']+[f' {{ "{name}" "{element}" "{hybrid}" }}' for name,(element,hybrid) in sorted(atom_type_rows.items())]+['}']
    commands=declarations+['loadamberparams parameters.frcmod','loadoff protein.lib',
              'saveamberparm protein rebuilt.parm7 rebuilt.rst7','quit']
    (directory/'rebuild.in').write_text('\n'.join(commands)+'\n' if not s.cmaps else
        '# CMAP is fully retained in system.parm7, system.xml and PARAMETERS.json.\n'
        '# Use: python rebuild.py NEW_DIRECTORY\n# A conventional frcmod alone would omit CMAP.\nquit\n')
    (directory/'PARAMETERS_README.md').write_text(f'''# Complete {force_field} parameters

`PARAMETERS.json` prints every atom charge/mass/LJ term, bond, angle, torsion,
1-4/exclusion entry, and CMAP map/index at full numeric precision. `ATOM_TYPES.csv`
maps exported two-character types to their original force-field identities.
`parameters.frcmod` contains every ordinary Amber parameter used by this system;
`protein.lib` contains all atoms, residue membership, charges, coordinates and
connectivity, including peptide bonds and disulfides. These are system-specific
files, not modifications to the original fitted library.

`system.parm7`/`system.rst7` are ready physical Amber inputs when Amber export is
requested. No source leaprc or installed Amber parameter data is needed to read
them. `python rebuild.py NEW_DIRECTORY` reconstructs them from `system.xml` and
the complete JSON with OpenMM and ParmEd. For systems without CMAP, `tleap -s -f
rebuild.in` is also provided; run it in a copy of this directory.

CMAP terms present: {len(s.cmaps)}. All are retained in XML/JSON/parm7. The ordinary
frcmod alone does not encode them. An Amber engine must support the physical terms
in its input: ff19SB requires a CMAP-capable version (Amber20+ recommended;
Amber18/AmberTools19 needed the authors' patch). No parameter dump can add CMAP
support to an older executable. Numerical engine/version tests remain separate
from file conversion checks. No force constants or RESP charges are refitted here.
''')
    return dict(full_parameters='PARAMETERS.json',library='protein.lib',ordinary_parameters='parameters.frcmod',
                printed_CMAP_terms=len(s.cmaps),explicit_zero_torsion_types=zero_terms,source_Amber_data_tree_required=False,
                type_mapping='ATOM_TYPES.csv',CMAP_rebuild='rebuild.py')
