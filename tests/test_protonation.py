"""Conservative decisions and real-PROPKA boundary checks without remote tools."""
from pathlib import Path
from types import SimpleNamespace as NS
import importlib.util
import tempfile
import unittest
from cationize_protein.protonation import assess_predictions, prepare_input, calculate, ProtonationError


def row(ordinal=2,group='HIS',state='HIE',pka=9.,**extra):
    return dict(residue_ordinal=ordinal,chain_id='B',chain_ordinal=2,original_id='42',
                group=group,atom='CG',source_state=state,pKa=pka,introduced_amine=False,
                modified_residue=False,**extra)


def fixture():
    # A dipeptide with sane geometry; all heavy bonds are explicit.
    atoms=[]; residues=[]; bonds=[]; xyz=[]
    data=[('ALA',[('N',7,(0.,0.,0.)),('CA',6,(1.45,0.,0.)),('C',6,(2.,1.4,0.)),('O',8,(1.4,2.4,0.)),('CB',6,(1.8,-.8,1.2))]),
          ('GLY',[('N',7,(3.30,1.5,0.)),('CA',6,(4.,2.7,0.)),('C',6,(5.5,2.5,0.)),('O',8,(6.,1.4,0.)),('OXT',8,(6.2,3.5,0.))])]
    for ri,(name,rows) in enumerate(data):
        native={}; indices=[]
        for name_a,z,coord in rows:
            i=len(atoms); atoms.append(NS(index=i,name=name_a,atomic_number=z));xyz.append(coord);native[name_a]=i;indices.append(i)
        residues.append(NS(ordinal=ri+1,input_name=name,canonical_name=name,parent=name,context='n_terminal' if ri==0 else 'c_terminal',modified_groups=(),chain_id='original-chain',chain_ordinal=1,original_id=str(42+ri),atom_indices=tuple(indices),native_atom_indices=native,native_atom_names={n:n for n in native},source_spec=None))
        for a,b in [('N','CA'),('CA','C'),('C','O')]+([('CA','CB')] if ri==0 else [('C','OXT')]):
            bonds.append((native[a],native[b]))
    bonds.append((residues[0].native_atom_indices['C'],residues[1].native_atom_indices['N']))
    return NS(atoms=tuple(atoms),residues=tuple(residues),coordinates=xyz,bonds=tuple(bonds),peptide_bonds=(bonds[-1],),chains=((1,2),))


def terminal_template_fixture(name):
    """Place one incoming ALA on a fitted C-terminal product for API checks."""
    import math
    from cationize_protein.off import read_units
    from cationize_protein.catalog import by_name
    from cationize_protein.resources import default_parameters
    spec=by_name(name);u=read_units(default_parameters()/f'{name}.lib')[name]
    def add(a,b):return tuple(x+y for x,y in zip(a,b))
    def sub(a,b):return tuple(x-y for x,y in zip(a,b))
    def mul(a,k):return tuple(x*k for x in a)
    def norm(a):return mul(a,1/math.sqrt(sum(x*x for x in a)))
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    amap={'N1':'N','C1':'CA','C2':'C','O1':'O'} if name=='LET' else {}
    native={amap.get(a['name'],a['name']):i for i,a in enumerate(u['atoms'])}
    n=u['positions'][native['N']]
    neighbors=[]
    for i,j in u['bonds']:
        if i==native['N']:neighbors.append(j)
        if j==native['N']:neighbors.append(i)
    v=norm(tuple(-sum(norm(sub(u['positions'][j],n))[k] for j in neighbors) for k in range(3)))
    seed=(1.,0.,0.) if abs(v[0])<.8 else (0.,1.,0.)
    perpendicular=norm(sub(seed,mul(v,dot(seed,v))))
    q=(v[1]*perpendicular[2]-v[2]*perpendicular[1],v[2]*perpendicular[0]-v[0]*perpendicular[2],v[0]*perpendicular[1]-v[1]*perpendicular[0])
    placements=[]
    for k in range(36):
        angle=2*math.pi*k/36;p=add(mul(perpendicular,math.cos(angle)),mul(q,math.sin(angle)))
        c=add(n,mul(v,1.33));ca=add(c,mul(add(mul(v,.5),mul(p,.8660254)),1.52))
        trial=[add(ca,mul(v,1.45)),ca,c,add(c,mul(sub(mul(v,.5),mul(p,.8660254)),1.23)),add(ca,mul(p,1.53))]
        clearance=min(math.dist(x,y) for i,x in enumerate(trial) for j,y in enumerate(u['positions']) if u['atoms'][j]['atomic_number']!=1 and (i,j)!=(2,native['N']))
        placements.append((clearance,trial))
    coords=max(placements,key=lambda x:x[0])[1]
    atoms=[NS(index=i,name=a,atomic_number=z) for i,(a,z) in enumerate([('N',7),('CA',6),('C',6),('O',8),('CB',6)])]
    offset=len(atoms)
    atoms.extend(NS(index=offset+i,name=a['name'],atomic_number=a['atomic_number']) for i,a in enumerate(u['atoms']))
    coords.extend(u['positions']);parent_names={a.name:a.index for a in atoms[:offset]}
    first=NS(ordinal=1,input_name='ALA',canonical_name='ALA',parent='ALA',context='n_terminal',modified_groups=(),chain_id='B',chain_ordinal=1,original_id='42',atom_indices=tuple(range(offset)),native_atom_indices=parent_names,native_atom_names={n:n for n in parent_names},source_spec=None)
    second=NS(ordinal=2,input_name=name,canonical_name=name,parent=spec['parent'],context='c_terminal',modified_groups=tuple(spec['modified_groups']),chain_id='B',chain_ordinal=1,original_id='43',atom_indices=tuple(range(offset,len(atoms))),native_atom_indices={n:i+offset for n,i in native.items()},native_atom_names={amap.get(a['name'],a['name']):a['name'] for a in u['atoms']},source_spec=spec)
    bonds=[(0,1),(1,2),(2,3),(1,4),(2,offset+native['N'])]+[(i+offset,j+offset) for i,j in u['bonds']]
    return NS(atoms=tuple(atoms),residues=(first,second),coordinates=coords,bonds=bonds,peptide_bonds=((2,offset+native['N']),),chains=((1,2),))


class ProtonationDecisionTests(unittest.TestCase):
    def test_two_controls_and_consecutive_proposals_required(self):
        first=assess_predictions([row()],[row(pka=8.5)])
        self.assertEqual(first['proposed_native_states'],{'2':'HIP'})
        self.assertEqual(first['applicable_native_states'],{})
        second=assess_predictions([row()],[row(pka=8.5)],first)
        self.assertEqual(second['applicable_native_states'],{'2':'HIP'})
        self.assertEqual(second['decisions'][0]['chain_id'],'B')

    def test_borderline_or_disagreeing_controls_reset_persistence(self):
        for pka in [7.9,5.0]:
            result=assess_predictions([row()],[row(pka=pka)],{'2':'HIP'})
            self.assertEqual(result['proposed_native_states'],{})
            self.assertEqual(result['applicable_native_states'],{})
            self.assertEqual(result['decisions'][0]['status'],'held_borderline_or_disagreement')

    def test_histidine_neutral_tautomer_is_retained_not_invented(self):
        for state,target in [('HID','HID'),('HIE','HIE'),('HIP','HIE')]:
            result=assess_predictions([row(state=state,pka=5.)],[row(state=state,pka=5.)],{'2':target})
            self.assertEqual(result['decisions'][0]['target_state'],target)
        result=assess_predictions([row(state='HIP',pka=5.)],[row(state='HIP',pka=5.)])
        self.assertIn('tautomer',result['decisions'][0]['reason'].lower())

    def test_unsupported_nonhistidine_change_is_explicit_and_can_be_capability_enabled(self):
        r=row(group='ASP',state='ASP',pka=9.)
        result=assess_predictions([r],[r],{'2':'ASH'})
        self.assertEqual(result['proposed_native_states'],{'2':'ASH'})
        self.assertEqual(result['applicable_native_states'],{})
        self.assertEqual(result['decisions'][0]['status'],'held_unsupported_state')
        result=assess_predictions([r],[r],{'2':'ASH'},supported_native_states={'ASH'})
        self.assertEqual(result['applicable_native_states'],{'2':'ASH'})

    def test_native_terminal_and_introduced_amine_groups_are_held(self):
        for group,introduced in [('N+',False),('C-',False),('N31',True)]:
            r=row(group=group,pka=1.);r['introduced_amine']=introduced
            result=assess_predictions([r],[r],{'2':'HIE'})
            self.assertEqual(result['proposed_native_states'],{})
            self.assertEqual(result['decisions'][0]['status'],'held_fixed_group')

    def test_missing_control_is_a_hard_error_not_a_silent_skip(self):
        with self.assertRaisesRegex(ProtonationError,'control'):
            assess_predictions([row()],[])

    def test_extra_native_group_in_one_control_is_rejected(self):
        with self.assertRaisesRegex(ProtonationError,'control'):
            assess_predictions([],[row()])

    def test_modified_native_group_change_requires_matching_validated_template(self):
        r=row(state='HET',pka=9.);r['modified_residue']=True;r['native_state']='HIE'
        result=assess_predictions([r],[r],{'2':'HIP'})
        self.assertEqual(result['decisions'][0]['status'],'held_modified_template')
        self.assertEqual(result['applicable_native_states'],{})

    def test_added_amine_without_three_hydrogens_fails_before_propka(self):
        info=fixture();r=info.residues[0];r.canonical_name='GLX';r.modified_groups=('side',)
        with self.assertRaisesRegex(ProtonationError,'ammonium nitrogen'):
            prepare_input(info)

    def test_mapping_uses_global_ordinals_and_preserves_original_identity(self):
        info=fixture(); info.residues[1].chain_ordinal=2;info.residues[1].chain_id='B';info.residues[1].original_id='42';info.chains=((1,),(2,));info.peptide_bonds=()
        text,metadata=prepare_input(info)
        self.assertEqual(metadata['residues'][0]['chain_id'],'original-chain')
        self.assertEqual(metadata['residues'][1]['original_id'],'42')
        lines=[s for s in text.splitlines() if s.startswith('ATOM')]
        self.assertEqual({s[22:26].strip() for s in lines},{'1','2'})
        self.assertEqual({s[21] for s in lines},{'A','B'})

    def test_new_and_reused_templates_keep_added_amine_distinct_from_native_n(self):
        for state,group,amine,nativeN in [('ASX','side','N3','N1'),('GLX','side','S4','N'),('ALT','terminal','T4','N')]:
            info=fixture(); r=info.residues[0]; r.canonical_name=state;r.input_name=state;r.modified_groups=(group,);r.source_spec={'reused_from':'NSL' if state=='ASX' else None}
            old=r.native_atom_indices['N'];info.atoms[old].name=nativeN;r.native_atom_names['N']=nativeN
            # All added ammonium groups have exactly one carbon and three hydrogens.
            indices=list(r.atom_indices);at=list(info.atoms);xyz=list(info.coordinates);edges=list(info.bonds)
            root=len(at);at.append(NS(index=root,name=amine,atomic_number=7));xyz.append((0.,0.,5.));indices.append(root)
            for n,z in [('Cnew',6),('Hn1',1),('Hn2',1),('Hn3',1)]:
                i=len(at);at.append(NS(index=i,name=n,atomic_number=z));xyz.append((0.,0.,6.+i));indices.append(i);edges.append((root,i))
            r.atom_indices=tuple(indices);info.atoms=tuple(at);info.coordinates=xyz;info.bonds=tuple(edges)
            text,metadata=prepare_input(info)
            self.assertEqual(metadata['introduced_amines'][0]['atom'],amine)
            self.assertEqual(metadata['introduced_amines'][0]['fixed_charge_e'],1)
            self.assertEqual(metadata['native_n_termini'][0]['atom'],'N')


@unittest.skipUnless(importlib.util.find_spec('propka'),'PROPKA not installed in this interpreter')
class RealPropkaTests(unittest.TestCase):
    def test_detailed_propka_log_stays_in_file(self):
        from cationize_protein.protonation import calculate
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'quiet'
            with self.assertNoLogs(level='INFO'):
                calculate(terminal_template_fixture('HET'),out)
            self.assertGreater((out/'propka.log').stat().st_size,0)

    def test_real_modified_n_terminal_groups_are_both_recognized(self):
        from cationize_protein.off import read_units
        from cationize_protein.catalog import by_name
        from cationize_protein.resources import default_parameters
        for name in ['ASX','GLX']:
            spec=by_name(name);unit=read_units(default_parameters()/f'{name}.lib')[name]
            at=[NS(index=i,name=a['name'],atomic_number=a['atomic_number']) for i,a in enumerate(unit['atoms'])]
            rename={'N1':'N','C1':'CA','C2':'C','O1':'O','C3':'CB','C4':'CG','O2':'OD1'} if name=='ASX' else {}
            native={rename.get(a.name,a.name):i for i,a in enumerate(at)}
            r=NS(ordinal=1,input_name=name,canonical_name=name,parent=spec['parent'],context='n_terminal',
                 modified_groups=('side',),chain_id='X',chain_ordinal=1,original_id='42',
                 atom_indices=tuple(range(len(at))),native_atom_indices=native,
                 native_atom_names={rename.get(a.name,a.name):a.name for a in at},source_spec=spec)
            info=NS(atoms=tuple(at),residues=(r,),coordinates=unit['positions'],bonds=unit['bonds'],peptide_bonds=(),chains=((1,),))
            with self.subTest(name=name),tempfile.TemporaryDirectory() as d:
                standard=calculate(info,Path(d)/'standard')
                fixed=calculate(info,Path(d)/'fixed',fixed=True)
                self.assertEqual(standard['recognized_added_amines'],1)
                self.assertEqual(fixed['recognized_added_amines'],1)
                self.assertIn('N+',{g['group'] for g in standard['predictions']})
                self.assertIn('N31',{g['group'] for g in standard['predictions']})
                self.assertNotIn('N31',{g['group'] for g in fixed['predictions']})

    def test_real_terminal_products_keep_amines_and_native_side_groups(self):
        expected={'LET':None,'HET':'HIS','HDT':'HIS','HPT':'HIS','AST':'ASP','ASZ':None,'ASD':None,'GLT':'GLU','GLZ':None,'GLD':None}
        for name,native_group in expected.items():
            info=terminal_template_fixture(name)
            with self.subTest(name=name),tempfile.TemporaryDirectory() as d:
                result=calculate(info,Path(d)/'standard')
                self.assertEqual(result['recognized_added_amines'],len(info.residues[1].modified_groups))
                if native_group:self.assertIn(native_group,{r['group'] for r in result['predictions']})
                fixed=calculate(info,Path(d)/'fixed',fixed=True)
                self.assertEqual(fixed['recognized_added_amines'],result['recognized_added_amines'])

    def test_real_inferred_extra_bond_blocks_calculation(self):
        info=fixture();info.coordinates[-1]=(1.5,1.,.5)
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ProtonationError,'graph differs'):
                calculate(info,Path(d)/'standard')
            self.assertTrue((Path(d)/'standard'/'graph_mismatch.json').is_file())

    def test_real_package_calculation_checks_graph_and_preserves_source(self):
        info=fixture(); before=list(info.coordinates)
        with tempfile.TemporaryDirectory() as d:
            result=calculate(info,Path(d)/'standard')
            self.assertTrue(result['exact_heavy_graph_match'])
            self.assertIn('N+',{r['group'] for r in result['predictions']})
            self.assertIn('C-',{r['group'] for r in result['predictions']})
            self.assertEqual(info.coordinates,before)
