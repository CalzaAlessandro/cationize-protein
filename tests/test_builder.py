"""Real ParmEd graph fixtures catch coordinate drift and lost crosslinks."""
import unittest
import numpy as np
try:
    import parmed as pmd
    from parmed.modeller import ResidueTemplate
except ImportError:
    pmd = None
from cationize_protein.builder import replace_residue, export_mutation_pdb


def atom(name, xyz, z=6, charge=0.):
    a = pmd.Atom(name=name, atomic_number=z, mass=float(z * 2), charge=charge, type='CT')
    a.xx, a.xy, a.xz = xyz
    return a


def fixtures(legacy=False):
    s = pmd.Structure()
    rows = [
        ('ALA', 5, 'A', [('C', (-1., 0, 0), 6), ('O', (-1., 1, 0), 8)]),
        ('ASP', 6, 'A', [('N', (0., 0, 0), 7), ('CA', (1., 0, 0), 6),
                         ('C', (1., 1, 0), 6), ('O', (1., 2, 0), 8),
                         ('CB', (2., 0, 0), 6), ('CG', (2., 0, 1), 6),
                         ('OD1', (2., 1, 1), 8), ('OD2', (3., 0, 1), 8)]),
        ('CYX', 7, 'A', [('N', (2., 1, 0), 7), ('CA', (3., 1, 0), 6),
                         ('C', (3., 2, 0), 6), ('O', (3., 3, 0), 8),
                         ('CB', (4., 1, 0), 6), ('SG', (4., 1, 1), 16)]),
        ('CYX', 7, 'B', [('SG', (4., 1, 3), 16)]),
    ]
    for rn, number, chain, atoms in rows:
        for an, xyz, z in atoms:
            s.add_atom(atom(an, xyz, z, charge=.123), rn, number, chain=chain)
    def bond(r1, a1, r2, a2, order=1.):
        d1 = {a.name:a for a in s.residues[r1].atoms}; d2 = {a.name:a for a in s.residues[r2].atoms}
        s.bonds.append(pmd.Bond(d1[a1], d2[a2], order=order))
    bond(0,'C',0,'O',2.); bond(0,'C',1,'N'); bond(1,'N',1,'CA')
    for a,b in [('CA','C'),('C','O'),('CA','CB'),('CB','CG'),('CG','OD1'),('CG','OD2')]:
        bond(1,a,1,b)
    bond(1,'C',2,'N')
    for a,b in [('N','CA'),('CA','C'),('C','O'),('CA','CB'),('CB','SG')]:
        bond(2,a,2,b)
    bond(2,'SG',3,'SG'); s.residues[2].ter=True
    s.box = [30.,31.,32.,90.,90.,90.]
    t = ResidueTemplate('ASC')
    coords = [('N',(0.,0,0),7),('CA',(1.,0,0),6),('C',(1.,1,0),6),('O',(1.,2,0),8),
              ('CB',(2.,0,0),6),('CG',(3.,0,0),6),('OD1',(3.,1,0),8),
              ('S1',(4.33,0,0),7),('S2',(4.33,0,1.),6)]
    mapping = dict(N='N1',CA='C1',C='C2',O='O1',CB='C3',CG='C4',OD1='O2') if legacy else None
    for name, xyz, z in coords:
        t.add_atom(atom(mapping.get(name,name) if mapping else name,xyz,z,charge=.456))
    for a,b in [('N','CA'),('CA','C'),('C','O'),('CA','CB'),('CB','CG'),('CG','OD1'),('CG','S1'),('S1','S2')]:
        t.add_bond(mapping.get(a,a) if mapping else a,mapping.get(b,b) if mapping else b)
    return s,t,mapping


def graph(structure):
    return {tuple(sorted(((b.atom1.residue.idx,b.atom1.name),(b.atom2.residue.idx,b.atom2.name)))):b.order
            for b in structure.bonds}


@unittest.skipIf(pmd is None,'ParmEd not installed in this interpreter')
class BuilderTests(unittest.TestCase):
    def test_inherited_and_outside_coordinates_and_external_graph_are_exact(self):
        s,t,_=fixtures(); before=s.coordinates.copy(); old_graph=graph(s)
        result=replace_residue(s,1,t)
        self.assertIsNot(result,s)
        self.assertEqual(result.residues[1].name,'ASC')
        self.assertNotIn('OD2',[a.name for a in result.residues[1].atoms])
        for ri in range(len(s.residues)):
            old={a.name:a for a in s.residues[ri].atoms}
            for a in result.residues[ri].atoms:
                if a.name in old:
                    np.testing.assert_array_equal([a.xx,a.xy,a.xz],[old[a.name].xx,old[a.name].xy,old[a.name].xz])
                    if ri!=1:self.assertEqual(a.charge,old[a.name].charge)
        for edge,order in old_graph.items():
            if all(ri!=1 for ri,name in edge) or edge[0][0]!=edge[1][0]:
                self.assertEqual(graph(result)[edge],order)
        self.assertTrue(result.residues[2].ter)
        self.assertEqual([r.chain for r in result.residues],['A','A','A','B'])
        np.testing.assert_array_equal(result.box,s.box)
        np.testing.assert_array_equal(s.coordinates,before)
        self.assertEqual(graph(s),old_graph)
        self.assertFalse(result.angles)

    def test_new_component_uses_local_carbonyl_frame_not_global_backbone_fit(self):
        s,t,_=fixtures()
        result=replace_residue(s,1,t)
        atoms={a.name:a for a in result.residues[1].atoms}
        np.testing.assert_allclose([atoms['S1'].xx,atoms['S1'].xy,atoms['S1'].xz],[2.,0.,2.33],atol=1e-12)
        np.testing.assert_allclose([atoms['S2'].xx,atoms['S2'].xy,atoms['S2'].xz],[1.,0.,2.33],atol=1e-12)
        self.assertAlmostEqual(np.linalg.norm(result.coordinates[atoms['S1'].idx]-result.coordinates[atoms['CG'].idx]),1.33)
        self.assertEqual(atoms['S1'].charge,.456)

    def test_legacy_atom_map_reconnects_original_peptide_bonds(self):
        s,t,mapping=fixtures(legacy=True)
        result=replace_residue(s,1,t,mapping)
        self.assertIn(((0,'C'),(1,'N1')),graph(result))
        self.assertIn(((1,'C2'),(2,'N')),graph(result))
        target={a.name:a for a in result.residues[1].atoms}
        for old,new in mapping.items():
            original=next(a for a in s.residues[1].atoms if a.name==old)
            np.testing.assert_array_equal([target[new].xx,target[new].xy,target[new].xz],[original.xx,original.xy,original.xz])

    def test_deleted_covalent_attachment_is_rejected(self):
        s,t,_=fixtures()
        od2=next(a for a in s.residues[1].atoms if a.name=='OD2')
        s.bonds.append(pmd.Bond(od2,s.residues[3].atoms[0]))
        with self.assertRaisesRegex(ValueError,'attachment'):
            replace_residue(s,1,t)

    def test_multiple_attachment_new_component_and_collinear_frames_are_rejected(self):
        s,t,_=fixtures();t.add_bond('S1','C')
        with self.assertRaisesRegex(ValueError,'attachment'):
            replace_residue(s,1,t)
        s,t,_=fixtures()
        for a in t.atoms:a.xy=a.xz=0.
        with self.assertRaisesRegex(ValueError,'frame|collinear'):
            replace_residue(s,1,t)

    def test_invalid_indices_and_atom_maps_are_rejected(self):
        s,t,_=fixtures()
        for ri in (-1,100,True,1.5):
            with self.subTest(ri=ri):
                with self.assertRaises(ValueError):replace_residue(s,ri,t)
        for mapping in ({'N':'CA','CA':'CA'},{'doesnotexist':'N'},{'N':'absent'}):
            with self.subTest(mapping=mapping):
                with self.assertRaises(ValueError):replace_residue(s,1,t,mapping)

    def test_export_has_every_bond_chain_break_and_full_precision_coordinates(self):
        s,t,_=fixtures();result=replace_residue(s,1,t)
        pdb,info=export_mutation_pdb(result)
        atom_lines=[line for line in pdb.splitlines() if line[:6] in ('ATOM  ','HETATM')]
        self.assertEqual(len(atom_lines),len(result.atoms))
        self.assertEqual([line[17:20].strip() for line in atom_lines], [a.residue.name for a in result.atoms])
        self.assertEqual(sum(line.startswith('CONECT') for line in pdb.splitlines()),len(result.bonds))
        self.assertGreaterEqual(sum(line.startswith('TER') for line in pdb.splitlines()),2)
        self.assertEqual(len(info['interresidue_bonds']),3)
        self.assertEqual(len(info['bonds']),len(result.bonds))
        np.testing.assert_array_equal(info['coordinates_A'],result.coordinates)
        self.assertEqual(info['residues'][2]['source_number'],7)
        self.assertEqual(info['residues'][3]['chain'],'B')
        self.assertTrue(info['requires_parameterization'])

    def test_dual_modifications_place_each_added_component_about_its_own_carbonyl(self):
        s,t,_=fixtures()
        t.add_atom(atom('T1',(2.33,1.,0.),7));t.add_atom(atom('T2',(2.33,1.,1.),6))
        t.add_bond('C','T1');t.add_bond('T1','T2')
        result=replace_residue(s,1,t)
        a={atom.name:atom for atom in result.residues[1].atoms}
        np.testing.assert_allclose([a['S1'].xx,a['S1'].xy,a['S1'].xz],[2.,0.,2.33],atol=1e-12)
        np.testing.assert_allclose([a['T1'].xx,a['T1'].xy,a['T1'].xz],[2.33,1.,0.],atol=1e-12)
        self.assertAlmostEqual(np.linalg.norm(result.coordinates[a['T1'].idx]-result.coordinates[a['C'].idx]),1.33)

    def test_replacing_disulfide_residue_retains_its_sulfur_crosslink(self):
        from copy import copy
        s,_,_=fixtures();t=ResidueTemplate('CYT')
        for a in s.residues[2].atoms:t.add_atom(copy(a))
        for b in s.bonds:
            if b.atom1.residue is s.residues[2] and b.atom2.residue is s.residues[2]:
                t.add_bond(b.atom1.name,b.atom2.name,order=b.order)
        t.add_atom(atom('T1',(4.33,2.,0.),7));t.add_bond('C','T1')
        result=replace_residue(s,2,t)
        self.assertIn(((2,'SG'),(3,'SG')),graph(result))
        original=next(a for a in s.residues[2].atoms if a.name=='SG')
        new=next(a for a in result.residues[2].atoms if a.name=='SG')
        np.testing.assert_array_equal([original.xx,original.xy,original.xz],[new.xx,new.xy,new.xz])

    def test_pdb_fields_roundtrip_with_real_parmed_without_amber(self):
        import io,json
        s,t,_=fixtures();result=replace_residue(s,1,t)
        text,info=export_mutation_pdb(result)
        parsed=pmd.formats.PDBFile.parse(io.StringIO(text),skip_bonds=True)
        self.assertEqual([a.name for a in parsed.atoms],[a.name for a in result.atoms])
        self.assertEqual([r.name for r in parsed.residues],['ALA','ASC','CYX','CYX'])
        self.assertEqual([r.chain for r in parsed.residues],['A','A','A','B'])
        np.testing.assert_allclose(parsed.coordinates,result.coordinates,atol=.00051)
        self.assertEqual(json.loads(json.dumps(info))['interresidue_bonds'],info['interresidue_bonds'])

    def test_nonfinite_coordinates_and_pdb_field_overflow_are_rejected(self):
        s,t,_=fixtures();t.atoms[-1].xx=float('nan')
        with self.assertRaisesRegex(ValueError,'finite'):replace_residue(s,1,t)
        s,t,_=fixtures();s.residues[0].name='TOOLONG'
        with self.assertRaisesRegex(ValueError,'field width'):export_mutation_pdb(s)

    def test_unresolved_inherited_off_elements_and_zero_masses_use_native_metadata(self):
        s,t,_=fixtures()
        for a in t.atoms:a.mass=0.
        unknown=next(a for a in t.atoms if a.name=='CB');unknown.atomic_number=-1
        result=replace_residue(s,1,t)
        byname={a.name:a for a in result.residues[1].atoms}
        self.assertEqual(byname['CB'].atomic_number,6)
        self.assertEqual(byname['CB'].mass,next(a.mass for a in s.residues[1].atoms if a.name=='CB'))
        self.assertGreater(byname['S1'].mass,0.)
        self.assertEqual(unknown.atomic_number,-1)
        self.assertEqual(unknown.mass,0.)

    def test_unknown_new_element_cannot_be_guessed_from_an_atom_name(self):
        s,t,_=fixtures();next(a for a in t.atoms if a.name=='S1').atomic_number=-1
        with self.assertRaisesRegex(ValueError,'element'):
            replace_residue(s,1,t)


if __name__=='__main__':unittest.main()
