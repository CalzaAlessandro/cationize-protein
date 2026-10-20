"""Native state changes preserve the molecular graph and fixed product chemistry."""
from types import SimpleNamespace as NS
import importlib.util
import unittest
from cationize_protein.native_states import plan_native_states, UnsupportedNativeState


def topology(names):
    residues=[NS(ordinal=i+1,canonical_name=n,parent=n,modified_groups=(),chain_id='B',original_id=str(42+i)) for i,n in enumerate(names)]
    return NS(residues=tuple(residues),chains=(tuple(range(1,len(names)+1)),),peptide_bonds=(),disulfide_bonds=(),atoms=())


class NativeStatePlanTests(unittest.TestCase):
    def test_histidine_change_and_acid_change_have_physical_formal_charge_delta(self):
        info=topology(['HIE','ASP','LYS','CYS','TYR'])
        plan=plan_native_states(info,{'1':'HIP','2':'ASH','3':'LYN','4':'CYM','5':'TYM'})
        self.assertEqual(plan['expected_charge_change_e'],-1)
        self.assertEqual(plan['changes'][0]['chain_id'],'B')
        self.assertEqual(plan['changes'][0]['original_id'],'42')

    def test_invalid_residue_and_cross_family_state_are_rejected(self):
        for changes in [{'0':'HIP'},{'x':'HIP'},{'1':'ASH'},{'1':'HIS'},{True:'HIP'}]:
            with self.subTest(changes=changes),self.assertRaises(UnsupportedNativeState):
                plan_native_states(topology(['HIE']),changes)

    def test_modified_product_is_held_until_matching_fitted_template_is_supported(self):
        info=topology(['HET']);info.residues[0].parent='HIE';info.residues[0].modified_groups=('terminal',)
        with self.assertRaisesRegex(UnsupportedNativeState,'modified'):
            plan_native_states(info,{'1':'HIP'})

    def test_noop_requests_are_not_applied(self):
        self.assertEqual(plan_native_states(topology(['HID']),{'1':'HID'})['changes'],[])


_RUNTIME=all(importlib.util.find_spec(n) for n in ('numpy','parmed','openmm'))


def native_peptide(state):
    """Use fitted source geometry only; assign every native charge afresh."""
    import numpy as np
    import parmed
    from openmm import Platform,unit
    from openmm.app import Modeller
    from cationize_protein.parameter_assignment import parameterize
    from cationize_protein.resources import default_parameters
    from test_protonation import terminal_template_fixture
    lib={'HIE':'HET','ASP':'AST','GLU':'GLT','LYS':'LYT','CYS':'CYT','TYR':'TYT'}[state]
    info=terminal_template_fixture(lib)
    # Transform ALA-X(terminal product) into ALA-X-GLY; the native target is
    # internal because ff14SB lacks terminal ASH/GLH/LYN/CYM templates.
    s=parmed.Structure();lookup={}
    for r in info.residues:
        for i in r.atom_indices:
            a=info.atoms[i]
            if r.ordinal==2 and (a.name.startswith('T') or a.name.startswith('HT')):continue
            new=parmed.Atom(name=a.name,atomic_number=a.atomic_number,mass=float(a.atomic_number*2))
            new.xx,new.xy,new.xz=info.coordinates[i]
            s.add_atom(new,'ALA' if r.ordinal==1 else ('HIS' if state=='HIE' else state),r.ordinal,chain='B')
            lookup[i]=new
    for i,j in info.bonds:
        if i in lookup and j in lookup:s.bonds.append(parmed.Bond(lookup[i],lookup[j]))
    r=s.residues[-1];names={a.name:a for a in r.atoms}
    xyz=lambda a:np.array([a.xx,a.xy,a.xz])
    c=xyz(names['C']);ca=xyz(names['CA']);o=xyz(names['O'])
    v=-((ca-c)/np.linalg.norm(ca-c)+(o-c)/np.linalg.norm(o-c));v/=np.linalg.norm(v)
    seed=np.array([1.,0,0]) if abs(v[0])<.8 else np.array([0.,1.,0]);p=seed-np.dot(seed,v)*v;p/=np.linalg.norm(p)
    n=c+1.33*v;gca=n+1.45*(.5*v+.8660254*p);gc=gca+1.52*v
    added={}
    for name,z,coord in [('N',7,n),('CA',6,gca),('C',6,gc),('O',8,gc+1.23*(.5*v+.8660254*p)),('OXT',8,gc+1.25*(.5*v-.8660254*p))]:
        a=parmed.Atom(name=name,atomic_number=z,mass=float(z*2));a.xx,a.xy,a.xz=coord;s.add_atom(a,'GLY',3,chain='B');added[name]=a
    for a,b in [('N','CA'),('CA','C'),('C','O'),('C','OXT')]:s.bonds.append(parmed.Bond(added[a],added[b]))
    s.bonds.append(parmed.Bond(names['C'],added['N']))
    m=Modeller(s.topology,np.asarray(s.coordinates)*unit.angstrom)
    m.addHydrogens(forcefield=None,pH=7.4,platform=Platform.getPlatformByName('Reference'))
    rebuilt=parmed.openmm.load_topology(m.topology,xyz=m.positions.value_in_unit(unit.angstrom))
    rebuilt.residues[1].name=state
    for r in rebuilt.residues:r.chain='B';r.number+=40
    typed,_=parameterize(rebuilt,default_parameters())
    return typed


@unittest.skipUnless(_RUNTIME,'OpenMM/ParmEd/numpy not installed in this interpreter')
class RealNativeStateTests(unittest.TestCase):
    def test_actual_hydrogens_charges_and_exact_heavy_coordinates_for_supported_families(self):
        import numpy as np
        from cationize_protein.native_states import apply_native_states
        from cationize_protein.resources import default_parameters
        for source,target,delta in [('HIE','HIP',1),('ASP','ASH',1),('GLU','GLH',1),('LYS','LYN',-1),('CYS','CYM',-1)]:
            with self.subTest(source=source,target=target):
                original=native_peptide(source);before=original.coordinates.copy();q=sum(a.charge for a in original.atoms)
                typed,system,record=apply_native_states(original,{'2':target},default_parameters())
                self.assertEqual(typed.residues[1].name,target)
                self.assertAlmostEqual(sum(a.charge for a in typed.atoms)-q,delta,places=5)
                self.assertEqual(record['expected_charge_change_e'],delta)
                self.assertTrue(record['heavy_coordinates_exact'])
                np.testing.assert_array_equal(original.coordinates,before)
                self.assertEqual(system.getNumParticles(),len(typed.atoms))
                restored,_,reverse=apply_native_states(typed,{'2':source},default_parameters())
                self.assertAlmostEqual(sum(a.charge for a in restored.atoms)-q,0.,places=5)
                self.assertEqual(restored.residues[1].name,source)

    def test_other_native_change_preserves_added_ammonium_hydrogens_and_charge(self):
        import numpy as np
        import parmed
        from parmed.modeller import ResidueTemplate
        from cationize_protein.builder import replace_residue
        from cationize_protein.off import read_units
        from cationize_protein.parameter_assignment import parameterize
        from cationize_protein.native_states import apply_native_states
        from cationize_protein.resources import default_parameters
        source=native_peptide('HIE');u=read_units(default_parameters()/'GLX.lib')['GLX'];t=ResidueTemplate('GLX')
        for a,xyz in zip(u['atoms'],u['positions']):
            atom=parmed.Atom(name=a['name'],atomic_number=a['atomic_number'],mass=float(a['atomic_number']*2),charge=a['initial_template_charge_e'],type=a['amber_type'])
            atom.xx,atom.xy,atom.xz=xyz;t.add_atom(atom)
        for i,j in u['bonds']:t.add_bond(i,j)
        replaced=replace_residue(source,0,t)
        modified,_=parameterize(replaced,default_parameters())
        before={a.name:np.array([a.xx,a.xy,a.xz]) for a in modified.residues[0].atoms}
        result,_,record=apply_native_states(modified,{'2':'HIP'},default_parameters())
        self.assertEqual(record['introduced_ammonium_groups_retained'],1)
        for a in result.residues[0].atoms:np.testing.assert_array_equal([a.xx,a.xy,a.xz],before[a.name])
        self.assertEqual(sum(a.atomic_number==1 for a in result.residues[0].atoms),sum(a.atomic_number==1 for a in modified.residues[0].atoms))

    def test_missing_forcefield_state_is_explicitly_rejected_before_rebuild(self):
        from cationize_protein.native_states import apply_native_states
        from cationize_protein.resources import default_parameters
        with self.assertRaisesRegex(UnsupportedNativeState,'template'):
            apply_native_states(native_peptide('TYR'),{'2':'TYM'},default_parameters())
