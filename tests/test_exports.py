import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

RUNTIME = importlib.util.find_spec('openmm') is not None and importlib.util.find_spec('parmed') is not None


class ExportAPITests(unittest.TestCase):
    def test_unknown_format_is_rejected_without_runtime(self):
        from cationize_protein.exports import export_system
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, 'format'):
                export_system(None, d, formats=('unsupported',))


@unittest.skipUnless(RUNTIME, 'OpenMM and ParmEd runtime required')
class PhysicalRoundTripTests(unittest.TestCase):
    def fixture(self, cmap=False):
        import numpy as np
        import parmed as pmd
        import openmm as mm
        from openmm import unit
        structure = pmd.Structure()
        for i in range(5):
            structure.add_atom(pmd.Atom(name='C'+str(i), atomic_number=6, mass=12.), 'TST', 1)
        for i in range(4):
            structure.bonds.append(pmd.Bond(structure.atoms[i], structure.atoms[i+1]))
        xyz = np.array([[0.,0.,0.], [1.5,.1,0.], [2.5,1.2,.2], [3.9,1.1,.8], [4.6,2.3,1.1]])
        system = mm.System()
        for _ in structure.atoms: system.addParticle(12.)
        bonds = mm.HarmonicBondForce()
        for i in range(4): bonds.addBond(i,i+1,.15,200000.)
        system.addForce(bonds)
        angles = mm.HarmonicAngleForce()
        for i in range(3): angles.addAngle(i,i+1,i+2,1.9,400.)
        system.addForce(angles)
        torsions = mm.PeriodicTorsionForce()
        for i in range(2): torsions.addTorsion(i,i+1,i+2,i+3,3,.2,2.)
        system.addForce(torsions)
        nb = mm.NonbondedForce()
        for q in [.1,-.2,.2,-.2,.1]: nb.addParticle(q,.32,.4)
        nb.createExceptionsFromBonds([(i,i+1) for i in range(4)], 1/1.2,.5)
        system.addForce(nb)
        if cmap:
            force = mm.CMAPTorsionForce()
            idx = force.addMap(4, [float(i)/10 for i in range(16)])
            force.addTorsion(idx,0,1,2,3,1,2,3,4);system.addForce(force)
        typed = pmd.openmm.load_topology(structure.topology, system, xyz=xyz*unit.angstrom)
        typed._cationize_system = system
        return typed, system

    def test_split_itp_export_survives_relocation_and_preserves_cmap_and_pdb_identity(self):
        import shutil
        import numpy as np
        import parmed as pmd
        from cationize_protein.exports import export_system
        typed,system=self.fixture(cmap=True)
        typed.residues[0].name='ASC'
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'original'; output.mkdir()
            record=export_system(typed,output,formats=('amber','gromacs','openmm'),force_field='ff19SB')
            self.assertTrue((output/'protein.itp').is_file(),'The protein molecule must be reusable as a separate ITP')
            self.assertTrue((output/'forcefield.itp').is_file(),'Global types must be separate from molecule definitions')
            top=(output/'system.top').read_text()
            self.assertIn('#include "forcefield.itp"',top)
            self.assertIn('#include "protein.itp"',top)
            self.assertNotIn(str(output),top)
            moved=Path(tmp)/'relocated';shutil.move(str(output),moved)
            restored=pmd.load_file(str(moved/'system.top'),xyz=str(moved/'system.gro'))
            self.assertTrue(all(a.type.startswith('CP_') for a in restored.atoms))
            self.assertEqual(record['exports']['gromacs']['validation']['cmap_terms'],1)
            from cationize_protein.exports import _validate
            self.assertTrue(_validate(typed,system,restored)['passed'])
            pdb=pmd.load_file(str(moved/'protein.pdb'))
            self.assertEqual([(a.residue.name,a.name) for a in pdb.atoms],
                             [(a.residue.name,a.name) for a in typed.atoms])
            np.testing.assert_allclose(pdb.coordinates,typed.coordinates,atol=.00051,rtol=0)
            global_sections=(moved/'forcefield.itp').read_text()
            self.assertIn('[ atomtypes ]',global_sections)
            self.assertNotIn('[ moleculetype ]',global_sections)
            molecule=(moved/'protein.itp').read_text()
            self.assertIn('[ moleculetype ]',molecule)
            self.assertNotIn('[ defaults ]',molecule)

    def test_openmm_only_export_includes_a_pdb_without_requiring_amber_files(self):
        from cationize_protein.exports import export_system
        typed,_=self.fixture()
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)
            record=export_system(typed,output,formats=('openmm',))
            self.assertEqual(set(record['exports']),{'openmm'})
            self.assertTrue((output/'protein.pdb').is_file())
            self.assertTrue((output/'protein.lib').is_file())
            self.assertTrue((output/'parameters.frcmod').is_file())
            self.assertFalse((output/'system.parm7').exists())

    def test_protein_itp_can_be_combined_with_an_independent_ligand(self):
        import numpy as np
        import openmm as mm
        from openmm import unit
        import parmed as pmd
        from cationize_protein.exports import export_system,_validate
        typed,system=self.fixture(cmap=True)
        # Build the expected physical System independently of the includes:
        # five protein atoms plus a neutral three-atom ligand fragment.
        for _ in range(3):system.addParticle(12.)
        bonds=next(f for f in system.getForces() if type(f).__name__=='HarmonicBondForce')
        bonds.addBond(5,6,.15,200000.)
        bonds.addBond(6,7,.15,200000.)
        nb=next(f for f in system.getForces() if type(f).__name__=='NonbondedForce')
        for _ in range(3):nb.addParticle(0.,.32,.4)
        for i,j in ((5,6),(5,7),(6,7)):nb.addException(i,j,0.,.32,0.)
        graph=pmd.Structure()
        for a in typed.atoms:
            graph.add_atom(pmd.Atom(name=a.name,atomic_number=a.atomic_number,mass=a.mass),a.residue.name,1)
        for name in ('L1','L2','L3'):
            graph.add_atom(pmd.Atom(name=name,atomic_number=6,mass=12.),'LIG',2)
        for b in typed.bonds:graph.bonds.append(pmd.Bond(graph.atoms[b.atom1.idx],graph.atoms[b.atom2.idx]))
        graph.bonds.append(pmd.Bond(graph.atoms[5],graph.atoms[6]))
        graph.bonds.append(pmd.Bond(graph.atoms[6],graph.atoms[7]))
        xyz=np.vstack([typed.coordinates,[[20.,0.,0.],[21.5,0.,0.],[22.5,1.,0.]]])
        expected=pmd.openmm.load_topology(graph.topology,system,xyz=xyz*unit.angstrom)
        original=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
        # Recreate the unmodified fixture rather than mutating the source used
        # by the export under test.
        protein,_=self.fixture(cmap=True)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            export_system(protein,root,formats=('gromacs',),force_field='ff19SB')
            self.assertTrue((root/'protein.itp').is_file(),'Reusable molecule ITP is missing')
            (root/'ligand_types.itp').write_text('[ atomtypes ]\nLC 6 12.0 0.0 A 0.32 0.4\n')
            (root/'ligand.itp').write_text('[ moleculetype ]\nLIG 3\n[ atoms ]\n'
                '1 LC 1 LIG L1 1 0.0 12.0\n2 LC 1 LIG L2 1 0.0 12.0\n3 LC 1 LIG L3 1 0.0 12.0\n'
                '[ bonds ]\n1 2 1 0.15 200000.0\n2 3 1 0.15 200000.0\n')
            top=(root/'system.top').read_text().replace('#include "protein.itp"',
                '#include "ligand_types.itp"\n#include "protein.itp"\n#include "ligand.itp"')
            (root/'complex.top').write_text(top+'\nLIG 1\n')
            expected.save(str(root/'complex.gro'),format='gro')
            reloaded=pmd.load_file(str(root/'complex.top'),xyz=str(root/'complex.gro'))
            self.assertEqual(len(reloaded.atoms),8)
            self.assertEqual(len(reloaded.residues),2)
            from cationize_protein.exports import _nonbonded
            from openmm.app import NoCutoff
            actual=reloaded.createSystem(nonbondedMethod=NoCutoff,constraints=None,rigidWater=False,removeCMMotion=False)
            self.assertEqual(set(_nonbonded(original)[1]),set(_nonbonded(actual)[1]))
            self.assertTrue(_validate(expected,original,reloaded)['passed'])

    def test_amber_gromacs_and_openmm_preserve_physical_parameters(self):
        from cationize_protein.exports import export_system
        typed, _ = self.fixture()
        with tempfile.TemporaryDirectory() as d:
            record = export_system(typed, d, formats=('amber','gromacs','openmm'))
            self.assertEqual(record['status'], 'validated_physical_exports')
            self.assertEqual(set(record['exports']), {'amber','gromacs','openmm'})
            for name in ('amber','gromacs'):
                self.assertTrue(record['exports'][name]['validation']['passed'])
                self.assertLess(record['exports'][name]['validation']['max_charge_error_e'], 1e-6)
                self.assertLess(record['exports'][name]['validation']['max_sigma_error_nm'], 1e-6)
                self.assertLess(record['exports'][name]['validation']['max_epsilon_error_kJ_mol'], 1e-6)
            self.assertTrue((Path(d)/'system.xml').exists())
            for name in ('parameters.frcmod','protein.lib','PARAMETERS.json','rebuild.py','ATOM_TYPES.csv'):
                self.assertTrue((Path(d)/name).is_file(),name)
            parameters=json.loads((Path(d)/'PARAMETERS.json').read_text())
            self.assertEqual(len(parameters['atoms']),len(typed.atoms))
            self.assertEqual(len(parameters['bonds']),len(typed.bonds))
            self.assertIn('exceptions',parameters)
            self.assertEqual(json.loads((Path(d)/'EXPORTS.json').read_text())['status'], record['status'])

    def test_cmap_survives_amber_gromacs_and_openmm(self):
        from cationize_protein.exports import export_system
        typed, _ = self.fixture(cmap=True)
        with tempfile.TemporaryDirectory() as d:
            record = export_system(typed, d, formats=('amber','gromacs','openmm'), force_field='ff19SB')
            for name in ('amber','gromacs'):
                self.assertEqual(record['exports'][name]['validation']['cmap_terms'], 1)
            self.assertEqual(record['exports']['openmm']['validation']['cmap_terms'],1)
            parameters=json.loads((Path(d)/'PARAMETERS.json').read_text())
            self.assertEqual(len(parameters['cmap_terms']),1)
            self.assertEqual(len(parameters['cmap_maps'][0]['grid_kcal_mol']),16)

    def test_zero_epsilon_particle_sigma_is_physically_irrelevant_in_roundtrip(self):
        import parmed as pmd
        from openmm import unit
        from cationize_protein.exports import export_system
        typed,system=self.fixture()
        nb=next(f for f in system.getForces() if type(f).__name__=='NonbondedForce')
        q,_,_=nb.getParticleParameters(2);nb.setParticleParameters(2,q,1.,0.)
        # Rebuild exceptions to include the zero-LJ particle consistently.
        for i in range(nb.getNumExceptions()):
            a,b,q,s,e=nb.getExceptionParameters(i)
            if 2 in (a,b):nb.setExceptionParameters(i,a,b,q,s,0.)
        source=pmd.openmm.load_topology(typed.topology,system,xyz=typed.coordinates*unit.angstrom)
        source._cationize_system=system
        with tempfile.TemporaryDirectory() as d:
            result=export_system(source,d,formats=('amber','gromacs'))
            self.assertTrue(result['exports']['amber']['validation']['passed'])

    def test_complete_frcmod_declares_zero_torsions_and_leap_atom_hybridization(self):
        import parmed as pmd
        from openmm import unit
        from cationize_protein.exports import export_system
        typed,system=self.fixture()
        for i,f in enumerate(system.getForces()):
            if type(f).__name__=='PeriodicTorsionForce':system.removeForce(i);break
        source=pmd.openmm.load_topology(typed.topology,system,xyz=typed.coordinates*unit.angstrom)
        source._cationize_system=system
        with tempfile.TemporaryDirectory() as d:
            export_system(source,d)
            params=pmd.amber.AmberParameterSet(str(Path(d)/'parameters.frcmod'))
            self.assertGreaterEqual(len(params.dihedral_types),2)
            self.assertIn('addAtomTypes',(Path(d)/'rebuild.in').read_text())

    def test_custom_restraint_is_rejected_instead_of_silently_dropped(self):
        import openmm as mm
        from cationize_protein.exports import export_system
        typed, system = self.fixture()
        system.addForce(mm.CustomExternalForce('x*x'))
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, 'CustomExternalForce'):
                export_system(typed, d)
            self.assertFalse((Path(d)/'EXPORTS.json').exists())

    def test_implicit_solvent_force_is_rejected_before_export(self):
        import openmm as mm
        from cationize_protein.exports import export_system
        typed, system = self.fixture()
        system.addForce(mm.GBSAOBCForce())
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, 'GBSAOBCForce'):
                export_system(typed, d)

    def test_parameter_bundle_uses_actual_uniform_one_four_scaling(self):
        from openmm import unit
        from cationize_protein.exports import export_system
        typed, system = self.fixture()
        force = next(f for f in system.getForces() if type(f).__name__ == 'NonbondedForce')
        for i in range(force.getNumExceptions()):
            a,b,q,s,e = force.getExceptionParameters(i)
            if abs(q.value_in_unit(unit.elementary_charge**2)) > 0:
                force.setExceptionParameters(i,a,b,q*(.7/(1/1.2)),s,e*.5)
        with tempfile.TemporaryDirectory() as d:
            export_system(typed, d, formats=('amber',))
            physical=json.loads((Path(d)/'PARAMETERS.json').read_text())
            self.assertTrue(all(abs(t['SCEE']-1/.7)<1e-10 and abs(t['SCNB']-4)<1e-10 for t in physical['torsions']))


if __name__ == '__main__': unittest.main()
