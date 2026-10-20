"""Initial preparation repairs atoms without changing supplied chemistry."""
import contextlib
import hashlib
import importlib
import io
from pathlib import Path
import tempfile
import unittest

from cationize_protein.preflight import check_pdb

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT/'examples/pdb_checks/reference_peptide.pdb'


class PreparationTests(unittest.TestCase):
    def prepare(self, *args):
        try:
            module = importlib.import_module('cationize_protein.preparation')
        except ModuleNotFoundError:
            self.fail('Initial preparation must provide the PDBFixer prepare_pdb function')
        return module.prepare_pdb(*args)

    def test_missing_sidechain_and_terminal_atom_are_repaired_with_source_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'input.pdb'; output=root/'prepared'
            original='\n'.join(line for line in REFERENCE.read_text().splitlines()
                if not (line.startswith('ATOM') and line[12:16].strip() in {'OD2','OXT'}))+'\n'
            source.write_text(original)
            result=self.prepare(source, output)
            self.assertEqual(source.read_text(),original)
            self.assertEqual((output/'original.pdb').read_bytes(),source.read_bytes())
            repaired=check_pdb(output/'protein.pdb')
            self.assertTrue(repaired['ok'], repaired['messages'])
            self.assertEqual(repaired['residue_count'],3)
            self.assertEqual(len(result['added_atoms']),3)
            self.assertEqual({a['atom_name'] for a in result['added_atoms']},{'OD2','OXT'})
            original_lines=[line for line in original.splitlines() if line.startswith('ATOM')]
            repaired_lines=(output/'protein.pdb').read_text().splitlines()
            for line in original_lines:
                self.assertIn(line,repaired_lines)
            self.assertIn('OD2',(output/'preparation.log').read_text())
            self.assertIn('PROPKA',(output/'preparation.log').read_text())

    def test_missing_backbone_atom_can_be_repaired_before_parameter_assignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'prepared'
            self.prepare(ROOT/'examples/pdb_checks/missing_backbone_atom.pdb',output)
            report=check_pdb(output/'protein.pdb')
            self.assertTrue(report['ok'],report['messages'])
            self.assertEqual(report['atom_count'],46)
            from cationize_protein.workflow import _pdb_structure
            from cationize_protein.resources import default_parameters
            structure,system=_pdb_structure(output/'protein.pdb',default_parameters())
            self.assertEqual(len(structure.residues),3)
            self.assertEqual(system.getNumParticles(),len(structure.atoms))

    def test_complete_structure_is_not_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'prepared';result=self.prepare(REFERENCE,output)
            self.assertEqual(result['added_atoms'],[])
            self.assertEqual((output/'protein.pdb').read_bytes(),REFERENCE.read_bytes())

    def test_invalid_columns_fail_without_creating_an_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'prepared'
            with self.assertRaisesRegex(ValueError,'truncated|fixed-column'):
                self.prepare(ROOT/'examples/pdb_checks/wrong_columns.pdb',output)
            self.assertFalse(output.exists())

    def test_missing_sequence_is_reported_without_building_residues(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'gap.pdb';output=root/'prepared'
            chain_a='\n'.join(line[:21]+'A'+line[22:] if line.startswith('ATOM') else line
                              for line in REFERENCE.read_text().splitlines())+'\n'
            source.write_text('SEQRES   1 A    4  ASP ASP LEU GLY\n'+chain_a)
            with self.assertRaisesRegex(ValueError,'Missing residues|missing residues'):
                self.prepare(source,output)
            self.assertFalse(output.exists())

    def test_custom_cationized_residues_are_never_replaced(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'modified.pdb';output=root/'prepared'
            source.write_text('\n'.join(line[:17]+'ASC'+line[20:]
                if line.startswith('ATOM') and line[22:26].strip()=='2' else line
                for line in REFERENCE.read_text().splitlines())+'\n')
            before=source.read_bytes()
            with self.assertRaisesRegex(ValueError,'native|cationized'):
                self.prepare(source,output)
            self.assertEqual(source.read_bytes(),before)
            self.assertFalse(output.exists())

    def test_alternates_are_not_silently_selected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'alternates.pdb';output=root/'prepared'
            lines=REFERENCE.read_text().splitlines()
            index=next(i for i,line in enumerate(lines) if line.startswith('ATOM'))
            lines[index]=lines[index][:16]+'B'+lines[index][17:]
            source.write_text('\n'.join(lines)+'\n')
            with self.assertRaisesRegex(ValueError,'alternate'):
                self.prepare(source,output)
            self.assertFalse(output.exists())

    def test_declared_states_and_disulfide_records_are_preserved(self):
        from test_preflight import atom
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'states.pdb';output=root/'prepared'
            # Separate chains avoid creating an artificial peptide bond.
            rows=[]
            for number,chain,offset in [(1,'A',0.),(2,'B',2.05)]:
                coords={'N':(-2.,3.,0.),'CA':(-2.,1.5,0.),'C':(-3.,.5,0.),
                        'O':(-4.,.5,0.),'OXT':(-2.8,-.7,0.),'CB':(-.7,1.,0.),'SG':(0.,0.,0.)}
                for name,xyz in coords.items():
                    rows.append(atom(len(rows)+1,name,'CYX',number,chain,
                        (xyz[0]+offset,xyz[1],xyz[2]),element='S' if name=='SG' else name[0]))
                rows.append('TER')
            decl=list(' '*80);decl[:6]='SSBOND';decl[15]='A';decl[17:21]='   1';decl[29]='B';decl[31:35]='   2'
            source.write_text(''.join(decl)+'\n'+'\n'.join(rows)+'\nEND\n')
            result=self.prepare(source,output)
            prepared=check_pdb(output/'protein.pdb')
            self.assertEqual([r['residue_name'] for r in prepared['residues']],['CYX','CYX'])
            self.assertEqual(len(prepared['ssbond_declarations']),1)
            self.assertEqual(result['added_atoms'],[])
            self.assertEqual((output/'protein.pdb').read_bytes(),source.read_bytes())

    def test_existing_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'prepared';output.mkdir();sentinel=output/'keep';sentinel.write_text('keep')
            with self.assertRaisesRegex(ValueError,'exist|new directory'):
                self.prepare(REFERENCE,output)
            self.assertEqual(sentinel.read_text(),'keep')
            self.assertEqual(sorted(p.name for p in output.iterdir()),['keep'])

    def test_prepare_command_reports_added_atoms_on_screen(self):
        from cationize_protein.cli import main,parser
        self.assertIn('prepare',parser()._subparsers._group_actions[0].choices,
                      'Expose initial PDBFixer preparation in the CLI')
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()) as stream:
            code=main(['prepare',str(ROOT/'examples/pdb_checks/missing_backbone_atom.pdb'),
                       '--output',str(Path(tmp)/'prepared')])
            self.assertEqual(code,0)
            self.assertIn('CA',stream.getvalue())
            self.assertIn('protein.pdb',stream.getvalue())

    def test_pdbfixer_is_installed_and_reported_as_a_runtime_dependency(self):
        from cationize_protein.dependencies import PYTHON_PACKAGES
        self.assertIn('pdbfixer',PYTHON_PACKAGES)
        self.assertIn('pdbfixer>=1.12', (ROOT/'requirements.txt').read_text())


if __name__=='__main__':
    unittest.main()
