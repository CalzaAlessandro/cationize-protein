"""Initial PROPKA chooses histidine charge before cationization."""
import importlib
from pathlib import Path
import tempfile
import unittest

from cationize_protein.preflight import check_pdb
from cationize_protein.preparation import prepare_pdb


class InitialHistidineTests(unittest.TestCase):
    def source(self, root, state='HIE', hydrogens=False):
        from test_native_states import native_peptide
        source=root/'input.pdb';structure=native_peptide('HIE')
        if state!='HIE':
            from cationize_protein.native_states import apply_native_states
            from cationize_protein.resources import default_parameters
            structure,_,_=apply_native_states(structure,{'2':state},default_parameters())
        structure.save(str(source),format='pdb')
        lines=[]
        for line in source.read_text().splitlines():
            if line.startswith(('ATOM','HETATM')):
                if not hydrogens and line[76:78].strip()=='H':continue
                if not hydrogens and line[17:20].strip() in {'HIE','HID','HIP'}:
                    line=line[:17]+'HIS'+line[20:]
            lines.append(line)
        source.write_text('\n'.join(lines)+'\n');return source

    def test_default_preparation_assigns_an_unlabelled_native_histidine(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=self.source(root);before=source.read_bytes()
            result=prepare_pdb(source,root/'prepared')
            report=check_pdb(root/'prepared/protein.pdb')
            self.assertTrue(report['ok'],report['messages'])
            self.assertIn(report['residues'][1]['residue_name'],{'HID','HIE','HIP'})
            self.assertEqual(len(result.get('histidines',[])),1,'Prepare must run PROPKA on native HIS')
            row=result['histidines'][0]
            self.assertEqual(row['state'],'HIP' if row['pKa']>7.4 else row['neutral_tautomer'])
            self.assertTrue((root/'prepared/histidines.csv').is_file())
            self.assertTrue((root/'prepared/protonation/PROPKA_COMPLETE.json').is_file())
            self.assertEqual(source.read_bytes(),before)
            before_heavy={(r['ordinal'],a['name']):a['coordinates_A'] for r in check_pdb(source)['residues'] for a in r['atoms'] if not a['hydrogen']}
            after_heavy={(r['ordinal'],a['name']):a['coordinates_A'] for r in report['residues'] for a in r['atoms'] if not a['hydrogen']}
            self.assertEqual(before_heavy,after_heavy)
            from cationize_protein.workflow import _pdb_structure
            from cationize_protein.resources import default_parameters
            typed,_=_pdb_structure(root/'prepared/protein.pdb',default_parameters())
            self.assertEqual(typed.residues[1].name,row['state'])

    def test_pH_changes_charge_and_rebuilds_ring_hydrogens(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=self.source(root,'HIP',hydrogens=True)
            low=prepare_pdb(source,root/'low',pH=1.)
            high=prepare_pdb(source,root/'high',pH=14.)
            self.assertEqual(low['histidines'][0]['state'],'HIP')
            self.assertIn(high['histidines'][0]['state'],{'HID','HIE'})
            for folder,count in [('low',2),('high',1)]:
                residue=check_pdb(root/folder/'protein.pdb')['residues'][1]
                self.assertEqual(len({a['name'] for a in residue['atoms']} & {'HD1','HE2'}),count)

    def test_fixed_option_preserves_supplied_labels_and_hydrogens(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=self.source(root,'HIP',hydrogens=True)
            prepare_pdb(source,root/'fixed',protonation='fixed')
            self.assertEqual(source.read_bytes(),(root/'fixed/protein.pdb').read_bytes())
            self.assertFalse((root/'fixed/protonation').exists())

    def test_invalid_pH_is_rejected_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=self.source(root)
            with self.assertRaisesRegex(ValueError,'pH'):
                prepare_pdb(source,root/'bad',pH=float('nan'))
            self.assertFalse((root/'bad').exists())

    def test_cli_defaults_to_propka_at_ph74(self):
        from cationize_protein.cli import parser
        args=parser().parse_args(['prepare','input.pdb','--output','prepared'])
        self.assertEqual(getattr(args,'protonation',None),'propka')
        self.assertEqual(getattr(args,'ph',None),7.4)
