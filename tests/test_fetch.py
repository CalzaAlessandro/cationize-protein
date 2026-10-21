"""Fetch preserves downloaded records and does not guess a structural model."""
import importlib
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]


class FetchTests(unittest.TestCase):
    def module(self):
        try:return importlib.import_module('cationize_protein.fetch')
        except ModuleNotFoundError:self.fail('The package must provide PDB and UniProt retrieval')

    def test_pdb_download_is_preserved_and_hashed(self):
        module=self.module();data=(ROOT/'examples/pdb_checks/reference_peptide.pdb').read_bytes()
        with tempfile.TemporaryDirectory() as tmp,patch.object(module,'download',return_value=data):
            dest=Path(tmp)/'fetch';record=module.fetch(dest,pdb_id='1n5u')
            self.assertEqual((dest/'1N5U.pdb').read_bytes(),data)
            self.assertEqual(record['pdb_id'],'1N5U')
            self.assertIn('1N5U.pdb',record['files_sha256'])
            self.assertTrue((dest/'fetch.log').is_file())

    def test_uniprot_keeps_full_sequence_and_separates_mature_chain(self):
        module=self.module()
        raw=json.dumps(dict(primaryAccession='P02768',sequence={'value':'MMMMDAHK'},
            features=[{'type':'Signal','location':{'start':{'value':1},'end':{'value':4}}},
                      {'type':'Chain','description':'Albumin','location':{'start':{'value':5},'end':{'value':8}}}],
            uniProtKBCrossReferences=[{'database':'PDB','id':'1N5U','properties':[{'key':'Method','value':'X-ray'},{'key':'Chains','value':'A=5-8'}]}])).encode()
        with tempfile.TemporaryDirectory() as tmp,patch.object(module,'download',return_value=raw):
            dest=Path(tmp)/'fetch';record=module.fetch(dest,uniprot_id='P02768')
            self.assertIn('MMMMDAHK',(dest/'sequence.fasta').read_text())
            self.assertIn('\nDAHK\n',(dest/'mature.fasta').read_text())
            self.assertIn('1N5U',(dest/'pdb_candidates.csv').read_text())
            self.assertEqual(record['selected_structure'],None)
            self.assertFalse(list(dest.glob('*.pdb')))

    def test_bad_ids_and_existing_output_are_rejected_before_network(self):
        module=self.module()
        with tempfile.TemporaryDirectory() as tmp,patch.object(module,'download',side_effect=AssertionError('No network')):
            dest=Path(tmp)/'new'
            for kwargs in ({'pdb_id':'../bad'},{'uniprot_id':'x/y'},{},):
                with self.assertRaises(ValueError):module.fetch(dest,**kwargs)
                self.assertFalse(dest.exists())
            dest.mkdir();(dest/'keep').write_text('keep')
            with self.assertRaisesRegex(ValueError,'exist|new directory'):module.fetch(dest,pdb_id='1N5U')
            self.assertEqual((dest/'keep').read_text(),'keep')

    def test_html_error_page_is_not_saved_as_a_structure(self):
        module=self.module()
        with tempfile.TemporaryDirectory() as tmp,patch.object(module,'download',return_value=b'<html>Error</html>'):
            dest=Path(tmp)/'bad'
            with self.assertRaisesRegex(ValueError,'PDB|coordinate'):module.fetch(dest,pdb_id='1N5U')
            self.assertFalse(dest.exists())

    def test_multiple_mature_chains_are_not_silently_selected(self):
        module=self.module()
        raw=json.dumps({'primaryAccession':'P02768','sequence':{'value':'ACDE'},'features':[
            {'type':'Chain','location':{'start':{'value':1},'end':{'value':2}}},
            {'type':'Chain','location':{'start':{'value':3},'end':{'value':4}}}]}).encode()
        with tempfile.TemporaryDirectory() as tmp,patch.object(module,'download',return_value=raw):
            dest=Path(tmp)/'fetch';module.fetch(dest,uniprot_id='P02768')
            self.assertFalse((dest/'mature.fasta').exists())
            self.assertIn('ACDE',(dest/'sequence.fasta').read_text())

    def test_cli_accepts_pdb_and_uniprot_together(self):
        from cationize_protein.cli import parser
        choices=parser()._subparsers._group_actions[0].choices
        self.assertIn('fetch',choices)
        args=parser().parse_args(['fetch','--pdb','1N5U','--uniprot','P02768','--output','native'])
        self.assertEqual(args.pdb,'1N5U');self.assertEqual(args.uniprot,'P02768')

    def test_uniprot_disulfide_table_retains_canonical_positions_and_evidence(self):
        module=self.module()
        raw=json.dumps({'primaryAccession':'TEST01','sequence':{'value':'MCCK'},'features':[
            {'type':'Disulfide bond','location':{'start':{'value':2,'modifier':'EXACT'},
                                               'end':{'value':3,'modifier':'EXACT'}},
             'evidences':[{'evidenceCode':'ECO:0000269','source':'PubMed','id':'123'}]}]}).encode()
        with tempfile.TemporaryDirectory() as tmp,patch.object(module,'download',return_value=raw):
            dest=Path(tmp)/'fetch';result=module.fetch(dest,uniprot_id='TEST01')
            self.assertTrue((dest/'disulfides.csv').is_file(), 'Export the UniProt disulfide annotations')
            with (dest/'disulfides.csv').open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual([(r['uniprot_start'],r['uniprot_end']) for r in rows],[('2','3')])
            self.assertIn('PubMed',rows[0]['evidence'])
            self.assertIn('disulfides.csv',result['files_sha256'])
