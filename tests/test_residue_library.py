"""Reusable residue templates are distinct from the complete protein unit."""
import importlib
from pathlib import Path
import tempfile
import unittest


class ResidueLibraryTests(unittest.TestCase):
    def export(self,*args):
        try:module=importlib.import_module('cationize_protein.residue_library')
        except ModuleNotFoundError:self.fail('Run outputs must include a reusable residue_templates library')
        return module.export_residue_library(*args)

    def test_copy_contains_all_32_separate_off_libraries_and_verifies_resume(self):
        from cationize_protein.resources import default_parameters
        source=default_parameters()
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'residue_templates'
            result=self.export(source,target)
            self.assertEqual(result['template_count'],32)
            self.assertEqual(len(list(target.glob('*.lib'))),32)
            for lib in source.glob('*.lib'):
                self.assertEqual(lib.read_bytes(),(target/lib.name).read_bytes())
            self.assertTrue((target/'README.md').is_file())
            self.assertEqual(result,self.export(source,target))
            (target/'ASC.lib').write_text('corrupted')
            with self.assertRaisesRegex(ValueError,'hash|changed|differs'):
                self.export(source,target)
            self.assertEqual((target/'ASC.lib').read_text(),'corrupted')

    def test_missing_or_modified_source_template_is_rejected_without_output(self):
        import hashlib,json
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'parameters';source.mkdir()
            (source/'manifest.json').write_text(json.dumps({'templates':[{'name':'ASC','status':'validated','library_sha256':hashlib.sha256(b'expected').hexdigest()}]}))
            for supplied in (None,b'bad'):
                if supplied is not None:(source/'ASC.lib').write_bytes(supplied)
                with self.assertRaisesRegex(ValueError,'missing|hash|differs'):
                    self.export(source,root/'output')
                self.assertFalse((root/'output').exists())
