"""The installed interface must expose only supported engines and formats."""
import contextlib
import importlib
import io
from pathlib import Path
import tempfile
import tomllib
import unittest


class PublicInterfaceTests(unittest.TestCase):
    def cli(self):
        try:
            return importlib.import_module('cationize_protein.cli')
        except ImportError:
            self.fail('The public Python package must be named cationize_protein')

    def test_distribution_entrypoint_and_conda_use_the_public_name(self):
        root=Path(__file__).resolve().parents[1]
        metadata=tomllib.loads((root/'pyproject.toml').read_text())
        self.assertEqual(metadata['project']['name'],'cationize-protein')
        self.assertEqual(metadata['project']['scripts'],{'cationize-protein':'cationize_protein.cli:main'})
        self.assertEqual((root/'environment.yml').read_text().splitlines()[0],'name: cationize-protein')

    def test_cli_accepts_only_the_three_supported_formats(self):
        cli=self.cli()
        command=cli.parser()
        self.assertEqual(command.prog,'cationize-protein')
        args=command.parse_args(['run','--pdb','native.pdb','--output','out',
                                 '--formats','amber','openmm','gromacs'])
        self.assertEqual(args.formats,['amber','openmm','gromacs'])
        for unsupported in ('namd','charmm'):
            with self.subTest(format=unsupported), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    command.parse_args(['run','--pdb','native.pdb','--output','out','--formats',unsupported])

    def test_removed_export_formats_fail_before_loading_a_structure(self):
        self.cli()
        from cationize_protein.exports import export_system
        with tempfile.TemporaryDirectory() as tmp:
            for unsupported in ('namd','charmm'):
                with self.subTest(format=unsupported), self.assertRaisesRegex(ValueError,'format'):
                    export_system(None,Path(tmp)/unsupported,formats=(unsupported,))
            self.assertEqual(list(Path(tmp).iterdir()),[])


if __name__=='__main__':
    unittest.main()
