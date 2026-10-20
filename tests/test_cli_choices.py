"""Interactive scientific choices must reach the run unchanged and resume safely."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cationize_protein import cli


class RuntimeChoiceTests(unittest.TestCase):
    def arguments(self, output, *options):
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                return cli.parser().parse_args([
                    'run', '--pdb', 'native.pdb', '--output', str(output), *options])
        except SystemExit:
            self.fail('A run must accept the default probe without physics flags')

    def availability(self):
        return {'engine': {'openmm': True, 'amber': True, 'sander': True, 'pmemd': False},
                'solver': {'apbs': True, 'pbsa': False, 'coulomb': True},
                'force_field': {'ff14SB': True, 'ff19SB': True, 'ff99SBildn': True}}

    def test_probe_has_a_documented_default_and_can_be_overridden(self):
        self.assertEqual(self.arguments('out').probe_radius, 1.369)
        self.assertEqual(self.arguments('out', '--probe-radius', '1.5').probe_radius, 1.5)

    def test_terminal_selections_override_installation_defaults_without_writing_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = self.arguments(root/'out')
            with patch.object(cli.sys.stdin, 'isatty', return_value=True), \
                 patch.object(cli, '_available_choices', return_value=self.availability()), \
                 patch('builtins.input', side_effect=['sander', 'coulomb', 'ff19SB']), \
                 contextlib.redirect_stdout(io.StringIO()):
                selected = cli.run_choices(args, {}, root)
            self.assertEqual(selected, {'engine':'sander', 'solver':'coulomb', 'force_field':'ff19SB'})
            self.assertEqual(list(root.iterdir()), [])

    def test_explicit_flags_never_prompt_even_in_a_terminal(self):
        args = self.arguments('out', '--engine', 'openmm', '--solver', 'coulomb',
                              '--force-field', 'ff19SB')
        with patch.object(cli.sys.stdin, 'isatty', return_value=True), \
             patch('builtins.input', side_effect=AssertionError('Unexpected prompt')):
            selected = cli.run_choices(args, {}, Path('.'))
        self.assertEqual(selected, {'engine':'openmm', 'solver':'coulomb', 'force_field':'ff19SB'})

    def test_noninteractive_execution_uses_saved_choices(self):
        args = self.arguments('out')
        stored = {'engine':'sander', 'solver':'pbsa', 'force_field':'ff99SBildn'}
        with patch.object(cli.sys.stdin, 'isatty', return_value=False), \
             patch('builtins.input', side_effect=AssertionError('Unexpected prompt')):
            selected = cli.run_choices(args, stored, Path('.'))
        self.assertEqual(selected, stored)

    def test_resume_recovers_original_choices_instead_of_changed_installation_defaults(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); output=root/'out'; output.mkdir()
            original={'schema_version':2, 'engine':'sander', 'solver':'pbsa', 'force_field':'ff19SB'}
            content=json.dumps(original); (output/'RUN.json').write_text(content)
            args=self.arguments(output)
            with patch.object(cli.sys.stdin, 'isatty', return_value=True), \
                 patch('builtins.input', side_effect=AssertionError('A resume must not prompt')):
                selected=cli.run_choices(args, {'engine':'openmm', 'solver':'apbs', 'force_field':'ff14SB'}, root)
            self.assertEqual(selected, {'engine':'sander', 'solver':'pbsa', 'force_field':'ff19SB'})
            self.assertEqual((output/'RUN.json').read_text(), content)

    def test_missing_solver_cannot_be_selected_or_silently_replaced(self):
        args=self.arguments('out', '--engine','openmm','--force-field','ff14SB')
        available=self.availability(); available['solver']['apbs']=False
        with patch.object(cli.sys.stdin,'isatty',return_value=True), \
             patch.object(cli,'_available_choices',return_value=available), \
             patch('builtins.input',return_value='apbs'), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, 'apbs.*unavailable'):
                cli.run_choices(args, {}, Path('.'))

    def test_numbered_menu_accepts_enter_for_defaults(self):
        args=self.arguments('out')
        with patch.object(cli.sys.stdin,'isatty',return_value=True), \
             patch.object(cli,'_available_choices',return_value=self.availability()), \
             patch('builtins.input',side_effect=['', '3', '']), contextlib.redirect_stdout(io.StringIO()):
            selected=cli.run_choices(args, {}, Path('.'))
        self.assertEqual(selected, {'engine':'openmm','solver':'coulomb','force_field':'ff14SB'})

    def test_cancelled_selection_does_not_create_an_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); args=self.arguments(root/'out')
            with patch.object(cli.sys.stdin,'isatty',return_value=True), \
                 patch.object(cli,'_available_choices',return_value=self.availability()), \
                 patch('builtins.input',side_effect=EOFError), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'cancelled'):
                    cli.run_choices(args, {}, root)
            self.assertEqual(list(root.iterdir()), [])


if __name__=='__main__':
    unittest.main()
