"""The bootstrap must preserve unrelated environments and distinguish capabilities."""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from cationize_protein import dependencies as deps


class DependencyTests(unittest.TestCase):
    def test_install_completes_without_external_programs_or_scientific_prompts(self):
        # The process boundary stands in for pip/venv; discovery and settings
        # serialization are real. Missing optional executables must not undo a
        # successfully installed Python package or trigger external downloads.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); environment=root/'.venv'; environment.mkdir()
            python=Path(self.supported_python())
            probe={'executable':str(python),'usable':True,'packages':{
                name:{'satisfies':True,'version':spec[1],'minimum':spec[1]}
                for name,spec in deps.PYTHON_PACKAGES.items()}}
            with patch.object(deps,'ensure_environment',return_value=python), \
                 patch.object(deps,'install_python_package',return_value=None), \
                 patch.object(deps,'_python_probe',return_value=probe), \
                 patch.object(deps,'_force_field_report',return_value={
                     name:{'supported':True} for name in deps.FORCE_FIELDS}), \
                 patch.object(deps,'find_executable',return_value=None), \
                 patch.object(deps,'install_apbs',side_effect=AssertionError('No external download requested')), \
                 patch.object(deps.sys.stdin,'isatty',return_value=True), \
                 patch('builtins.input',side_effect=AssertionError('Choose scientific settings at runtime')), \
                 contextlib.redirect_stdout(io.StringIO()) as output:
                code=deps.installer_main([], root=root)
            self.assertEqual(code, 0)
            self.assertIn('Python package installation complete', output.getvalue())
            self.assertIn('cationize-protein --help',output.getvalue())
            self.assertIn('Activate',output.getvalue())
            self.assertIn('GROMACS', output.getvalue())
            self.assertIn('not found', output.getvalue())
            self.assertEqual(json.loads((environment/'cationize-protein-settings.json').read_text())['engine'],'openmm')
            self.assertEqual(json.loads((environment/'cationize-protein-settings.json').read_text())['solver'],'coulomb')
            self.assertIn('Selected workflow dependencies are available.',output.getvalue())
            self.assertFalse((root/'.tools').exists())

    def test_gromacs_is_detected_again_after_a_later_installation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); environment={'PATH':str(root)}
            first=deps.inspect_dependencies(root, python=sys.executable, environ=environment)
            binary=root/'gmx'; binary.write_text('#!/bin/sh\nexit 0\n'); binary.chmod(0o755)
            second=deps.inspect_dependencies(root, python=sys.executable, environ=environment)
            self.assertFalse(first['capabilities']['gromacs']['available'])
            self.assertTrue(second['capabilities']['gromacs']['available'])
            self.assertEqual(second['executables']['gmx'], str(binary))

    def supported_python(self):
        try:
            return deps._select_python(None)
        except deps.DependencyError:
            self.skipTest('No supported Python 3.10-3.14 available for a real venv test')

    def test_python_discovery_prefers_supported_named_interpreters(self):
        versions = {'/available/python3.13': [3, 13, 1], '/available/python3.11': [3, 11, 9],
                    sys.executable: [3, 15, 0]}
        names = {'python3.13': '/available/python3.13', 'python3.11': '/available/python3.11'}
        def result(command, **kwargs):
            return subprocess.CompletedProcess(command, 0, json.dumps(versions[command[0]]), '')
        with patch.object(deps.shutil, 'which', side_effect=lambda name: names.get(name)), \
             patch.object(deps.subprocess, 'run', side_effect=result):
            self.assertEqual(deps._select_python(None), '/available/python3.13')

    def test_explicit_unsupported_python_fails_before_creating_environment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad = '/unsupported/python3.15'
            result = subprocess.CompletedProcess([], 0, '[3,15,0]', '')
            output = io.StringIO()
            with patch.object(deps.subprocess, 'run', return_value=result), contextlib.redirect_stdout(output):
                code = deps.installer_main(['--python', bad, '--json'], root=root)
            self.assertEqual(code, 1)
            self.assertIn('3.10', json.loads(output.getvalue())['error'])
            self.assertEqual(list(root.iterdir()), [])

    def test_relative_path_executables_remain_usable_after_working_directory_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / 'apbs'
            executable.write_text('#!/bin/sh\nexit 0\n')
            executable.chmod(0o755)
            original = Path.cwd()
            try:
                os.chdir(root)
                found = deps.find_executable('apbs', root, {'PATH': '.'})
            finally:
                os.chdir(original)
            self.assertTrue(Path(found).is_absolute())
            result = subprocess.run([found], cwd=original, capture_output=True)
            self.assertEqual(result.returncode, 0)

    def test_missing_amber_never_becomes_available_from_openmm(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = deps.inspect_dependencies(Path(tmp), python=sys.executable, environ={'PATH': ''})
            self.assertFalse(report['capabilities']['amber_minimization']['available'])
            self.assertFalse(report['capabilities']['pbsa']['available'])
            self.assertFalse(report['capabilities']['apbs']['available'])
            self.assertEqual(report['executables']['tleap'], None)

    def test_cuda_only_amber_is_not_a_serial_minimizer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);binary=root/'amber/bin/pmemd.cuda';binary.parent.mkdir(parents=True)
            binary.write_text('#!/bin/sh\nexit 0\n');binary.chmod(0o755)
            report=deps.inspect_dependencies(root,python=sys.executable,environ={'AMBERHOME':str(root/'amber'),'PATH':''},engine='amber',solver='coulomb')
            self.assertIsNotNone(report['executables']['pmemd.cuda'])
            self.assertFalse(report['selected']['engine_available'])

    def test_amberhome_and_private_apbs_are_detected_separately(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            amber = root / 'amber'
            for name in ['sander', 'tleap', 'pbsa']:
                binary = amber / 'bin' / name
                binary.parent.mkdir(parents=True, exist_ok=True)
                binary.write_text('#!/bin/sh\nexit 0\n')
                binary.chmod(0o755)
            apbs = root / '.tools/apbs/bin/apbs'
            apbs.parent.mkdir(parents=True)
            apbs.write_text('#!/bin/sh\nexit 0\n')
            apbs.chmod(0o755)
            report = deps.inspect_dependencies(root, python=sys.executable,
                                              environ={'AMBERHOME': str(amber), 'PATH': ''})
            self.assertTrue(report['capabilities']['amber_minimization']['available'])
            self.assertTrue(report['capabilities']['pbsa']['available'])
            self.assertEqual(report['executables']['apbs'], str(apbs))
            self.assertIsNone(report['executables']['pmemd'])

    def test_forcefield_selection_does_not_claim_unvalidated_compatibility(self):
        self.assertEqual(deps.validate_force_field('ff14SB'), 'ff14SB')
        self.assertEqual(deps.validate_force_field('ff19SB'), 'ff19SB')
        self.assertEqual(deps.validate_force_field('ff99SBildn'), 'ff99SBildn')
        with self.assertRaisesRegex(deps.DependencyError, 'validated'):
            deps.validate_force_field('charmm36')

    def test_unknown_directory_and_project_root_are_never_used_as_venv(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            existing = root / 'important'
            existing.mkdir()
            sentinel = existing / 'keep.txt'
            sentinel.write_text('original bytes')
            for path in [root, root.parent, existing, root / 'src']:
                with self.subTest(path=path), self.assertRaises(deps.DependencyError):
                    deps.ensure_environment(root, path, sys.executable)
            self.assertEqual(sentinel.read_text(), 'original bytes')
            self.assertEqual(sorted(p.name for p in existing.iterdir()), ['keep.txt'])

    def test_symlink_target_is_never_followed_for_install(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            other = root / 'other'
            other.mkdir()
            alias = root / '.venv'
            alias.symlink_to(other, target_is_directory=True)
            with self.assertRaises(deps.DependencyError):
                deps.ensure_environment(root, alias, sys.executable)
            self.assertEqual(list(other.iterdir()), [])

    def test_real_new_environment_is_private_and_can_be_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # No network or pip is needed to exercise real venv creation here.
            with patch.object(deps, '_create_venv', side_effect=lambda py, dst: subprocess.run(
                    [py, '-m', 'venv', '--without-pip', str(dst)], check=True, capture_output=True)):
                python = deps.ensure_environment(root, root / '.venv', self.supported_python())
            result = subprocess.run([str(python), '-c', 'import sys; print(sys.prefix != sys.base_prefix)'],
                                    capture_output=True, text=True, check=True)
            self.assertEqual(result.stdout.strip(), 'True')
            self.assertEqual(deps.ensure_environment(root, root / '.venv', self.supported_python()), python)

    def test_failed_pip_is_reported_without_leaking_credentials(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def failure(*args, **kwargs):
                return subprocess.CompletedProcess(args[0], 1,
                    stdout='https://user:secret@example.invalid/index', stderr='token=secret')
            with patch.object(deps.subprocess, 'run', side_effect=failure):
                with self.assertRaises(deps.DependencyError) as caught:
                    deps.install_python_package(root, Path('/private/env/bin/python'))
            self.assertIn('pip', str(caught.exception))
            self.assertNotIn('secret', str(caught.exception))
            self.assertNotIn('user:', str(caught.exception))

    def test_failed_venv_creation_is_not_reported_as_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with patch.object(deps.subprocess, 'run', return_value=subprocess.CompletedProcess([], 1, '', 'failure')):
                with self.assertRaisesRegex(deps.DependencyError, 'create.*environment'):
                    deps.ensure_environment(root, root / '.venv', sys.executable)
            self.assertTrue((root / '.venv' / deps.ENV_MARKER).is_file())
            self.assertFalse((root / '.venv/bin/python').exists())

    def test_apbs_library_paths_are_set_only_for_the_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / 'apbs/bin/apbs'
            executable.parent.mkdir(parents=True)
            (root / 'apbs/lib').mkdir()
            original = {'PATH': '/existing/bin', 'LD_LIBRARY_PATH': '/existing/lib'}
            child = deps.external_program_environment(executable, original)
            self.assertEqual(original, {'PATH': '/existing/bin', 'LD_LIBRARY_PATH': '/existing/lib'})
            self.assertEqual(child['LD_LIBRARY_PATH'].split(os.pathsep)[0], str(root / 'apbs/lib'))

    def test_apbs_bad_digest_never_installs_an_executable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'config').mkdir()
            asset = {'url': 'https://github.com/Electrostatics/apbs/releases/download/v3.4.1/test.zip',
                     'sha256': '0' * 64, 'bundle': 'bundle'}
            (root / 'config/install-manifest.json').write_text(json.dumps({'apbs': {'assets': {'Linux-x86_64': asset}}}))
            incoming = io.BytesIO(b'wrong bytes')
            incoming.geturl = lambda: asset['url']
            with patch.object(deps.platform, 'system', return_value='Linux'), \
                 patch.object(deps.platform, 'machine', return_value='x86_64'), \
                 patch.object(deps, 'urlopen', return_value=incoming):
                with self.assertRaisesRegex(deps.DependencyError, 'SHA256'):
                    deps.install_apbs(root)
            self.assertFalse((root / '.tools/apbs').exists())
            self.assertEqual(list((root / '.tools').iterdir()), [])

    def test_unknown_apbs_architecture_is_reported_without_fallback_or_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'config').mkdir()
            (root / 'config/install-manifest.json').write_text('{"apbs": {"assets": {}}}')
            with patch.object(deps.platform, 'system', return_value='Darwin'), \
                 patch.object(deps.platform, 'machine', return_value='arm64'):
                with self.assertRaisesRegex(deps.DependencyError, 'No pinned official APBS binary.*arm64'):
                    deps.install_apbs(root)
            self.assertFalse((root / '.tools').exists())

    def test_check_uses_saved_choices_without_modifying_settings(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            environment = root / '.venv'
            python = environment / 'bin/python'
            python.parent.mkdir(parents=True)
            python.symlink_to(sys.executable)
            settings = environment / 'cationize-protein-settings.json'
            content = '{"schema_version":1,"engine":"sander","solver":"pbsa","force_field":"ff19SB"}'
            settings.write_text(content)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                code = deps.installer_main(['--check', '--json'], root=root)
            report = json.loads(output.getvalue())
            self.assertEqual(report['selected']['engine'], 'sander')
            self.assertEqual(report['selected']['solver'], 'pbsa')
            self.assertEqual(report['selected']['force_field'], 'ff19SB')
            self.assertEqual(settings.read_text(), content)
            self.assertEqual(code, 1)

    def test_settings_symlink_is_rejected_without_overwriting_its_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            environment = root / '.venv'
            with patch.object(deps, '_create_venv', side_effect=lambda py, dst: subprocess.run(
                    [py, '-m', 'venv', '--without-pip', str(dst)], check=True, capture_output=True)):
                deps.ensure_environment(root, environment, self.supported_python())
            other = root / 'unrelated.json'
            other.write_text('{"preserve":"original"}')
            (environment / 'cationize-protein-settings.json').symlink_to(other)
            output = io.StringIO()
            with contextlib.redirect_stdout(output), patch.object(deps, 'install_python_package', return_value=None):
                code = deps.installer_main(['--json', '--solver', 'coulomb'], root=root)
            self.assertEqual(code, 1)
            self.assertIn('error', json.loads(output.getvalue()))
            self.assertEqual(other.read_text(), '{"preserve":"original"}')

    def test_archive_traversal_symlink_and_drive_paths_are_rejected_before_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, mode in [('../escape', 0o100644), ('/absolute', 0o100644),
                               ('C:/escape', 0o100644), ('a\\..\\escape', 0o100644),
                               ('bundle/link', 0o120777)]:
                archive = root / 'bad.zip'
                with zipfile.ZipFile(archive, 'w') as stream:
                    info = zipfile.ZipInfo(name)
                    info.external_attr = mode << 16
                    stream.writestr(info, 'payload')
                destination = root / 'extracted'
                with self.subTest(name=name), self.assertRaises(deps.DependencyError):
                    deps.extract_apbs_archive(archive, destination, 'bundle')
                self.assertFalse(destination.exists())
            self.assertFalse((root / 'escape').exists())

    def test_archive_copies_bundle_contents_and_preserves_executable_bit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / 'good.zip'
            with zipfile.ZipFile(archive, 'w') as stream:
                info = zipfile.ZipInfo('bundle/bin/apbs')
                info.external_attr = 0o100755 << 16
                stream.writestr(info, '#!/bin/sh\nexit 0\n')
                stream.writestr('bundle/share/data.txt', 'needed companion file')
            destination = root / 'unpacked'
            deps.extract_apbs_archive(archive, destination, 'bundle')
            self.assertTrue(os.access(destination / 'bin/apbs', os.X_OK))
            self.assertEqual((destination / 'share/data.txt').read_text(), 'needed companion file')
            self.assertFalse((destination / 'bundle').exists())

    def test_check_from_copied_folder_is_readonly_and_uses_script_location(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'different checkout name'
            (root / 'src/cationize_protein').mkdir(parents=True)
            source = Path(__file__).resolve().parents[1]
            for name in ['__init__.py', 'dependencies.py', 'optional_modelers.py']:
                shutil.copy2(source / 'src/cationize_protein' / name, root / 'src/cationize_protein' / name)
            shutil.copy2(source / 'install.py', root / 'install.py')
            result = subprocess.run([sys.executable, str(root / 'install.py'), '--check', '--json'],
                                    cwd=tmp, capture_output=True, text=True)
            self.assertEqual(result.returncode, 1, result.stderr)
            report = json.loads(result.stdout)
            self.assertEqual(report['project_root'], str(root.resolve()))
            self.assertFalse((root / '.venv').exists())
            self.assertFalse((root / 'src/cationize_protein/__pycache__').exists())
            self.assertIn('selected', report)


if __name__ == '__main__':
    unittest.main()
