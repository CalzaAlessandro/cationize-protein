"""Contract tests for supplied-charge PB scoring and OpenDX interpolation."""
import json
import stat
import tempfile
import unittest
from pathlib import Path

import numpy as np

from cationize_protein.electrostatics import (
    APBSConfig, APBSError, DXGrid, build_candidate_grid, rank_candidates,
    read_dx, run_apbs,
)


DX = '''object 1 class gridpositions counts 2 2 2
origin 1 2 3
delta 0 2 0
delta -3 0 0
delta 0 0 4
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
0 1 2
3 4 5
6 7
attribute "dep" string "positions"
'''


class DXTests(unittest.TestCase):
    def test_rotated_orthogonal_grid_interpolates_and_rejects_outside(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'field.dx'
            path.write_text(DX)
            grid = read_dx(path)
        np.testing.assert_allclose(grid.sample([[1, 2, 3], [-.5, 3, 5]]), [0, 3.5])
        with self.assertRaises(ValueError):
            grid.sample([[1, 2, 7.1]])

    def test_truncated_nonfinite_and_skew_grids_are_rejected(self):
        for contents in [DX.replace('6 7', '6'), DX.replace('6 7', '6 nan'),
                         DX.replace('delta -3 0 0', 'delta -3 1 0'),
                         DX.replace('items 8', 'items 9')]:
            with self.subTest(contents=contents), tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / 'field.dx'
                path.write_text(contents)
                with self.assertRaises(ValueError):
                    read_dx(path)

    def test_grid_size_limit_is_checked_before_data_allocation(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'field.dx'
            path.write_text(DX)
            with self.assertRaises(ValueError):
                read_dx(path, max_items=7)


class CandidateTests(unittest.TestCase):
    def test_spheres_keep_chain_group_identity_and_exclude_ion_overlap(self):
        grid = build_candidate_grid({'B:42:side': [[0, 0, 0]],
                                     'A:42:term': [[0, 0, 0], [10, 0, 0]]},
                                    [[0, 0, 0]], [1.0], 1.0, spacing_A=5)
        self.assertEqual(set(grid.indices), {'A:42:term', 'B:42:side'})
        origin = int(np.flatnonzero(np.all(grid.points == 0, axis=1))[0])
        self.assertTrue(grid.blocked[origin])
        self.assertEqual(len(grid.indices['B:42:side']), 7)
        self.assertEqual(len(grid.indices['A:42:term']), 13)
        self.assertTrue(np.all(np.linalg.norm(grid.points[grid.indices['B:42:side']], axis=1) <= 5))

    def test_ranking_uses_accessible_minimum_and_stable_id_ties(self):
        scores = np.array([-100., -3., -3., -1.])
        rows, unavailable = rank_candidates(scores, {'z': [0, 3], 'b': [2], 'a': [1], 'none': [0]}, [True, False, False, False])
        self.assertEqual([row.candidate_id for row in rows], ['a', 'b', 'z'])
        self.assertEqual([row.point_index for row in rows], [1, 2, 3])
        self.assertEqual(unavailable, ['none'])
        self.assertEqual(rows[2].score_kBT_per_e, -1)

    def test_sampling_inputs_reject_invalid_radii_and_nonfinite_values(self):
        for radii, ion in [([-1.], 1.), ([1.], -1.), ([float('nan')], 1.)]:
            with self.subTest(radii=radii, ion=ion), self.assertRaises(ValueError):
                build_candidate_grid({'a': [[0, 0, 0]]}, [[0, 0, 0]], radii, ion)


class APBSRunnerTests(unittest.TestCase):
    def fake_solver(self, directory, log='CONVERGED\nThanks for using APBS!\n', exit_code=0):
        """Small real subprocess double; parser and filesystem effects remain real."""
        path = Path(directory) / 'apbs-test'
        script = '#!/usr/bin/env python3\nfrom pathlib import Path\nimport sys\n'
        script += 'Path("potential.dx").write_text(' + repr(DX) + ')\n'
        script += 'print(' + repr(log) + ')\nsys.exit(' + str(exit_code) + ')\n'
        path.write_text(script)
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    def test_runner_preserves_supplied_charges_config_and_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            exe = self.fake_solver(temp)
            out = Path(temp) / 'run'
            config = APBSConfig(executable=exe, temperature_K=310, solute_dielectric=4,
                                solvent_dielectric=70, salt_M=.05)
            scores, metadata = run_apbs([[0, 0, 0], [2, 0, 0]], [.25, -.5], [1.2, 1.5],
                                        [[1, 2, 3]], out, config,
                                        provenance={'topology_sha256': 'source-hash'}, keep_map=False)
            np.testing.assert_allclose(scores, [0])
            self.assertFalse((out / 'potential.dx').exists())
            self.assertTrue((out / 'PB_COMPLETE.json').is_file())
            self.assertEqual(metadata['protein_charge_e'], -.25)
            self.assertEqual(metadata['inputs']['topology_sha256'], 'source-hash')
            self.assertEqual(metadata['config']['temperature_K'], 310)
            self.assertEqual(metadata['potential_units'], 'kBT/e')
            pqr = (out / 'molecule.pqr').read_text().splitlines()
            self.assertEqual([float(line.split()[-2]) for line in pqr], [.25, -.5])
            input_text = (out / 'apbs.in').read_text()
            self.assertIn('temp 310', input_text)
            self.assertIn('sdie 70', input_text)
            self.assertIn('conc 0.05', input_text)
            self.assertEqual(json.loads((out / 'PB_COMPLETE.json').read_text())['DX_sha256'], metadata['DX_sha256'])

    def test_real_errors_nonconvergence_and_exit_failure_never_write_completion(self):
        cases = [('ERROR: unable to read molecule', 0), ('failed to converge', 0),
                 ('convergence not achieved', 0), ('CONVERGED', 2),
                 ('WARNING: unconverged solution\nThanks for using APBS!', 0),
                 ('Vpmg_fillco: Error setting dielectric map values\nThanks for using APBS!', 0)]
        for log, code in cases:
            with self.subTest(log=log), tempfile.TemporaryDirectory() as temp:
                exe = self.fake_solver(temp, log, code)
                out = Path(temp) / 'run'
                with self.assertRaises(APBSError):
                    run_apbs([[0, 0, 0]], [0], [1.], [[1, 2, 3]], out, APBSConfig(executable=exe))
                self.assertFalse((out / 'PB_COMPLETE.json').exists())

    def test_optional_map_retention_and_benign_error_tolerance_text(self):
        with tempfile.TemporaryDirectory() as temp:
            exe = self.fake_solver(temp, 'Error tolerance = 1e-6\nrelative error tolerance = 1e-6\nCONVERGED')
            out = Path(temp) / 'run'
            run_apbs([[0, 0, 0]], [0], [1.], [[1, 2, 3]], out, APBSConfig(executable=exe), keep_map=True)
            self.assertTrue((out / 'potential.dx').is_file())

    def test_broken_symlink_output_is_not_followed(self):
        with tempfile.TemporaryDirectory() as temp:
            exe = self.fake_solver(temp)
            out = Path(temp) / 'run'
            target = Path(temp) / 'elsewhere'
            out.symlink_to(target, target_is_directory=True)
            with self.assertRaises(FileExistsError):
                run_apbs([[0, 0, 0]], [0], [1.], [[1, 2, 3]], out, APBSConfig(executable=exe))
            self.assertFalse(target.exists())

    def test_invalid_topology_arrays_and_configuration_fail_before_launch(self):
        with self.assertRaises(ValueError):
            APBSConfig(executable='apbs', temperature_K=-1)
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / 'run'
            with self.assertRaises(ValueError):
                run_apbs([[0, 0, 0]], [0, 1], [1], [[1, 2, 3]], out, APBSConfig(executable='missing-apbs'))
            self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
