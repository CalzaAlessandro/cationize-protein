"""Generic graph constraints and deterministic geometric overlap search."""
import unittest
import numpy as np

from cationize_protein.sterics import RotatableBond, derive_rotations, overlap_score, choose_rotamer


class StericTests(unittest.TestCase):
    def test_distal_component_keeps_backbone_fixed(self):
        axes = derive_rotations([(0, 1), (1, 2), (2, 3)], [RotatableBond(1, 2)],
                                ['A:42'] * 4, protected_atoms={0, 1})
        self.assertEqual(axes[0].moving_atoms, (2, 3))

    def test_ring_interresidue_disulfide_and_protected_movement_are_rejected(self):
        cases = [([(0, 1), (1, 2), (2, 0)], ['a'] * 3, set(), ()),
                 ([(0, 1), (1, 2)], ['a', 'a', 'b'], set(), ()),
                 ([(0, 1), (1, 2)], ['a'] * 3, {2}, ()),
                 ([(0, 1), (1, 2)], ['a'] * 3, set(), [(1, 2)])]
        for bonds, residues, protected, disulfides in cases:
            with self.subTest(bonds=bonds, residues=residues, protected=protected, disulfides=disulfides):
                with self.assertRaises(ValueError):
                    derive_rotations(bonds, [RotatableBond(0, 1)], residues,
                                     protected_atoms=protected, disulfide_bonds=disulfides)

    def test_overlap_excludes_through_three_bonds_but_scores_unbonded_atom(self):
        xyz = np.array([[0., 0, 0], [.1, 0, 0], [.2, 0, 0], [.3, 0, 0], [.4, 0, 0]])
        chain = [(0, 1), (1, 2), (2, 3)]
        self.assertEqual(overlap_score(xyz[:4], chain, [1.] * 4, [0]), 0.)
        self.assertGreater(overlap_score(xyz, chain, [1.] * 5, [0]), 1.)

    def test_search_resolves_clash_preserves_bonds_and_protected_coordinates(self):
        # Atom 3 starts on the unrelated atom 4; rotation about 1-2 moves it clear.
        xyz = np.array([[-1., 0, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 0]])
        bonds = [(0, 1), (1, 2), (2, 3)]
        args = (xyz, bonds, [.7] * 5, [RotatableBond(1, 2)])
        kwargs = dict(residue_ids=['a'] * 4 + ['b'], protected_atoms={0, 1, 4}, beam_width=4)
        result, report = choose_rotamer(*args, **kwargs)
        second, second_report = choose_rotamer(*args, **kwargs)
        self.assertLess(report['final_steric_score'], report['initial_steric_score'])
        np.testing.assert_array_equal(result[[0, 1, 4]], xyz[[0, 1, 4]])
        np.testing.assert_array_equal(second, result)
        self.assertEqual(second_report, report)
        for i, j in bonds:
            self.assertAlmostEqual(np.linalg.norm(result[i] - result[j]), np.linalg.norm(xyz[i] - xyz[j]), places=12)
        np.testing.assert_array_equal(xyz[3], [1, 1, 0])

    def test_flat_score_ties_choose_zero_rotation(self):
        xyz = np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0]])
        result, report = choose_rotamer(xyz, [(0, 1), (1, 2)], [0.] * 3,
                                        [RotatableBond(0, 1)], residue_ids=['a'] * 3,
                                        protected_atoms={0})
        np.testing.assert_array_equal(result, xyz)
        self.assertEqual(report['rotations_degrees'], [0])

    def test_invalid_graph_and_degenerate_axis_fail_explicitly(self):
        with self.assertRaises(ValueError):
            derive_rotations([(0, 8)], [RotatableBond(0, 1)], ['a', 'a'])
        with self.assertRaises(ValueError):
            choose_rotamer([[0, 0, 0], [0, 0, 0]], [(0, 1)], [1, 1],
                            [RotatableBond(0, 1)], residue_ids=['a', 'a'])


if __name__ == '__main__':
    unittest.main()
