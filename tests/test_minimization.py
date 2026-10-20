import unittest
from cationize_protein.minimization import group_restraints, amber_input, parse_amber_minimum


class MinimizationTests(unittest.TestCase):
    def test_group_restraints_keep_exact_indices_without_mask_length_limit(self):
        indices=list(range(0,6000,3))+[5999]
        text=group_restraints(indices,10.)
        recovered=[]
        for line in text.splitlines():
            if line.startswith('ATOM '):
                numbers=list(map(int,line.split()[1:]))
                for a,b in zip(numbers[::2],numbers[1::2]):recovered.extend(range(a-1,b))
        self.assertEqual(recovered,sorted(indices))
        self.assertTrue(text.endswith('END\nEND\n'))

    def test_amber_minimization_has_no_dynamics_and_uses_explicit_group_restraints(self):
        text=amber_input([0,3],restraint=10.,tolerance=.1,max_iterations=20000)
        self.assertIn('imin=1',text)
        self.assertIn('igb=8',text)
        self.assertIn('saltcon=0.0',text)
        self.assertIn("restraintmask=''",text)
        self.assertIn('ATOM 1 1 4 4',text)

    def test_amber_parser_requires_final_complete_finite_gradient(self):
        text='FINAL RESULTS\n NSTEP ENERGY RMS GMAX NAME NUMBER\n 88 -123.50 0.010 0.050 CA 5\n Total wall time: 1.0 seconds\n'
        result=parse_amber_minimum(text,.1)
        self.assertEqual(result['steps'],88)
        self.assertEqual(result['energy_kcal_mol'],-123.50)
        self.assertEqual(parse_amber_minimum(text.replace('Total wall time: 1.0 seconds','5. TIMINGS\n| Run done at 15:25:17.077 on 10/20/2026'),.1)['steps'],88)
        for bad in [text.replace('Total wall time:','interrupted'),text.replace('0.010','NaN'),text.replace('0.010','1.5')]:
            with self.subTest(bad=bad),self.assertRaises(RuntimeError):parse_amber_minimum(bad,.1)


if __name__=='__main__':unittest.main()
