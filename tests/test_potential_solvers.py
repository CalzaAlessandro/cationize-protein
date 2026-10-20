import unittest,numpy as np
from cationize_protein.potential_solvers import screened_coulomb,pbsa_input,pbsa_to_thermal
class Potentials(unittest.TestCase):
 def test_unscreened_single_charge_and_salt_reduce_potential(self):
  a=screened_coulomb([[0,0,0]],[1],[[10,0,0]],temperature_K=298.,dielectric=78.5,salt_M=0.)
  self.assertAlmostEqual(float(a[0]),.7143,delta=.002)
  b=screened_coulomb([[0,0,0]],[1],[[10,0,0]],temperature_K=298.,dielectric=78.5,salt_M=.15)
  self.assertGreater(b[0],0);self.assertLess(b[0],a[0])
  np.testing.assert_allclose(screened_coulomb([[0,0,0]],[-1],[[10,0,0]],salt_M=0),-a,rtol=1e-10)
 def test_invalid_field_and_singular_point_rejected(self):
  with self.assertRaises(ValueError):screened_coulomb([[0,0,0]],[1],[[0,0,0]])
  with self.assertRaises(ValueError):screened_coulomb([[0,0,0]],[1],[[1,0,0]],salt_M=-1)
 def test_pbsa_units_and_input_are_explicit(self):
  text=pbsa_input(temperature_K=298.,salt_M=.15,spacing_A=.75,solute_dielectric=2.,solvent_dielectric=78.5,ion_radius_A=2.,fillratio=4.)
  for token in ('radiopt=0','npbopt=0','solvopt=1','phiform=2','istrng=150','pbtemp=298'):self.assertIn(token,text)
  self.assertAlmostEqual(pbsa_to_thermal(np.array([.00198720425864083*298]),298)[0],1.)
