import unittest
from cationize_protein.resp import safe_ac_names,validate_resp_elements
class Elements(unittest.TestCase):
 def test_surrogate_names_preserve_elements(self):
  atoms=[dict(name='S1',atomic_number=7),dict(name='T2',atomic_number=6),dict(name='HG',atomic_number=1)]
  self.assertEqual(safe_ac_names(atoms),['N001','C002','H003'])
 def test_wrong_atomic_number_is_rejected(self):
  text='header\n &cntrl\n nmol=2\n &end\n 1.0\nmodel\n 2 3\n 7 0\n 6 0\n 1 0\n 1.0\nmodel\n 2 3\n 7 0\n 6 0\n 1 0\n'
  validate_resp_elements(text,[7,6,1],2)
  with self.assertRaises(ValueError):validate_resp_elements(text.replace(' 7 0',' 16 0',1),[7,6,1],2)
  with self.assertRaises(ValueError):validate_resp_elements(text,[7,6,1],1)
if __name__=='__main__':unittest.main()
