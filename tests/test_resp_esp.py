import tempfile,unittest
from pathlib import Path
import numpy as np
from cationize_protein.resp import read_esp
class ESPFormat(unittest.TestCase):
 def test_amber_espgen_optional_zero_flag(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'amber.esp';p.write_text('   2  1    0\n 1 2 3\n 4 5 6\n .25 7 8 9\n')
   xyz,grid,pot=read_esp(p)
   np.testing.assert_array_equal(xyz,[[1,2,3],[4,5,6]])
   np.testing.assert_array_equal(grid,[[7,8,9]])
   np.testing.assert_array_equal(pot,[.25])
 def test_unsupported_multipole_flag_and_infinite_coordinates_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'amber.esp';p.write_text('1 1 1\n 1 2 3\n .25 7 8 9\n')
   with self.assertRaises(ValueError):read_esp(p)
   p.write_text('1 1 0\n nan 2 3\n .25 7 8 9\n')
   with self.assertRaises(ValueError):read_esp(p)
