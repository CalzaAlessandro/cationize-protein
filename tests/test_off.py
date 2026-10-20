import unittest,tempfile
from pathlib import Path
from cationize_protein.off import read_units
class OFF(unittest.TestCase):
 def test_atoms_bonds_and_positions(self):
  text='''!!index array str
 "TST"
!entry.TST.unit.atoms table str name str type int typex int resx int flags int seq int elmnt dbl chg
 "O" "OW" 0 1 0 1 8 -0.8
 "H1" "HW" 0 1 0 2 -1 0.4
 "H2" "HW" 0 1 0 3 1 0.4
!entry.TST.unit.connectivity table int atom1x int atom2x int flags
 1 2 1
 1 3 1
!entry.TST.unit.positions table dbl x dbl y dbl z
 0 0 0
 1 0 0
 0 1 0
'''
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'test.lib';p.write_text(text);u=read_units(p)['TST']
  self.assertEqual(u['bonds'],[(0,1),(0,2)]);self.assertEqual(u['atoms'][1]['atomic_number'],1)
  self.assertAlmostEqual(sum(a['initial_template_charge_e'] for a in u['atoms']),0)
if __name__=='__main__':unittest.main()
