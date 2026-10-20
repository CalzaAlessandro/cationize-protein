import unittest
from cationize_protein.catalog import templates,resolve
class Catalog(unittest.TestCase):
 def test_states(self):
  specs=templates();self.assertEqual(len(specs),32);self.assertEqual(len({x['name'] for x in specs}),32)
  self.assertEqual(resolve('ASP','internal',{'side'})['name'],'ASC')
  self.assertEqual(resolve('GLU','n_terminal',{'side'})['name'],'GLX')
  self.assertEqual(resolve('ASP','c_terminal',{'side'})['name'],'ASZ')
  self.assertEqual(resolve('ASP','c_terminal',{'terminal'})['name'],'AST')
  self.assertEqual(resolve('ASP','c_terminal',{'side','terminal'})['name'],'ASD')
  self.assertEqual(resolve('HIP','c_terminal',{'terminal'})['formal_charge'],2)
  self.assertEqual(resolve('CYM','c_terminal',{'terminal'})['formal_charge'],0)
 def test_unsupported(self):
  with self.assertRaises(ValueError):resolve('ALA','internal',{'side'})
  with self.assertRaises(ValueError):resolve('LEU','both_terminal',{'terminal'})
if __name__=='__main__':unittest.main()
