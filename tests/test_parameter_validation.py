import unittest
import numpy as np
try:
 import parmed
except ImportError:
 parmed=None
from cationize_protein.parameter_validation import capped_structure

@unittest.skipIf(parmed is None,'ParmEd unavailable')
class CappedStructure(unittest.TestCase):
 def mapping(self):
  return dict(atoms=[dict(index=1,name='N',atomic_number=7,amber_type='N',residue='ALT',residue_index=2),dict(index=2,name='C',atomic_number=6,amber_type='C',residue='ACE',residue_index=1),dict(index=3,name='T1',atomic_number=7,amber_type='N',residue='ALT',residue_index=2)],bonds=[[1,2],[1,3]],bond_orders=[1.,1.])
 def test_group_residues_without_losing_coordinate_atom_correspondence(self):
  xyz=np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
  s,original=capped_structure(self.mapping(),xyz,[.1,.2,.3])
  self.assertEqual([r.name for r in s.residues],['ACE','ALT'])
  self.assertEqual([a.name for a in s.atoms],['C','N','T1'])
  self.assertEqual(original,[2,1,3])
  np.testing.assert_array_equal(s.coordinates,xyz[[1,0,2]])
  np.testing.assert_allclose([a.charge for a in s.atoms],[.2,.1,.3])
  self.assertEqual({tuple(sorted((b.atom1.idx,b.atom2.idx))) for b in s.bonds},{(0,1),(1,2)})
  self.assertTrue(all(a.mass>0 for a in s.atoms))
 def test_invalid_mapping_or_dimensions_rejected(self):
  m=self.mapping()
  with self.assertRaises(ValueError):capped_structure(m,np.zeros((2,3)),[0.,0.,0.])
  m['atoms'][2]['index']=2
  with self.assertRaises(ValueError):capped_structure(m,np.zeros((3,3)),[0.,0.,0.])
