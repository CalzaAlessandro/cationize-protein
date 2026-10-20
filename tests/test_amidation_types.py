import unittest
from cationize_protein.models import amidation_atom_types
class AmidationTypes(unittest.TestCase):
 def test_carboxylate_becomes_amide_not_an_anionic_carbonyl(self):
  atoms=[dict(name='CD',atomic_number=6,amber_type='CO'),dict(name='OE1',atomic_number=8,amber_type='O2'),dict(name='CG',atomic_number=6,amber_type='2C')]
  amidation_atom_types(atoms,0,1)
  self.assertEqual([a['amber_type'] for a in atoms],['C','O','2C'])
 def test_wrong_elements_are_rejected(self):
  with self.assertRaises(ValueError):amidation_atom_types([dict(atomic_number=7),dict(atomic_number=8)],0,1)
