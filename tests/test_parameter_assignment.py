import importlib.util
import unittest


class AssignmentAPITests(unittest.TestCase):
    def test_unsupported_force_field_never_falls_back_to_ff14(self):
        from cationize_protein.parameter_assignment import parameterize
        with self.assertRaisesRegex(ValueError, 'Unknown protein force field'):
            parameterize(None, '.', force_field='OpenFF')

    @unittest.skipUnless(importlib.util.find_spec('openmm') is not None and importlib.util.find_spec('parmed') is not None,'Numerical runtime required')
    def test_missing_unrelated_product_template_does_not_block_native_alternate_field(self):
        from test_native_states import native_peptide
        from cationize_protein.parameter_assignment import parameterize
        from cationize_protein.resources import default_parameters
        native=native_peptide('HIE')
        for field in ('ff19SB','ff99SBildn'):
            with self.subTest(force_field=field):
                typed,system=parameterize(native,default_parameters(),force_field=field)
                self.assertEqual(len(typed.atoms),system.getNumParticles())


if __name__ == '__main__': unittest.main()
