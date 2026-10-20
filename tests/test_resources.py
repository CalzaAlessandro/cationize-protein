import unittest
from cationize_protein.resources import default_parameters
from cationize_protein.cli import parser
class BundledParameters(unittest.TestCase):
 def test_source_checkout_finds_manifest_and_default_cli_is_serial(self):
  path=default_parameters();self.assertTrue((path/'manifest.json').is_file())
  args=parser().parse_args(['run','--pdb','protein.pdb','--output','series','--probe-radius','1.369'])
  self.assertEqual(args.parameters,path);self.assertEqual(args.platform,'CPU')
