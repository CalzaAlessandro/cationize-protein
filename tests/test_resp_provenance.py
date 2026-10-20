import hashlib,tempfile,unittest
from pathlib import Path
from cationize_protein.parameter_provenance import known_resp_source
from cationize_protein.resp import esp_header,rename_off_unit
class SourceAndFormats(unittest.TestCase):
 def test_fixed_width_espgen_large_grid_header(self):
  self.assertEqual(esp_header('   4210179    0'),(42,10179))
  self.assertEqual(esp_header('   5011773    0'),(50,11773))
  self.assertEqual(esp_header('   41 9440    0'),(41,9440))
 def test_unit_variable_rename_preserves_numeric_data(self):
  text='!!index array str\n "u"\n!entry.u.unit.atoms table str name\n "N" "N" 0 1 0 1 7 -0.4\n!entry.u.unit.positions table dbl x dbl y dbl z\n 0.1 0.2 0.3\n'
  renamed=rename_off_unit(text,'u','SET')
  self.assertIn('!entry.SET.',renamed);self.assertIn(' "SET"\n',renamed)
  self.assertEqual(rename_off_unit(renamed,'SET','u'),text)
 def test_historical_source_requires_actual_hash_matching_bytes(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);a=p/'parameters/provenance';a.mkdir(parents=True)
   data=b'original scientific fitting code';sha=hashlib.sha256(data).hexdigest();f=a/f'resp_{sha}.py';f.write_bytes(data)
   self.assertTrue(known_resp_source(p,sha));f.write_bytes(b'changed');self.assertFalse(known_resp_source(p,sha))
