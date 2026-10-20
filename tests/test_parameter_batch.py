import hashlib,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from cationize_protein.parameter_batch import template_state, collect

class BatchTests(unittest.TestCase):
 def test_missing_fit_is_not_validated(self):
  with tempfile.TemporaryDirectory() as d:
   self.assertEqual(template_state(Path(d),'GLX')['status'],'awaiting_QM_RESP_validation')
 def test_incomplete_library_without_marker_is_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'parameters').mkdir();(p/'parameters/GLX.lib').write_text('incomplete')
   r=template_state(p,'GLX');self.assertEqual(r['status'],'incomplete_or_invalid');self.assertIn('marker',r['reason'])
 def test_reused_template_hash_must_match_and_no_status_is_silently_promoted(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'parameters').mkdir();lib=p/'parameters/ASC.lib';lib.write_text('validated local fit')
   manifest=dict(schema=1,templates=[dict(name='ASC',reused_from='ASL',status='reused_validated_HSA_template',library_sha256=hashlib.sha256(lib.read_bytes()).hexdigest())])
   (p/'parameters/manifest.json').write_text(json.dumps(manifest))
   r=collect(p);self.assertEqual(r['templates'][0]['status'],'reused_validated_HSA_template')
   lib.write_text('changed')
   r=collect(p);self.assertEqual(r['templates'][0]['status'],'incomplete_or_invalid')
