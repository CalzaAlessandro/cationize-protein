import tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from cationize_protein.parameter_serial import run_serial
from cationize_protein.cli import main

class SerialExecution(unittest.TestCase):
 def test_parameterization_runs_sequentially_without_scheduler(self):
  with tempfile.TemporaryDirectory() as d:
   calls=[]
   with patch('cationize_protein.parameter_jobs.worker',side_effect=lambda p,n:calls.append(n)),patch('cationize_protein.parameter_batch.collect',return_value={'templates':[]}),patch('subprocess.run',side_effect=AssertionError('Scheduler/external process unexpectedly called')):
    result=run_serial(Path(d),['ALT','ART'])
   self.assertEqual(calls,['ALT','ART']);self.assertEqual(result['completed'],calls)
 def test_invalid_or_reused_requests_do_not_start_calculations(self):
  for names in (['ALT','ALT'],['ZZZ'],['ASC']):
   with patch('cationize_protein.parameter_jobs.worker') as worker:
    with self.assertRaises(ValueError):run_serial(Path('.'),names)
    worker.assert_not_called()
 def test_serial_parameterization_cli_does_not_require_apbs_or_slurm(self):
  with patch('cationize_protein.parameter_serial.run_serial',return_value={'completed':['ALT']} ) as run:
   self.assertEqual(main(['parameterize','--project','example','--names','ALT']),0)
   run.assert_called_once_with(Path('example'),['ALT'])
