"""Exclusive, restart-aware Gaussian CPU model jobs with immutable attempts."""
import hashlib,json,os,shutil,subprocess,tempfile,time
from pathlib import Path

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def completed(source,output):
 source,output=Path(source),Path(output);expected=digest(source/'opt_esp.gjf')
 for p in sorted(output.glob('attempt_*/QM_COMPLETE.json'),reverse=True):
  r=json.loads(p.read_text())
  if r['input_sha256']==expected and digest(p.parent/'gaussian.log')==r['log_sha256']:return p.parent
 return None

def run(source,output,executable='g16'):
 source,output=Path(source).resolve(),Path(output).resolve()
 old=completed(source,output)
 if old:return old
 output.mkdir(parents=True,exist_ok=True)
 count=max([int(p.name.split('_')[-1]) for p in output.glob('attempt_*') if p.is_dir()]+[0])+1
 attempt=output/f'attempt_{count:03d}';attempt.mkdir();shutil.copy2(source/'opt_esp.gjf',attempt/'input.gjf')
 binary=shutil.which(executable)
 if not binary:raise FileNotFoundError(executable)
 env=os.environ.copy();start=time.time()
 with tempfile.TemporaryDirectory(prefix='cationize_qm_') as scratch:
  env['GAUSS_SCRDIR']=scratch
  with (attempt/'input.gjf').open() as inp,(attempt/'gaussian.log').open('x') as log:
   result=subprocess.run([binary],cwd=attempt,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT)
 text=(attempt/'gaussian.log').read_text(errors='replace')
 valid=result.returncode==0 and all(s in text for s in ('Normal termination of Gaussian','Optimization completed.','ESP Fit Center'))
 record=dict(returncode=result.returncode,wall_seconds=time.time()-start,input_sha256=digest(source/'opt_esp.gjf'),log_sha256=digest(attempt/'gaussian.log'),executable=binary,job_id=os.environ.get('SLURM_JOB_ID'),partition=os.environ.get('SLURM_JOB_PARTITION'),validated=valid)
 (attempt/('QM_COMPLETE.json' if valid else 'QM_FAILED.json')).write_text(json.dumps(record,indent=2)+'\n')
 if not valid:raise RuntimeError(f'Gaussian/ESP failed: {attempt}')
 return attempt
