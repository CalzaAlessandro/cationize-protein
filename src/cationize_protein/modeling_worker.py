"""Standalone adapters run only in the explicitly selected optional Python.

No cationize-protein imports: MODELLER and folding software can use their own
environments without installing the cationization runtime there.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys


def modeller_model(request):
    from modeller import Environ, log
    from modeller.automodel import AutoModel, assess

    class NativeModel(AutoModel):
        def special_patches(self, aln):
            for i, j in request['disulfides']:
                self.patch(residue_type='DISU', residues=(self.residues[f'{i+1}:A'], self.residues[f'{j+1}:A']))

    log.verbose()
    env = Environ(rand_seed=request['seed'])
    env.io.atom_files_directory = ['.']
    model = NativeModel(env, alnfile='alignment.ali', knowns='template', sequence='target', assess_methods=(assess.DOPE,))
    model.starting_model, model.ending_model = 1, 5
    model.make()
    candidates = [entry for entry in model.outputs if entry.get('failure') is None
                  and math.isfinite(float(entry.get('DOPE score', float('nan'))))]
    if not candidates:
        raise RuntimeError('MODELLER produced no complete, scored candidate')
    best = min(candidates, key=lambda entry: entry['DOPE score'])
    Path('candidate.pdb').write_bytes(Path(best['name']).read_bytes())
    with Path('candidates.csv').open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['pdb','dope_score','selected'])
        for entry in candidates:
            writer.writerow([entry['name'],entry['DOPE score'],entry['name']==best['name']])
    return dict(selected=best['name'], dope_score=float(best['DOPE score']), models=len(candidates),
                model_seed=request['seed'], license_validated=True)


def esmfold_model(request):
    import torch
    import transformers
    from transformers import EsmFold2Model
    device = request['device']
    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    if device == 'cuda' and not torch.cuda.is_available():
        raise RuntimeError('CUDA was requested but is unavailable in the ESMFold2 Python')
    # Transformers 5.20.1 ESMFold2 normalizes activations in float32. Loading
    # every weight in bfloat16 produces a LayerNorm dtype error on real CUDA.
    # Keep the model consistently float32 until a mixed-precision profile is qualified.
    dtype = torch.float32
    torch.manual_seed(request['seed'])
    print(f'ESMFold2: checkpoint={request["checkpoint"]}; device={device}; dtype={dtype}', flush=True)
    print('Loading local cached weights or downloading the public checkpoint on this explicit modeling request.', flush=True)
    model = EsmFold2Model.from_pretrained(request['checkpoint'], revision=request['revision'],
                                         dtype=dtype, device_map=device).eval()
    with torch.inference_mode():
        pdb = model.infer_protein_as_pdb(request['sequence'])
    Path('candidate.pdb').write_text(pdb)
    return dict(torch=torch.__version__, transformers=transformers.__version__, device=device, dtype=str(dtype),
                checkpoint=request['checkpoint'], resolved_revision=getattr(model.config, '_commit_hash', None),
                selection='highest predicted pTM among the diffusion samples', confidence_scale='pLDDT 0-1')


def main():
    request = json.loads(Path(sys.argv[1]).read_text())
    metadata = modeller_model(request) if request['method']=='modeller' else esmfold_model(request)
    Path('worker_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
