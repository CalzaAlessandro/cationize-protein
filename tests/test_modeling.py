"""Optional modelers never become required runtime dependencies."""
import contextlib
import importlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT/'examples/pdb_checks/reference_peptide.pdb'


class ModelingTests(unittest.TestCase):
    def module(self):
        try:
            return importlib.import_module('cationize_protein.modeling')
        except ModuleNotFoundError:
            self.fail('Optional native modeling must provide a separate modeling interface')

    def test_missing_method_has_actionable_failure_and_creates_no_output(self):
        module = self.module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root/'sequence.fasta'
            fasta.write_text('>peptide\nDDL\n')
            for method in ('modeller', 'esmfold'):
                with self.subTest(method=method), patch.object(module, 'inspect_modelers', return_value={
                        method: {'available': False, 'python': sys.executable, 'reason': 'not installed'}}):
                    with self.assertRaisesRegex(ValueError, 'not installed'):
                        module.build_model(fasta, root/method, method=method, template=REFERENCE if method=='modeller' else None)
                    self.assertFalse((root/method).exists())

    def test_fasta_and_template_are_validated_before_any_model_download(self):
        module = self.module()
        with tempfile.TemporaryDirectory() as tmp, patch.object(module, 'inspect_modelers', side_effect=AssertionError('No prediction')):
            root = Path(tmp); fasta = root/'sequence.fasta'
            for text in ('>a\nAXL\n', '>a\nDDL\n>b\nDDL\n', 'DDL\n', '>a\n\n'):
                fasta.write_text(text)
                with self.assertRaisesRegex(ValueError, 'FASTA|sequence|canonical'):
                    module.build_model(fasta, root/'out', method='esmfold')
                self.assertFalse((root/'out').exists())
            fasta.write_text('>a\nDDL\n')
            with self.assertRaisesRegex(ValueError, 'template'):
                module.build_model(fasta, root/'out', method='modeller')

    def test_result_is_a_hashed_candidate_with_original_inputs_and_explicit_next_step(self):
        module = self.module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root/'sequence.fasta'; fasta.write_text('>peptide\nDDL\n')
            original = REFERENCE.read_bytes()
            def worker(python, request, directory):
                (directory/'candidate.pdb').write_bytes(original)
                return subprocess.CompletedProcess([], 0, 'external adapter fixture\n')
            with patch.object(module, 'inspect_modelers', return_value={'esmfold': {
                    'available': True, 'python': sys.executable, 'reason': 'code detected'}}), \
                 patch.object(module, '_run_worker', side_effect=worker):
                result = module.build_model(fasta, root/'model', method='esmfold')
            self.assertEqual((root/'model/protein.pdb').read_bytes(), original)
            self.assertEqual((root/'model/sequence.fasta').read_bytes(), fasta.read_bytes())
            record = json.loads((root/'model/MODEL.json').read_text())
            self.assertEqual(record['method'], 'esmfold')
            self.assertEqual(record['status'], 'structural_candidate')
            self.assertEqual(record['sequence'], 'DDL')
            self.assertEqual(record['checkpoint'], 'biohub/ESMFold2-Fast')
            self.assertIn('protein.pdb', record['files_sha256'])
            self.assertIn('prepare', '\n'.join(result['messages']))
            self.assertEqual(REFERENCE.read_bytes(), original)
            self.assertNotIn('modeller', sys.modules)

    def test_wrong_sequence_or_invalid_geometry_never_produces_a_success_marker(self):
        module = self.module()
        for kind in ('sequence', 'geometry'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp); fasta = root/'sequence.fasta'
                fasta.write_text('>peptide\n' + ('DDA' if kind=='sequence' else 'DDL') + '\n')
                raw = REFERENCE.read_text()
                if kind == 'geometry':
                    rows = raw.splitlines()
                    i = next(i for i, row in enumerate(rows) if row.startswith('ATOM') and row[12:16].strip()=='C')
                    rows[i] = rows[i][:30] + ' 999.000' + rows[i][38:]
                    raw = '\n'.join(rows)+'\n'
                def worker(python, request, directory):
                    (directory/'candidate.pdb').write_text(raw)
                    return subprocess.CompletedProcess([], 0, '')
                with patch.object(module, 'inspect_modelers', return_value={'esmfold': {
                        'available': True, 'python': sys.executable}}), patch.object(module, '_run_worker', side_effect=worker):
                    with self.assertRaisesRegex(ValueError, 'sequence|geometry|validation'):
                        module.build_model(fasta, root/'model', method='esmfold')
                self.assertFalse((root/'model/MODEL.json').exists())
                self.assertFalse((root/'model/protein.pdb').exists())
                self.assertTrue((root/'model/MODELING_FAILED.json').exists())

    def test_external_failure_keeps_diagnostic_log_but_never_claims_success(self):
        module = self.module()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); fasta=root/'sequence.fasta'; fasta.write_text('>a\nDDL\n')
            with patch.object(module, 'inspect_modelers', return_value={'esmfold': {
                    'available': True, 'python': sys.executable}}), patch.object(module, '_run_worker', return_value=
                    subprocess.CompletedProcess([], 1, 'GPU allocation failed\n')):
                with self.assertRaisesRegex(ValueError, 'modeling.log'):
                    module.build_model(fasta, root/'model', method='esmfold')
            self.assertIn('GPU allocation failed', (root/'model/modeling.log').read_text())
            self.assertFalse((root/'model/MODEL.json').exists())

    def test_existing_output_and_bad_method_are_rejected(self):
        module=self.module()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); fasta=root/'sequence.fasta'; fasta.write_text('>a\nDDL\n')
            with self.assertRaisesRegex(ValueError, 'method'):
                module.build_model(fasta, root/'new', method='other')
            output=root/'model'; output.mkdir(); (output/'keep').write_text('original')
            with self.assertRaisesRegex(ValueError, 'exist|new directory'):
                module.build_model(fasta, output, method='esmfold')
            self.assertEqual((output/'keep').read_text(), 'original')

    def test_source_disulfides_cannot_disappear_from_a_modeled_candidate(self):
        module=self.module()
        from test_preflight import atom
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'candidate.pdb'
            rows=[]
            for number,offset in ((1,0.),(2,10.)):
                for name,xyz in {'N':(0.,0.,0.),'CA':(1.,0.,0.),'C':(2.,0.,0.),
                                 'O':(3.,0.,0.),'CB':(1.,1.,0.),'SG':(1.,2.,0.)}.items():
                    rows.append(atom(len(rows)+1,name,'CYS',number,'A',
                                     (xyz[0]+offset,xyz[1],xyz[2]), element='S' if name=='SG' else name[0]))
            path.write_text('\n'.join(rows)+'\nEND\n')
            with self.assertRaisesRegex(ValueError, 'disulfide'):
                module._validate_candidate(path, 'CC', [[0,1]])

    def test_nonfinite_disulfide_atom_is_a_reported_validation_error(self):
        module=self.module()
        from test_preflight import atom
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'candidate.pdb'; rows=[]
            for number,offset in ((1,0.),(2,2.05)):
                for name,xyz in {'N':(0.,0.,0.),'CA':(1.,0.,0.),'C':(2.,0.,0.),
                                 'O':(3.,0.,0.),'CB':(1.,1.,0.),'SG':(1.,2.,0.)}.items():
                    coords=(float('nan') if name=='SG' and number==1 else xyz[0]+offset,xyz[1],xyz[2])
                    rows.append(atom(len(rows)+1,name,'CYS',number,'A',coords,
                                     element='S' if name=='SG' else name[0]))
            path.write_text('\n'.join(rows)+'\nEND\n')
            with self.assertRaisesRegex(ValueError, 'validation'):
                module._validate_candidate(path, 'CC', [[0,1]])

    def test_written_disulfide_records_pass_the_same_preparation_checker(self):
        module=self.module()
        from test_preflight import atom
        from cationize_protein.preflight import check_pdb
        coordinates=[{'N':(0.,0.,0.),'CA':(1.45,0.,0.),'C':(2.6,.8,0.),
                      'O':(3.7,.4,0.),'CB':(1.45,-1.5,0.),'SG':(1.45,-3.3,0.)},
                     {'N':(3.8,1.4,0.),'CA':(5.2,1.4,0.),'C':(6.2,2.4,0.),
                      'O':(7.4,2.2,0.),'CB':(4.3,-1.,0.),'SG':(3.5,-3.3,0.)}]
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'candidate.pdb';rows=[]
            for number,atoms in enumerate(coordinates,1):
                for name,xyz in atoms.items():
                    rows.append(atom(len(rows)+1,name,'CYS',number,'A',xyz,
                                     element='S' if name=='SG' else name[0]))
            path.write_text('\n'.join(rows)+'\nEND\n')
            result,_=module._validate_candidate(path,'CC',[[0,1]])
            final=Path(tmp)/'protein.pdb';final.write_text(result)
            checked=check_pdb(final)
            self.assertTrue(checked['ok'],checked['messages'])
            self.assertEqual(checked['ssbond_declarations'][0]['residue_ordinals'],[1,2])

    def test_template_alignment_retains_explicit_sequence_mapping(self):
        module=self.module()
        self.assertEqual(module.template_mapping('DL','DGL'), [0,2])
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            module.template_mapping('DL','DDL')
        with self.assertRaisesRegex(ValueError, 'sequence'):
            module.template_mapping('DDK','DDL')

    def test_cli_defaults_keep_pdbfixer_and_do_not_enable_modeling_implicitly(self):
        from cationize_protein.cli import parser
        choices=parser()._subparsers._group_actions[0].choices
        self.assertIn('model', choices)
        args=parser().parse_args(['model','--method','esmfold','--fasta','native.fasta','--output','model'])
        self.assertEqual(args.esmfold_model, 'fast')
        self.assertEqual(args.device, 'auto')
        prepared=parser().parse_args(['prepare','native.pdb','--output','prepared'])
        self.assertEqual(prepared.protonation, 'propka')


class ModelerDiscoveryTests(unittest.TestCase):
    def module(self):
        try:
            return importlib.import_module('cationize_protein.optional_modelers')
        except ModuleNotFoundError:
            self.fail('Optional modeling installations need read-only discovery')

    def test_invalid_explicit_interpreter_does_not_fall_back(self):
        module=self.module()
        result=module.inspect_modelers(environ={'PATH':'', 'CATIONIZE_MODELLER_PYTHON':'/does/not/exist',
            'CATIONIZE_ESMFOLD_PYTHON':'/also/missing'})
        self.assertFalse(result['modeller']['available'])
        self.assertFalse(result['esmfold']['available'])
        self.assertIn('interpreter', result['modeller']['reason'])

    def test_separate_interpreters_are_inspected_again_without_installation(self):
        module=self.module()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); binary=root/'separate-python'
            binary.write_text('#!/bin/sh\nprintf \'%s\\n\' \'{"available":true,"version":"test-code"}\'\n')
            binary.chmod(0o755)
            env={'PATH':'', 'CATIONIZE_MODELLER_PYTHON':str(binary), 'CATIONIZE_ESMFOLD_PYTHON':str(binary)}
            found=module.inspect_modelers(environ=env)
            self.assertTrue(found['modeller']['available'])
            self.assertTrue(found['esmfold']['available'])
            self.assertIsNone(found['modeller']['license_validated'])
            self.assertIsNone(found['esmfold']['weights_available'])
            binary.unlink()
            missing=module.inspect_modelers(environ=env)
            self.assertFalse(missing['modeller']['available'])

    def test_read_only_probe_does_not_import_optional_engines_in_parent(self):
        module=self.module()
        before=set(sys.modules)
        found=module.inspect_modelers()
        self.assertIn('modeller', found)
        self.assertIn('esmfold', found)
        self.assertFalse({'modeller','torch','transformers'} & (set(sys.modules)-before))

    def test_broken_external_interpreter_report_cannot_break_the_main_installation(self):
        module=self.module()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); binary=root/'broken-python'
            binary.write_text('#!/bin/sh\nprintf \'%s\\n\' \'[]\'\n'); binary.chmod(0o755)
            found=module.inspect_modelers(environ={'PATH':'',
                'CATIONIZE_MODELLER_PYTHON':str(binary),'CATIONIZE_ESMFOLD_PYTHON':str(binary)})
            self.assertFalse(found['modeller']['available'])
            self.assertFalse(found['esmfold']['available'])
            self.assertIn('probe',found['esmfold']['reason'])

    def test_optional_modeler_can_be_found_in_a_different_python_on_path(self):
        module=self.module()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); binary=root/'python3.15'
            binary.write_text('#!/bin/sh\nprintf \'%s\\n\' \'{"available":true,"version":"separate-runtime"}\'\n')
            binary.chmod(0o755)
            found=module.inspect_modelers(python='/missing/default-python',environ={'PATH':str(root)})
            self.assertTrue(found['modeller']['available'])
            self.assertEqual(found['modeller']['python'],str(binary))

    def test_discovered_python_stays_usable_after_the_worker_changes_directory(self):
        module=self.module()
        import os
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); binary=root/'runtime/bin/python3.15';binary.parent.mkdir(parents=True)
            binary.write_text('#!/bin/sh\nprintf \'%s\\n\' \'{"available":true,"version":"relative-runtime"}\'\n')
            binary.chmod(0o755); staging=root/'worker';staging.mkdir()
            previous=Path.cwd()
            try:
                os.chdir(root)
                found=module.inspect_modelers(python='/missing/default-python',environ={'PATH':'runtime/bin'})
                result=subprocess.run([found['modeller']['python'],'-c','ignored'],cwd=staging,
                                      capture_output=True,text=True)
                self.assertEqual(result.returncode,0)
                self.assertTrue(json.loads(result.stdout)['available'])
            finally:os.chdir(previous)


class WorkerContractTests(unittest.TestCase):
    """API fixtures exercise adapters; they are not physical folding tests."""
    def worker(self):
        try:
            return importlib.import_module('cationize_protein.modeling_worker')
        except ModuleNotFoundError:
            self.fail('Optional external environments need a standalone adapter')

    def test_esmfold_uses_the_new_api_and_records_resolved_checkpoint(self):
        worker=self.worker()
        raw=REFERENCE.read_text()
        class Model:
            config=SimpleNamespace(_commit_hash='fixed-model-revision')
            @classmethod
            def from_pretrained(cls, name, **kwargs):
                if name != 'biohub/ESMFold2-Fast' or kwargs.get('revision') != 'pinned':
                    raise AssertionError('Wrong checkpoint')
                return cls()
            def eval(self): return self
            def infer_protein_as_pdb(self, sequence):
                if sequence != 'DDL': raise AssertionError('Wrong sequence')
                return raw
        torch=SimpleNamespace(__version__='fixture', cuda=SimpleNamespace(is_available=lambda:False),
            float32='float32', manual_seed=lambda seed:None, inference_mode=contextlib.nullcontext)
        transformers=SimpleNamespace(__version__='fixture',EsmFold2Model=Model)
        with tempfile.TemporaryDirectory() as tmp, patch.dict(sys.modules, {'torch':torch,'transformers':transformers}):
            previous=Path.cwd()
            try:
                import os
                os.chdir(tmp)
                with contextlib.redirect_stdout(io.StringIO()):
                    metadata=worker.esmfold_model(dict(device='auto',seed=0,checkpoint='biohub/ESMFold2-Fast',
                                                      revision='pinned',sequence='DDL'))
                self.assertEqual(Path('candidate.pdb').read_text(),raw)
                self.assertEqual(metadata['device'],'cpu')
                self.assertEqual(metadata['resolved_revision'],'fixed-model-revision')
            finally:
                os.chdir(previous)

    def test_modeller_selects_the_best_finite_successful_dope_candidate(self):
        worker=self.worker()
        class AutoModel:
            def __init__(self,*args,**kwargs):
                self.outputs=[]
            def make(self):
                for name,score,failure in [('first.pdb',-5,None),('best.pdb',-20,None),('failed.pdb',-50,'failed')]:
                    Path(name).write_text(name)
                    self.outputs.append({'name':name,'DOPE score':score,'failure':failure})
        modeller=SimpleNamespace(Environ=lambda **kwargs:SimpleNamespace(io=SimpleNamespace()),
                                 log=SimpleNamespace(verbose=lambda:None))
        automodel=SimpleNamespace(AutoModel=AutoModel,assess=SimpleNamespace(DOPE='DOPE'))
        with tempfile.TemporaryDirectory() as tmp, patch.dict(sys.modules, {'modeller':modeller,'modeller.automodel':automodel}):
            previous=Path.cwd()
            try:
                import os
                os.chdir(tmp)
                metadata=worker.modeller_model({'disulfides':[],'seed':-99})
                self.assertEqual(Path('candidate.pdb').read_text(),'best.pdb')
                self.assertEqual(metadata['dope_score'],-20)
                self.assertEqual(metadata['models'],2)
                self.assertEqual(metadata['model_seed'],-99)
                self.assertIn('best.pdb',Path('candidates.csv').read_text())
            finally:
                os.chdir(previous)

    def test_esmfold_cuda_keeps_consistent_float32_for_normalization(self):
        worker=self.worker()
        class Model:
            config=SimpleNamespace(_commit_hash='revision')
            @classmethod
            def from_pretrained(cls,name,**kwargs):
                if kwargs['dtype'] != 'float32':
                    raise AssertionError('The tested ESMFold2 profile requires float32 normalization weights')
                return cls()
            def eval(self):return self
            def infer_protein_as_pdb(self,sequence):return REFERENCE.read_text()
        torch=SimpleNamespace(__version__='fixture',cuda=SimpleNamespace(
            is_available=lambda:True,is_bf16_supported=lambda:True),float32='float32',bfloat16='bfloat16',
            manual_seed=lambda seed:None,inference_mode=contextlib.nullcontext)
        transformers=SimpleNamespace(__version__='fixture',EsmFold2Model=Model)
        with tempfile.TemporaryDirectory() as tmp,patch.dict(sys.modules,{'torch':torch,'transformers':transformers}):
            previous=Path.cwd()
            try:
                import os
                os.chdir(tmp)
                with contextlib.redirect_stdout(io.StringIO()):
                    metadata=worker.esmfold_model(dict(device='cuda',seed=0,checkpoint='biohub/ESMFold2-Fast',
                                                      revision='pinned',sequence='DDL'))
                self.assertEqual(metadata['dtype'],'float32')
                self.assertTrue(Path('candidate.pdb').is_file())
            finally:os.chdir(previous)
