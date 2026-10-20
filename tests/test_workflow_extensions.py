"""Public choices and charge-aware protonation integration."""
from dataclasses import replace
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

from cationize_protein.cli import parser
from cationize_protein.workflow import WorkflowConfig, validate_inputs
from test_workflow import WorkflowTests
from test_topology import structure, atom_named
from cationize_protein.topology import inspect_structure, UnsupportedTopology


class ExtensionTests(unittest.TestCase):
    def test_declared_disulfide_must_exist_in_prepared_pdb_graph(self):
        from cationize_protein.workflow import _verify_declared_disulfides,WorkflowError
        residues=[NS(index=0),NS(index=1)]
        a,b=[NS(name='SG',residue=r) for r in residues]
        report={'ssbond_declarations':[{'residue_ordinals':[1,2]}]}
        pdb=NS(topology=NS(bonds=lambda:[]))
        with self.assertRaisesRegex(WorkflowError,'disulfide'):_verify_declared_disulfides(pdb,report)
        pdb.topology.bonds=lambda:[(a,b)]
        _verify_declared_disulfides(pdb,report)
    def test_stage_export_recovers_after_copy_interruption_without_exposing_partial_stage(self):
        import hashlib,shutil
        from cationize_protein.workflow import _export_stage
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);checkpoint=root/'internal/steps/a/checkpoint';checkpoint.mkdir(parents=True)
            (checkpoint/'protein.pdb').write_text('coordinates')
            record=dict(protein_charge_e=0.,artifacts={'protein.pdb':hashlib.sha256(b'coordinates').hexdigest()})
            (checkpoint/'STEP_COMPLETE.json').write_text(json.dumps(record))
            with patch('cationize_protein.workflow.shutil.copy2',side_effect=OSError('interrupted copy')):
                with self.assertRaises(OSError):_export_stage(root,'25',checkpoint,{'step':1})
            self.assertFalse((root/'25').exists())
            _export_stage(root,'25',checkpoint,{'step':1})
            self.assertTrue((root/'25/STAGE_COMPLETE.json').is_file())

    def test_baseline_failure_is_preserved_and_automatically_retried(self):
        import hashlib
        from cationize_protein.workflow import run
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(WorkflowTests().config(Path(tmp)),percentages=(0,))
            residue=NS(ordinal=1,modified_groups=(),native_atom_indices={},source_spec=None)
            info=NS(residues=[residue],candidates=[],coordinates=[],protected_atoms=[],charges_e=NS(sum=lambda:0.))
            attempts=[]
            def save(s,system,directory,metadata,**kw):
                directory.mkdir(parents=True);attempts.append(directory)
                if len(attempts)==1:
                    (directory/'partial.txt').write_text('preserve this');raise OSError('interrupted baseline')
                (directory/'protein.pdb').write_text('coordinates')
                record=dict(metadata,artifacts={'protein.pdb':hashlib.sha256(b'coordinates').hexdigest()})
                (directory/'STEP_COMPLETE.json').write_text(json.dumps(record));return record
            def read(directory):return NS(),NS(),json.loads((directory/'STEP_COMPLETE.json').read_text())
            with patch('cationize_protein.workflow._initial_structure',return_value=(NS(),NS())), \
                 patch('cationize_protein.workflow.inspect_structure',return_value=info), \
                 patch('cationize_protein.workflow._save_checkpoint',side_effect=save), \
                 patch('cationize_protein.workflow._read_checkpoint',side_effect=read):
                with self.assertRaises(OSError):run(cfg)
                result=run(cfg)
            self.assertEqual(result['final_modified_groups'],0)
            self.assertEqual(len(list(cfg.output_dir.rglob('partial.txt'))),1)

    def test_resume_rejects_parameter_bundle_changes_even_before_template_first_use(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg=WorkflowTests().config(Path(tmp))
            signature=validate_inputs(cfg);cfg.output_dir.mkdir()
            (cfg.output_dir/'RUN.json').write_text(json.dumps(signature))
            (cfg.parameters_dir/'manifest.json').write_text('{"templates":[],"revision":2}')
            from cationize_protein.workflow import WorkflowError
            with self.assertRaises(WorkflowError):validate_inputs(cfg)
    def test_cli_accepts_explicit_physics_and_exports(self):
        args=parser().parse_args(['run','--pdb','input.pdb','--output','out','--probe-radius','1.369',
            '--force-field','ff19SB','--engine','openmm','--solver','coulomb',
            '--formats','amber','gromacs','--ph','7.4'])
        self.assertEqual(args.force_field,'ff19SB')
        self.assertEqual(args.formats,['amber','gromacs'])

    def test_coulomb_does_not_require_any_external_solver_and_fingerprints_choices(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg=WorkflowTests().config(Path(tmp))
            cfg=replace(cfg,solver='coulomb',force_field='ff19SB',protonation='fixed')
            with patch('cationize_protein.workflow.shutil.which',return_value=None):
                signature=validate_inputs(cfg)
            self.assertEqual(signature['force_field'],'ff19SB')
            self.assertEqual(signature['solver'],'coulomb')
            self.assertFalse(cfg.output_dir.exists())

    def test_protonated_acid_retains_identity_and_targetable_group(self):
        for state,oxygen,hydrogen in [('ASH','OD2','HD2'),('GLH','OE2','HE2')]:
            with self.subTest(state=state):
                s=structure([['ALA',state,'LEU']]); r=s.residues[1]
                h=NS(name=hydrogen,idx=len(s.atoms),residue=r,atomic_number=1,charge=.1,solvent_radius=1.2,rmin=.1)
                r.atoms.append(h);s.atoms.append(h);s.bonds.append(NS(atom1=atom_named(s,1,oxygen),atom2=h))
                import numpy as np
                s.coordinates=np.vstack([s.coordinates,[[1.,2.,3.]]])
                info=inspect_structure(s)
                self.assertEqual(info.residues[1].canonical_name,state)
                self.assertEqual([c.group for c in info.candidates if c.residue_ordinal==2],['side'])
                # A protonated acid with missing H must not silently become COO-.
                s.bonds.pop()
                with self.assertRaises(UnsupportedTopology):inspect_structure(s)

    def test_protonation_rebuild_minimize_then_audit_without_counting_as_another_chemical_step(self):
        from cationize_protein.workflow import _protonate
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(WorkflowTests().config(Path(tmp)),protonation='propka')
            s=NS(); min_s=NS(); events=[]
            info=NS(residues=[NS(ordinal=1,context='internal')])
            records=[dict(applicable_native_states={'1':'HIP'},proposed_native_states={'1':'HIP'},decisions=[]),
                     dict(applicable_native_states={},proposed_native_states={},decisions=[])]
            def audit(*a,**kw):events.append('audit');return records.pop(0)
            def apply(*a,**kw):events.append('rebuild');return s,NS(),{'expected_charge_change_e':1}
            def minimize(*a,**kw):events.append('minimize');return min_s,{'RMS_force_kcal_mol_A':.01}
            with patch('cationize_protein.protonation.audit',side_effect=audit), \
                 patch('cationize_protein.native_states.apply_native_states',side_effect=apply), \
                 patch('cationize_protein.native_states.supported_native_contexts',return_value={'HIP':['internal']}), \
                 patch('cationize_protein.workflow.inspect_structure',return_value=info), \
                 patch('cationize_protein.workflow._minimize',side_effect=minimize), \
                 patch('cationize_protein.workflow._reference_positions',return_value=[]):
                result,record=_protonate(s,Path(tmp)/'audit',cfg,{},previous={'proposed_native_states':{'1':'HIP'}})
            self.assertIs(result,min_s)
            self.assertEqual(events,['audit','rebuild','minimize','audit'])
            self.assertEqual(record['state_changes']['expected_charge_change_e'],1)
            self.assertEqual(record['history']['proposed_native_states'],{})


if __name__=='__main__':unittest.main()
