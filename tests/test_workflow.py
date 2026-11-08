"""Progressive group counts, immutable artifacts and precise library gates."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

from cationize_protein.electrostatics import APBSConfig
from cationize_protein.workflow import (
    WorkflowConfig, WorkflowError, stage_counts, validate_inputs,
    require_validated_template, verify_artifacts, inherited_atom_map, run, _rotatable_bonds, _pdb_structure,
)


class WorkflowTests(unittest.TestCase):
    def test_group_counts_use_decimal_half_up_including_independent_terminal_group(self):
        self.assertEqual(stage_counts(3), [('0', 0), ('25', 1), ('50', 2), ('75', 2), ('100', 3)])
        self.assertEqual(stage_counts(5, ['10', '30', '50', '70', '90']),
                         [('0', 0), ('10', 1), ('30', 2), ('50', 3), ('70', 4), ('90', 5)])
        self.assertEqual(stage_counts(0, [100]), [('0', 0), ('100', 0)])

    def test_custom_percentages_sort_and_reject_duplicates_or_invalid_input(self):
        self.assertEqual(stage_counts(10, [100, 12.5]), [('0', 0), ('12.5', 1), ('100', 10)])
        for percentages in [[25, '25.0'], [-1], [101], ['NaN']]:
            with self.subTest(percentages=percentages), self.assertRaises(ValueError):
                stage_counts(4, percentages)
        with self.assertRaises(ValueError):
            stage_counts(-1)

    def config(self, root):
        top, rst, params = root / 'input.parm7', root / 'input.rst7', root / 'parameters'
        top.write_text('topology')
        rst.write_text('coordinates')
        params.mkdir()
        (params / 'manifest.json').write_text('{"templates": []}')
        return WorkflowConfig(output_dir=root / 'output', parameters_dir=params,
                              apbs=APBSConfig(executable=sys.executable), probe_radius_A=1.369,
                              topology_path=top, coordinates_path=rst,protonation='fixed')

    def test_input_guards_are_readonly_and_require_one_complete_source(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cfg = self.config(root)
            record = validate_inputs(cfg)
            self.assertIn('topology_sha256', record['inputs'])
            self.assertFalse(cfg.output_dir.exists())
            with self.assertRaises(ValueError):
                validate_inputs(replace(cfg, coordinates_path=None))
            with self.assertRaises(ValueError):
                validate_inputs(replace(cfg, pdb_path=cfg.topology_path))
            with self.assertRaises(ValueError):
                validate_inputs(replace(cfg, probe_radius_A=-1))
            self.assertFalse(cfg.output_dir.exists())

    def test_matching_run_is_idempotent_and_changed_inputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg = self.config(Path(temp))
            signature = validate_inputs(cfg)
            cfg.output_dir.mkdir()
            (cfg.output_dir / 'RUN.json').write_text(json.dumps(signature))
            self.assertEqual(validate_inputs(cfg), signature)
            cfg.topology_path.write_text('changed topology')
            with self.assertRaises(WorkflowError):
                validate_inputs(cfg)

    def test_existing_unknown_directory_and_symlinks_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg = self.config(Path(temp))
            cfg.output_dir.symlink_to(Path(temp) / 'uncreated', target_is_directory=True)
            with self.assertRaises(WorkflowError):
                validate_inputs(cfg)
            self.assertFalse((Path(temp) / 'uncreated').exists())

    def test_selected_pending_template_fails_precisely_without_substitution(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            spec = dict(name='GLX', parent='GLU', context='n_terminal', modified_groups=['side'])
            (root / 'manifest.json').write_text(json.dumps({'templates': [dict(spec, status='awaiting_QM_RESP_validation')]}))
            candidate = NS(candidate_id='r1:chain=A:group=side', target_spec=spec)
            with self.assertRaisesRegex(WorkflowError, 'r1:chain=A:group=side.*GLX.*awaiting_QM_RESP_validation'):
                require_validated_template(root, candidate)

    def test_validated_library_requires_matching_content_hash(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            lib = root / 'ASC.lib'
            lib.write_text('library bytes')
            spec = dict(name='ASC', parent='ASP', context='internal', modified_groups=['side'])
            entry = dict(spec, status='reused_validated_HSA_template', library_sha256=hashlib.sha256(lib.read_bytes()).hexdigest())
            (root / 'manifest.json').write_text(json.dumps({'templates': [entry]}))
            candidate = NS(candidate_id='a', target_spec=spec)
            self.assertEqual(require_validated_template(root, candidate), lib)
            lib.write_text('modified')
            with self.assertRaisesRegex(WorkflowError, 'hash'):
                require_validated_template(root, candidate)

    def test_invalid_prepared_chemistry_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg = self.config(Path(temp))
            with patch('cationize_protein.workflow._initial_structure', side_effect=WorkflowError('unsupported chemistry')):
                with self.assertRaisesRegex(WorkflowError, 'unsupported chemistry'):
                    run(cfg)
            self.assertFalse(cfg.output_dir.exists())

    def test_native_pdb_rejects_unsupported_protonation_before_normalization(self):
        for name in ['ASH', 'GLH', 'CYM']:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temp:
                pdb = Path(temp) / 'input.pdb'
                pdb.write_text(f'ATOM  {1:5d} {"N":>4s} {name:>3s} A{1:4d}    0.000   0.000   0.000\n')
                with self.assertRaisesRegex(WorkflowError, name):
                    _pdb_structure(pdb, Path(temp))

    def test_guanidinium_resonance_bonds_are_never_rotatable(self):
        # CA-CD-NE-CZ(NH1)(NH2), with one terminal H on each N.
        elements = [6, 6, 7, 6, 7, 7, 1, 1]
        atoms = [NS(idx=i, atomic_number=z, residue_ordinal=1) for i, z in enumerate(elements)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (3, 5), (4, 6), (5, 7)]
        structure = NS(bonds=[NS(atom1=atoms[i], atom2=atoms[j], order=1.) for i, j in edges])
        info = NS(residues=[NS(atom_indices=tuple(range(8)))], atoms=atoms,
                  bonds=edges, protected_atoms=(0,), disulfide_bonds=())
        axes = _rotatable_bonds(structure, info, 1)
        selected = {frozenset((axis.base, axis.tip)) for axis in axes}
        self.assertIn(frozenset((1, 2)), selected)
        self.assertTrue(selected.isdisjoint({frozenset((2, 3)), frozenset((3, 4)), frozenset((3, 5))}))

    def test_inherited_atom_mapping_keeps_native_hydrogens_across_legacy_names(self):
        native_names = [('N', 7), ('H', 1), ('CA', 6), ('HA', 1), ('CB', 6), ('HB2', 1), ('HB3', 1)]
        product_names = [('N1', 7), ('H1', 1), ('C1', 6), ('H2', 1), ('C3', 6), ('H3', 1), ('H4', 1)]
        residue = NS(atoms=[])
        old = [NS(name=name, atomic_number=z, residue=residue) for name, z in native_names]
        new = [NS(name=name, atomic_number=z) for name, z in product_names]
        residue.atoms = old
        links = [(0, 1), (0, 2), (2, 3), (2, 4), (4, 5), (4, 6)]
        st = NS(atoms=old, residues=[residue], bonds=[NS(atom1=old[i], atom2=old[j]) for i, j in links])
        template = NS(atoms=new, bonds=[NS(atom1=new[i], atom2=new[j]) for i, j in links])
        state = NS(ordinal=1, native_atom_indices={name: i for i, (name, _) in enumerate(native_names)})
        mapping = inherited_atom_map(st, state, template, {'reused_from': 'ASL'})
        self.assertEqual(mapping, {'N': 'N1', 'CA': 'C1', 'CB': 'C3', 'H': 'H1', 'HA': 'H2', 'HB2': 'H3', 'HB3': 'H4'})

    def test_baseline_is_exported_before_a_later_step_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg = self.config(Path(temp))
            cfg = replace(cfg, percentages=(100,))
            residue = NS(ordinal=1, modified_groups=(), native_atom_indices={}, source_spec=None)
            info = NS(residues=[residue], candidates=[NS()], coordinates=[], protected_atoms=[], charges_e=NS(sum=lambda: 0.))
            def save(structure, system, directory, metadata, **kwargs):
                directory.mkdir(parents=True)
                hashes = {}
                for name in ['protein.pdb', 'system.xml', 'state.json']:
                    (directory / name).write_text(name)
                    hashes[name] = hashlib.sha256(name.encode()).hexdigest()
                record = dict(metadata, artifacts=hashes)
                (directory / 'STEP_COMPLETE.json').write_text(json.dumps(record))
                return record
            def read(directory):
                return NS(), NS(), json.loads((directory / 'STEP_COMPLETE.json').read_text())
            with patch('cationize_protein.workflow._initial_structure', return_value=(NS(), NS())), \
                 patch('cationize_protein.workflow.inspect_structure', return_value=info), \
                 patch('cationize_protein.workflow._save_checkpoint', side_effect=save), \
                 patch('cationize_protein.workflow._read_checkpoint', side_effect=read), \
                 patch('cationize_protein.workflow._one_step', side_effect=WorkflowError('pending template')):
                with self.assertRaisesRegex(WorkflowError, 'pending'):
                    run(cfg)
                marker = cfg.output_dir / '0' / 'STAGE_COMPLETE.json'
                contents, mtime = marker.read_bytes(), marker.stat().st_mtime_ns
                with self.assertRaisesRegex(WorkflowError, 'pending'):
                    run(cfg)
                self.assertEqual(marker.read_bytes(), contents)
                self.assertEqual(marker.stat().st_mtime_ns, mtime)
            reference = cfg.output_dir / 'REFERENCE.json'
            reference.write_text(reference.read_text() + ' ')
            with self.assertRaisesRegex(WorkflowError, 'reference'):
                run(cfg)

    def test_resume_artifacts_are_verified_and_path_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data = root / 'protein.pdb'
            data.write_text('coordinate bytes')
            hashes = {'protein.pdb': hashlib.sha256(data.read_bytes()).hexdigest()}
            verify_artifacts(root, hashes)
            data.write_text('corrupted')
            with self.assertRaises(WorkflowError):
                verify_artifacts(root, hashes)
            with self.assertRaises(WorkflowError):
                verify_artifacts(root, {'../elsewhere': 'bad'})

    def test_initial_structure_loads_amber_gromacs_charmm_and_openmm_topologies(self):
        import parmed
        from cationize_protein.exports import export_system
        from cationize_protein.resources import default_parameters
        from cationize_protein.workflow import _initial_structure
        root_dir = Path(__file__).resolve().parents[1]
        ref_pdb = root_dir / 'examples/pdb_checks/reference_peptide.pdb'
        params = default_parameters()
        with tempfile.TemporaryDirectory() as temp:
            tmp = Path(temp)
            structure, _ = _pdb_structure(ref_pdb, params)
            exp_dir = tmp / 'exports'
            export_system(structure, exp_dir, formats=('amber', 'openmm', 'gromacs'), force_field='ff14SB')
            st_parm = parmed.load_file(str(exp_dir / 'system.parm7'), xyz=str(exp_dir / 'system.rst7'))
            st_parm.save(str(tmp / 'input.psf'))
            st_parm.save(str(tmp / 'input.crd'))
            for top, crd in [
                (exp_dir / 'system.parm7', exp_dir / 'system.rst7'),
                (exp_dir / 'system.top', exp_dir / 'system.gro'),
                (tmp / 'input.psf', tmp / 'input.crd'),
                (exp_dir / 'system.xml', exp_dir / 'state.json'),
                (exp_dir / 'system.xml', exp_dir / 'protein.pdb'),
            ]:
                cfg = WorkflowConfig(
                    output_dir=tmp / 'out', parameters_dir=params,
                    apbs=APBSConfig(executable=sys.executable), probe_radius_A=1.369,
                    topology_path=top, coordinates_path=crd, protonation='fixed',
                )
                loaded, sys_obj = _initial_structure(cfg)
                self.assertEqual(len(loaded.atoms), len(structure.atoms))
                self.assertEqual(sys_obj.getNumParticles(), len(structure.atoms))


if __name__ == '__main__':
    unittest.main()
