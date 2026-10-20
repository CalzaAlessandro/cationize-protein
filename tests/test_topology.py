"""Topology-derived chains and independent carboxyl group identities."""
from types import SimpleNamespace as NS
from pathlib import Path
import tempfile
import unittest

import numpy as np

from cationize_protein.topology import inspect_structure, load_amber, UnsupportedTopology


LEGACY_NAMES = {
    'ASL': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'CB': 'C3', 'CG': 'C4', 'OD1': 'O2'},
    'ASC': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'CB': 'C3', 'CG': 'C4', 'OD1': 'O2'},
    'NSL': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'CB': 'C3', 'CG': 'C4', 'OD1': 'O2'},
    'ASX': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'CB': 'C3', 'CG': 'C4', 'OD1': 'O2'},
    'GLS': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'CB': 'C3', 'CG': 'C4', 'CD': 'C5', 'OE1': 'O2'},
    'GLC': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1', 'CB': 'C3', 'CG': 'C4', 'CD': 'C5', 'OE1': 'O2'},
    'LCL': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1'},
    'LET': {'N': 'N1', 'CA': 'C1', 'C': 'C2', 'O': 'O1'},
}


def structure(chains):
    """Minimal bonded ParmEd-compatible fixtures, without importing ParmEd."""
    atoms, residues, bonds = [], [], []
    for chain_number, sequence in enumerate(chains):
        last_c = None
        for local, name in enumerate(sequence):
            residue = NS(name=name, idx=len(residues), chain=chr(65 + chain_number),
                         number=42 + local, insertion_code='', atoms=[])
            names = ['N', 'CA', 'C', 'O', 'CB']
            links = [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB')]
            acidic = name in {'ASP', 'ASH', 'ASL', 'ASC', 'NSL', 'ASX', 'ASZ', 'AST', 'ASD'}
            glutamic = name in {'GLU', 'GLH', 'GLS', 'GLC', 'GLX', 'GLZ', 'GLT', 'GLD'}
            side_modified = name in {'ASL', 'ASC', 'NSL', 'ASX', 'ASZ', 'ASD', 'GLS', 'GLC', 'GLX', 'GLZ', 'GLD'}
            terminal_modified = name in {'LCL', 'LET', 'ALT', 'AST', 'ASD', 'GLT', 'GLD', 'CXT', 'CYT'}
            if acidic or glutamic:
                names += ['CG']
                links += [('CB', 'CG')]
                carbon, ox = 'CG', 'OD'
                if glutamic:
                    names += ['CD']
                    links += [('CG', 'CD')]
                    carbon, ox = 'CD', 'OE'
                names += [ox + '1']
                links += [(carbon, ox + '1')]
                if not side_modified:
                    names += [ox + '2']
                    links += [(carbon, ox + '2')]
            if name in {'CYS', 'CYX', 'CYM', 'CXT', 'CYT'}:
                names += ['SG']
                links += [('CB', 'SG')]
            if local == len(sequence) - 1 and not terminal_modified:
                names += ['OXT']
                links += [('C', 'OXT')]
            native = {}
            for native_name in names:
                atom_name = LEGACY_NAMES.get(name, {}).get(native_name, native_name)
                atomic_number = 8 if native_name.startswith('O') else 7 if native_name == 'N' else 16 if native_name == 'SG' else 6
                atom = NS(name=atom_name, idx=len(atoms), residue=residue, charge=.1,
                          solvent_radius=1.5, rmin=1.7, atomic_number=atomic_number)
                atoms.append(atom)
                residue.atoms.append(atom)
                native[native_name] = atom
            residues.append(residue)
            for a, b in links:
                bonds.append(NS(atom1=native[a], atom2=native[b]))
            if last_c is not None:
                bonds.append(NS(atom1=last_c, atom2=native['N']))
            last_c = native['C']
    return NS(atoms=atoms, residues=residues, bonds=bonds,
              coordinates=np.arange(len(atoms) * 3, dtype=float).reshape(-1, 3))


def atom_named(st, residue, name):
    return next(atom for atom in st.residues[residue].atoms if atom.name == name)


class TopologyTests(unittest.TestCase):
    def test_chain_termini_and_two_independent_acidic_terminal_groups(self):
        result = inspect_structure(structure([['ASP', 'GLU']]))
        self.assertEqual([r.context for r in result.residues], ['n_terminal', 'c_terminal'])
        self.assertEqual([(c.residue_ordinal, c.group, c.target_spec['name']) for c in result.candidates],
                         [(1, 'side', 'ASX'), (2, 'side', 'GLZ'), (2, 'terminal', 'GLT')])
        self.assertEqual(len({c.candidate_id for c in result.candidates}), 3)
        self.assertTrue(all(len(c.oxygen_indices) == 2 for c in result.candidates))

    def test_modified_terminal_states_select_only_remaining_group(self):
        for source, group, target in [('ASZ', 'terminal', 'ASD'), ('AST', 'side', 'ASD'),
                                      ('GLZ', 'terminal', 'GLD'), ('GLT', 'side', 'GLD')]:
            with self.subTest(source=source):
                result = inspect_structure(structure([['ALA', source]]))
                self.assertEqual([(c.group, c.target_spec['name']) for c in result.candidates], [(group, target)])
        for source in ['ASD', 'GLD']:
            with self.subTest(source=source):
                self.assertEqual(inspect_structure(structure([['ALA', source]])).candidates, ())

    def test_aliases_resolve_existing_backbone_names_without_duplicate_side(self):
        for source in ['ASL', 'ASC', 'GLS', 'GLC']:
            with self.subTest(source=source):
                result = inspect_structure(structure([['ALA', source, 'LEU']]))
                middle = result.residues[1]
                self.assertEqual(middle.native_atom_names['N'], 'N1')
                self.assertEqual(middle.native_atom_names['CA'], 'C1')
                self.assertFalse(any(c.residue_ordinal == 2 for c in result.candidates))
                self.assertIn(middle.native_atom_indices['CA'], result.protected_atoms)
        result = inspect_structure(structure([['NSL', 'LCL']]))
        self.assertEqual([r.canonical_name for r in result.residues], ['ASX', 'LET'])
        self.assertEqual(result.candidates, ())

    def test_chain_and_original_residue_identity_survive_repeated_numbers(self):
        result = inspect_structure(structure([['ASP', 'LEU'], ['ASP', 'LEU']]))
        self.assertEqual(len(result.chains), 2)
        self.assertEqual(len({c.candidate_id for c in result.candidates}), 4)
        self.assertEqual(result.residues[0].original_id, '42')
        self.assertEqual(result.residues[2].chain_id, 'B')
        self.assertIn('A', result.candidates[0].candidate_id)

    def test_disulfides_do_not_merge_peptide_chains_and_drive_cyx_state(self):
        st = structure([['CYS', 'LEU'], ['CYS', 'LEU']])
        a, b = atom_named(st, 0, 'SG'), atom_named(st, 2, 'SG')
        st.bonds.append(NS(atom1=a, atom2=b))
        result = inspect_structure(st)
        self.assertEqual(len(result.chains), 2)
        self.assertEqual(result.disulfide_bonds, ((a.idx, b.idx),))
        self.assertEqual(result.residues[0].parent, 'CYX')
        self.assertIn(a.idx, result.protected_atoms)
        self.assertIn(b.idx, result.protected_atoms)

    def test_disulfide_c_terminal_candidate_preserves_cyx_template(self):
        st = structure([['ALA', 'CYS'], ['ALA', 'CYS']])
        a, b = atom_named(st, 1, 'SG'), atom_named(st, 3, 'SG')
        st.bonds.append(NS(atom1=a, atom2=b))
        result = inspect_structure(st)
        self.assertEqual([c.target_spec['name'] for c in result.candidates], ['CXT', 'CXT'])

    def test_unsupported_protonation_and_nonprotein_residues_are_rejected(self):
        for name in ['ASH', 'GLH', 'HOH', 'HIS']:
            with self.subTest(name=name), self.assertRaises(UnsupportedTopology):
                inspect_structure(structure([[name, 'LEU']]))

    def test_modified_context_mismatch_and_unbonded_cyx_are_rejected(self):
        for sequence in [['ASL', 'LEU'], ['ALA', 'ASX'], ['CYX', 'LEU']]:
            with self.subTest(sequence=sequence), self.assertRaises(UnsupportedTopology):
                inspect_structure(structure([sequence]))

    def test_single_residue_chain_and_unknown_crosslink_are_rejected(self):
        with self.assertRaises(UnsupportedTopology):
            inspect_structure(structure([['ASP']]))
        st = structure([['ALA', 'LEU']])
        st.bonds.append(NS(atom1=atom_named(st, 0, 'CB'), atom2=atom_named(st, 1, 'CB')))
        with self.assertRaises(UnsupportedTopology):
            inspect_structure(st)

    def test_branching_peptide_and_cyclic_peptide_are_rejected(self):
        for kind in ['branch', 'cycle']:
            st = structure([['ALA', 'LEU'], ['ALA', 'LEU']])
            if kind == 'branch':
                st.bonds.append(NS(atom1=atom_named(st, 0, 'C'), atom2=atom_named(st, 2, 'N')))
            else:
                st.bonds.append(NS(atom1=atom_named(st, 1, 'C'), atom2=atom_named(st, 0, 'N')))
            with self.subTest(kind=kind), self.assertRaises(UnsupportedTopology):
                inspect_structure(st)

    def test_missing_carboxyl_oxygen_bond_and_nonfinite_input_are_rejected(self):
        for defect in ['oxygen_bond', 'charge', 'radius', 'coordinates']:
            st = structure([['ASP', 'LEU']])
            if defect == 'oxygen_bond':
                atom = atom_named(st, 0, 'OD2')
                st.bonds = [bond for bond in st.bonds if bond.atom1 is not atom and bond.atom2 is not atom]
            elif defect == 'charge':
                st.atoms[0].charge = float('nan')
            elif defect == 'radius':
                st.atoms[0].solvent_radius = -1
            else:
                st.coordinates[0, 0] = float('inf')
            with self.subTest(defect=defect), self.assertRaises(UnsupportedTopology):
                inspect_structure(st)

    def test_protonated_or_misidentified_carboxyl_oxygen_is_rejected(self):
        for defect in ['hydrogen', 'element']:
            st = structure([['ASP', 'LEU']])
            oxygen = atom_named(st, 0, 'OD2')
            if defect == 'element':
                oxygen.atomic_number = 7
            else:
                hydrogen = NS(name='HD2', idx=len(st.atoms), residue=st.residues[0],
                              charge=.1, solvent_radius=1.2, rmin=0., atomic_number=1)
                st.atoms.append(hydrogen)
                st.residues[0].atoms.append(hydrogen)
                st.bonds.append(NS(atom1=oxygen, atom2=hydrogen))
                st.coordinates = np.vstack([st.coordinates, [[0, 0, 0]]])
            with self.subTest(defect=defect), self.assertRaises(UnsupportedTopology):
                inspect_structure(st)

    def test_zero_pb_radii_are_not_accepted_as_parameterized_protein(self):
        st = structure([['ALA', 'LEU']])
        st.atoms[0].solvent_radius = 0
        with self.assertRaises(UnsupportedTopology):
            inspect_structure(st)

    def test_non_sulfur_disulfide_endpoint_is_rejected(self):
        st = structure([['CYS', 'LEU'], ['CYS', 'LEU']])
        a, b = atom_named(st, 0, 'SG'), atom_named(st, 2, 'SG')
        a.atomic_number = 6
        st.bonds.append(NS(atom1=a, atom2=b))
        with self.assertRaises(UnsupportedTopology):
            inspect_structure(st)

    def test_inspection_copies_data_and_preserves_original_input_identity(self):
        st = structure([['CYS', 'LEU'], ['CYS', 'LEU']])
        a, b = atom_named(st, 0, 'SG'), atom_named(st, 2, 'SG')
        st.bonds.append(NS(atom1=a, atom2=b))
        st.residues[0].insertion_code = 'A'
        result = inspect_structure(st, provenance={'source': 'test'})
        self.assertEqual(result.residues[0].original_id, '42A')
        self.assertEqual(st.residues[0].name, 'CYS')
        result.coordinates[0, 0] = -999
        self.assertEqual(st.coordinates[0, 0], 0)
        self.assertEqual(result.provenance, {'source': 'test'})

    def test_load_checks_input_paths_before_optional_dependency(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileNotFoundError):
                load_amber(Path(temp) / 'missing.parm7', Path(temp) / 'missing.rst7')


if __name__ == '__main__':
    unittest.main()
