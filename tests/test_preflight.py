"""PDB source diagnostics must catch ambiguity before a parser discards it."""
from pathlib import Path
import tempfile
import unittest
from cationize_protein.preflight import check_pdb


def atom(serial, name, residue='ALA', number=1, chain='A', xyz=(0., 0., 0.), altloc='', element=None):
    element = element or name[0]
    return f'ATOM  {serial:5d} {name:>4s}{altloc:1s}{residue:>3s} {chain:1s}{number:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}'


def residue(number=1, chain='A', offset=0.):
    return [atom(i+1, n, number=number, chain=chain, xyz=(offset+x,y,z)) for i,(n,x,y,z) in enumerate([
        ('N',0.,0.,0.),('CA',1.45,0.,0.),('C',2.0,1.4,0.),('O',1.4,2.4,0.),('CB',1.8,-.8,1.2)])]


class PreflightTests(unittest.TestCase):
    def check(self, lines):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'source.pdb'; path.write_text('\n'.join(lines)+'\nEND\n')
            original=path.read_bytes(); output=Path(d)/'preflight.json'
            result=check_pdb(path, output_path=output)
            self.assertEqual(original,path.read_bytes())
            self.assertTrue(output.is_file())
            return result

    def test_complete_heavy_residue_can_be_prepared_without_hydrogens(self):
        result=self.check(residue())
        self.assertTrue(result['ok'])
        self.assertIn('missing_hydrogens', {i['code'] for i in result['issues']})
        self.assertIn('missing_terminal_oxygen', {i['code'] for i in result['issues']})
        self.assertEqual(result['residues'][0]['chain_id'], 'A')

    def test_backbone_and_sidechain_missing_atoms_are_explicit_blockers(self):
        result=self.check([a for a in residue() if a[12:16].strip() not in {'C','CB'}])
        self.assertFalse(result['ok'])
        missing=next(i for i in result['issues'] if i['code']=='missing_heavy_atoms')
        self.assertEqual(missing['missing_atoms'], ['C','CB'])
        self.assertIn('backbone',missing['action'].lower())

    def test_alternate_coordinates_duplicate_atoms_and_models_are_not_collapsed(self):
        lines=residue()+[atom(6,'CA',altloc='B')]+[residue()[0]]
        result=self.check(['MODEL        1']+lines+['ENDMDL','MODEL        2']+residue()+['ENDMDL'])
        codes={i['code'] for i in result['issues']}
        self.assertFalse(result['ok'])
        self.assertTrue({'multiple_models','alternate_locations','duplicate_atom'} <= codes)

    def test_multichain_repeated_numbers_are_valid_but_reused_identity_is_blocked(self):
        result=self.check(residue(9,'A')+['TER']+residue(9,'B',10.))
        self.assertTrue(result['ok'])
        self.assertEqual([r['chain_id'] for r in result['residues']], ['A','B'])
        result=self.check(residue(9,'A')+['TER']+residue(9,'A',10.))
        self.assertFalse(result['ok'])
        self.assertIn('duplicate_residue_identity',{i['code'] for i in result['issues']})

    def test_numbering_gap_and_broken_peptide_geometry_are_distinct(self):
        result=self.check(residue(1)+residue(4,offset=20.))
        codes={i['code'] for i in result['issues']}
        self.assertFalse(result['ok'])
        self.assertTrue({'residue_number_gap','peptide_chain_gap'} <= codes)

    def test_invalid_coordinates_and_unsupported_ligand_are_reported(self):
        line=atom(1,'N'); line=line[:30]+'     nan'+line[38:]
        result=self.check([line, atom(2,'FE',residue='HEM',number=2,element='FE')])
        codes={i['code'] for i in result['issues']}
        self.assertFalse(result['ok'])
        self.assertTrue({'invalid_coordinates','unsupported_residue'} <= codes)

    def test_declared_disulfide_with_missing_partner_is_an_explicit_error(self):
        declaration=list(' '*80);declaration[0:6]='SSBOND';declaration[15]='A';declaration[17:21]='   1';declaration[29]='A';declaration[31:35]='   7'
        result=self.check([''.join(declaration)]+[atom(1,'SG','CYS',1,element='S')])
        self.assertIn('invalid_disulfide_declaration',{i['code'] for i in result['issues']})
        self.assertEqual(len(result['ssbond_declarations']),1)

    def test_explicit_thiolate_or_disulfide_state_cannot_keep_a_thiol_hydrogen(self):
        result=self.check([atom(1,'SG','CYM',element='S'),atom(2,'HG','CYM',element='H')])
        self.assertIn('inconsistent_protonation_hydrogens',{i['code'] for i in result['issues']})

    def test_histidine_state_and_disulfide_are_reviewable(self):
        lines=[atom(1,'SG','CYS',1,xyz=(0.,0.,0.),element='S'),atom(2,'SG','CYS',2,xyz=(2.05,0.,0.),element='S')]
        result=self.check(lines)
        self.assertEqual(len(result['disulfide_candidates']),1)
        self.assertIn('disulfide_review',{i['code'] for i in result['issues']})
