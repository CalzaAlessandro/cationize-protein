"""UniProt bonds must map by sequence and survive preparation as real bonds."""
import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from cationize_protein.preflight import check_pdb
from cationize_protein.preparation import prepare_pdb
from test_preflight import atom


def record(sequence='MCCK', pairs=((2, 3),)):
    return dict(primaryAccession='TEST01', sequence={'value': sequence}, features=[
        dict(type='Disulfide bond', location={
            'start': {'value': a, 'modifier': 'EXACT'},
            'end': {'value': b, 'modifier': 'EXACT'}},
            description='', evidences=[{'evidenceCode': 'ECO:0000269',
                                        'source': 'PubMed', 'id': '123'}]) for a, b in pairs])


def cysteine_pair(*, hydrogen=False, insertion=False):
    coordinates=[{'N':(0.,0.,0.),'CA':(1.45,0.,0.),'C':(2.6,.8,0.),
                  'O':(3.7,.4,0.),'CB':(1.45,-1.5,0.),'SG':(1.45,-3.3,0.)},
                 {'N':(3.8,1.4,0.),'CA':(5.2,1.4,0.),'C':(6.2,2.4,0.),
                  'O':(7.4,2.2,0.),'CB':(4.3,-1.,0.),'SG':(3.5,-3.3,0.),
                  'OXT':(6.,3.6,0.)}]
    if hydrogen:
        coordinates[0]['HG']=(1.45,-4.6,0.)
    rows=[]
    for index, atoms in enumerate(coordinates):
        for name, xyz in atoms.items():
            row=atom(len(rows)+1,name,'CYS',105 if insertion else 105+index,'A',xyz,
                     element='S' if name=='SG' else name[0])
            if insertion:row=row[:26]+('A' if index==0 else 'B')+row[27:]
            rows.append(row)
    return '\n'.join(rows)+'\nEND\n'


class UniProtPreparationTests(unittest.TestCase):
    def inputs(self, root, data=None, pdb=None):
        source=root/'input.pdb'; source.write_text(pdb or cysteine_pair())
        reference=root/'uniprot.json'; reference.write_text(json.dumps(data or record()))
        return source, reference, root/'prepared'

    def test_sequence_offset_is_mapped_and_the_output_contains_a_real_bond(self):
        from openmm.app import PDBFile
        with tempfile.TemporaryDirectory() as tmp:
            source, reference, output=self.inputs(Path(tmp),pdb=cysteine_pair(hydrogen=True))
            original=source.read_bytes()
            result=prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertEqual(source.read_bytes(),original)
            self.assertEqual((output/'original.pdb').read_bytes(),original)
            self.assertEqual((output/'uniprot.json').read_bytes(),reference.read_bytes())
            checked=check_pdb(output/'protein.pdb')
            self.assertTrue(checked['ok'],checked['messages'])
            self.assertEqual([r['residue_name'] for r in checked['residues']],['CYX','CYX'])
            self.assertFalse(any(a['name']=='HG' for r in checked['residues'] for a in r['atoms']))
            self.assertEqual(checked['ssbond_declarations'][0]['residue_ordinals'],[1,2])
            pdb=PDBFile(str(output/'protein.pdb'))
            self.assertEqual([tuple(sorted((a.residue.id,b.residue.id))) for a,b in pdb.topology.bonds()
                              if a.name==b.name=='SG'],[('105','106')])
            with (output/'disulfides.csv').open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual([(r['uniprot_start'],r['uniprot_end'],r['pdb_residue_1'],
                               r['pdb_residue_2'],r['status']) for r in rows],
                             [('2','3','105','106','verified')])
            self.assertIn('UniProt', '\n'.join(result['messages']))
            heavy=lambda report:[(r['ordinal'],a['name'],a['coordinates_A'])
                for r in report['residues'] for a in r['atoms'] if not a['hydrogen']]
            self.assertEqual(heavy(check_pdb(source)),heavy(checked))

    def test_insertion_codes_are_preserved_in_disulfide_declarations(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),pdb=cysteine_pair(insertion=True))
            prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            checked=check_pdb(output/'protein.pdb')
            self.assertTrue(checked['ok'],checked['messages'])
            self.assertEqual([r['original_id'] for r in checked['residues']],['105A','105B'])
            self.assertEqual([p['original_id'] for p in checked['ssbond_declarations'][0]['endpoints']],
                             ['105A','105B'])

    def test_ambiguous_sequence_mapping_is_rejected_without_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),record('MCCKCCN'))
            with self.assertRaisesRegex(ValueError,'ambiguous'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

    def test_a_different_sequence_is_not_mapped_using_residue_numbers(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),record('MCDCK',((2,4),)))
            with self.assertRaisesRegex(ValueError,'sequence|match'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

    def test_distorted_sulfur_geometry_is_not_repaired_by_renaming(self):
        rows=cysteine_pair().splitlines()
        rows=[line[:30]+f'{8.:8.3f}'+line[38:]
              if line.startswith('ATOM') and line[12:16].strip()=='SG' and line[22:26].strip()=='106'
              else line for line in rows]
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),pdb='\n'.join(rows)+'\n')
            with self.assertRaisesRegex(ValueError,'distance|geometry'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

    def test_missing_partner_within_a_truncated_construct_requires_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),record('CCMC',((2,4),)))
            with self.assertRaisesRegex(ValueError,'partner|partial|missing'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

    def test_bonds_entirely_outside_the_construct_are_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),record('CACMCC',((1,3),)))
            prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            with (output/'disulfides.csv').open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual([r['status'] for r in rows],['outside_construct'])

    def test_uncertain_or_interchain_annotations_are_not_applied_as_intrachain_bonds(self):
        for kind in ('uncertain','interchain'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as tmp:
                data=record()
                if kind=='uncertain':data['features'][0]['location']['start']['modifier']='UNSURE'
                else:data['features'][0]['description']='Interchain'
                source,reference,output=self.inputs(Path(tmp),data)
                with self.assertRaisesRegex(ValueError,'exact|interchain|Interchain|uncertain'):
                    prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
                self.assertFalse(output.exists())

    def test_a_sulfur_cannot_acquire_two_annotated_partners(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),record('MCCCK',((2,3),(3,4))))
            with self.assertRaisesRegex(ValueError,'partner|conflict|multiple'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

    def test_cli_accepts_the_local_record_and_reports_the_mapping(self):
        from cationize_protein.cli import main
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()) as stream:
            source,reference,output=self.inputs(Path(tmp))
            code=main(['prepare',str(source),'--output',str(output),'--protonation','fixed',
                       '--uniprot-record',str(reference)])
            self.assertEqual(code,0)
            self.assertIn('UniProt',stream.getvalue())
            self.assertIn('105',stream.getvalue())

    def test_a_null_json_record_cannot_silently_disable_requested_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp))
            reference.write_text('null\n')
            with self.assertRaisesRegex(ValueError,'UniProt|JSON'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

    def test_alternative_or_conditional_annotation_text_requires_review(self):
        for description in ('or C-3 with C-5', 'In oxidized form', 'Between processed chains A and B'):
            with self.subTest(description=description),tempfile.TemporaryDirectory() as tmp:
                data=record('MCCKC');data['features'][0]['description']=description
                source,reference,output=self.inputs(Path(tmp),data,pdb=cysteine_pair(hydrogen=True))
                with self.assertRaisesRegex(ValueError,'qualified|Qualified|review'):
                    prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
                self.assertFalse(output.exists())

    def test_malformed_annotation_fields_produce_a_clear_validation_error(self):
        for field in ('features','feature','location','description'):
            with self.subTest(field=field),tempfile.TemporaryDirectory() as tmp:
                data=record()
                if field=='features':data['features']=None
                elif field=='feature':data['features']=[None]
                else:data['features'][0][field]=None
                source,reference,output=self.inputs(Path(tmp),data)
                with self.assertRaisesRegex(ValueError,'UniProt|annotation'):
                    prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
                self.assertFalse(output.exists())

    def test_duplicate_atom_serials_cannot_delete_an_unrelated_hydrogen(self):
        rows=cysteine_pair(hydrogen=True).splitlines()
        rows.insert(-1,atom(7,'H','CYS',106,'A',(3.8,1.4,1.),element='H'))
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),pdb='\n'.join(rows)+'\n')
            original=source.read_bytes()
            with self.assertRaisesRegex(ValueError,'unique.*serial|duplicate.*serial'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertEqual(source.read_bytes(),original)
            self.assertFalse(output.exists())

    def test_noncontiguous_blocks_with_the_same_chain_id_are_not_joined(self):
        original=[line for line in cysteine_pair().splitlines() if line.startswith('ATOM')]
        first=[line for line in original if line[22:26].strip()=='105']
        first.append(atom(99,'OXT','CYS',105,'A',(2.4,1.9,0.)))
        coords={'N':(10.,0.,0.),'CA':(11.4,0.,0.),'C':(12.6,.8,0.),'O':(13.7,.4,0.),
                'OXT':(12.4,1.9,0.),'CB':(11.4,-1.5,0.),'CG':(11.4,-3.,0.),
                'CD':(11.4,-4.5,0.),'CE':(11.4,-6.,0.),'NZ':(11.4,-7.5,0.)}
        middle=[atom(99,name,'LYS',1,'B',xyz) for name,xyz in coords.items()]
        last=[line for line in original if line[22:26].strip()=='106']
        rows=[line[:6]+f'{i:5d}'+line[11:] for i,line in enumerate(first+middle+last,1)]
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),pdb='\n'.join(rows)+'\nEND\n')
            with self.assertRaisesRegex(ValueError,'chain|ambiguous'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

    def test_existing_matching_ssbond_is_preserved_without_duplication(self):
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp))
            prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            second=Path(tmp)/'second'
            prepare_pdb(output/'protein.pdb',second,protonation='fixed',uniprot_record=reference)
            self.assertEqual((second/'protein.pdb').read_bytes(),(output/'protein.pdb').read_bytes())
            self.assertEqual(len(check_pdb(second/'protein.pdb')['ssbond_declarations']),1)

    def test_homomer_chains_receive_separate_intrachain_pairs(self):
        from openmm.app import PDBFile
        rows=[line for line in cysteine_pair().splitlines() if line.startswith('ATOM')]
        second=[]
        for line in rows:
            second.append(line[:6]+f'{int(line[6:11])+20:5d}'+line[11:21]+'B'+line[22:46]+f'{10.:8.3f}'+line[54:])
        pdb='\n'.join(rows+['TER']+second+['TER','END'])+'\n'
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),pdb=pdb)
            prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            parsed=PDBFile(str(output/'protein.pdb'))
            pairs=[{a.residue.chain.id,b.residue.chain.id} for a,b in parsed.topology.bonds()
                   if a.name==b.name=='SG']
            self.assertEqual(pairs,[{'A'},{'B'}])

    def test_existing_ssbond_with_a_different_partner_is_not_replaced(self):
        rows=[line for line in cysteine_pair().splitlines()
              if line.startswith('ATOM') and line[12:16].strip()!='OXT']
        coords={'N':(7.3,3.1,0.),'CA':(8.6,3.1,0.),'C':(9.8,3.9,0.),
                'O':(10.8,3.5,0.),'OXT':(9.6,5.,0.),'CB':(8.6,1.6,0.),'SG':(1.45,-3.3,2.05)}
        for name,xyz in coords.items():
            rows.append(atom(len(rows)+1,name,'CYS',107,'A',xyz,element='S' if name=='SG' else name[0]))
        declaration=list(' '*80);declaration[:6]='SSBOND';declaration[7:10]='  1'
        for start,number in ((11,105),(25,107)):
            declaration[start:start+3]='CYS';declaration[start+4]='A';declaration[start+6:start+10]=f'{number:4d}'
        pdb=''.join(declaration)+'\n'+'\n'.join(rows)+'\nEND\n'
        with tempfile.TemporaryDirectory() as tmp:
            source,reference,output=self.inputs(Path(tmp),record('MCCCK'),pdb=pdb)
            with self.assertRaisesRegex(ValueError,'SSBOND.*conflict'):
                prepare_pdb(source,output,protonation='fixed',uniprot_record=reference)
            self.assertFalse(output.exists())

