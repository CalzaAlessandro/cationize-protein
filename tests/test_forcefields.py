import tempfile,unittest,xml.etree.ElementTree as E
from pathlib import Path
from cationize_protein.forcefields import profile,custom_xml
OFF='''!!index array str
 "ALT"
!entry.ALT.unit.atoms table str name str type int typex int resx int flags int seq int elmnt dbl chg
 "N" "N" 0 1 0 1 7 -0.4
 "CA" "CX" 0 1 0 2 6 0.2
 "C" "C" 0 1 0 3 6 0.6
 "O" "O" 0 1 0 4 8 -0.4
!entry.ALT.unit.connectivity table int atom1x int atom2x int flags
 1 2 1
 2 3 1
 3 4 1
!entry.ALT.unit.positions table dbl x dbl y dbl z
 0 0 0
 1 0 0
 1 1 0
 1 2 0
'''
class Fields(unittest.TestCase):
 def fake(self,p,ff19=True):
  types=[('n','N','N',14),('ca','CT','C',12),('c','C','C',12),('o','O','O',16)]
  r=E.Element('ForceField');ats=E.SubElement(r,'AtomTypes')
  for name,cls,element,mass in types:E.SubElement(ats,'Type',name=name,**{'class':cls},element=element,mass=str(mass))
  residues=E.SubElement(r,'Residues');res=E.SubElement(residues,'Residue',name='CALA')
  for an,t in [('N','n'),('CA','ca'),('C','c'),('O','o')]:E.SubElement(res,'Atom',name=an,type=t)
  nb=E.SubElement(r,'NonbondedForce',coulomb14scale='.833333',lj14scale='.5')
  if ff19:E.SubElement(nb,'UseAttributeFromResidue',name='charge')
  for name,cls,element,mass in types:E.SubElement(nb,'Atom',type=name,sigma='.3',epsilon='.2',**({} if ff19 else {'charge':'.123'}))
  (p/'parameters').mkdir();(p/'parameters/ALT.lib').write_text(OFF)
  path=p/profile('ff19SB' if ff19 else 'ff99SBildn')['xml'];path.parent.mkdir(parents=True,exist_ok=True);E.ElementTree(r).write(path)
 def test_ff19_product_keeps_native_backbone_type_identity_for_cmap(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);self.fake(p);r=E.fromstring(custom_xml(p/'parameters','ff19SB',['ALT'],data_dir=p))
   atoms=r.find('Residues/Residue').findall('Atom');self.assertEqual([a.get('type') for a in atoms],['n','ca','c','o'])
   self.assertAlmostEqual(sum(float(a.get('charge')) for a in atoms),0)
 def test_legacy_field_receives_unique_charge_types(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);self.fake(p,False);r=E.fromstring(custom_xml(p/'parameters','ff99SBildn',['ALT'],data_dir=p))
   nb=r.find('NonbondedForce');self.assertEqual([float(a.get('charge')) for a in nb],[-.4,.2,.6,-.4])
   self.assertEqual(len({a.get('type') for a in nb}),4)
 def test_unknown_field_is_not_silently_replaced(self):
  with self.assertRaises(ValueError):profile('openff')
