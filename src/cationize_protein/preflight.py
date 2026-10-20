"""Read-only PDB adequacy checks before OpenMM/ParmEd can discard ambiguity.

This checks supplied coordinates and heavy-atom completeness, not sequence
coverage or biological identity. No missing backbone/loop is reconstructed.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path

from .catalog import templates

_SIDE = {
    'ALA':'CB', 'ARG':'CB CG CD NE CZ NH1 NH2', 'ASN':'CB CG OD1 ND2',
    'ASP':'CB CG OD1 OD2', 'CYS':'CB SG', 'GLN':'CB CG CD OE1 NE2',
    'GLU':'CB CG CD OE1 OE2', 'GLY':'', 'HIS':'CB CG ND1 CD2 CE1 NE2',
    'ILE':'CB CG1 CG2 CD1', 'LEU':'CB CG CD1 CD2', 'LYS':'CB CG CD CE NZ',
    'MET':'CB CG SD CE', 'PHE':'CB CG CD1 CD2 CE1 CE2 CZ',
    'PRO':'CB CG CD', 'SER':'CB OG', 'THR':'CB OG1 CG2',
    'TRP':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR':'CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL':'CB CG1 CG2',
}
_NATIVE_PARENT = dict(HID='HIS',HIE='HIS',HIP='HIS',ASH='ASP',GLH='GLU',
                      CYX='CYS',CYM='CYS',LYN='LYS',TYM='TYR')
_ALIASES = dict(ASL='ASC',GLS='GLC',NSL='ASX',LCL='LET')


def _heavy_names(name):
    canonical = _ALIASES.get(name,name)
    spec = next((s for s in templates() if s['name']==canonical), None)
    if spec:
        # OFF metadata is authoritative for reused noncanonical atom names.
        from .resources import default_parameters
        from .off import read_units
        path=default_parameters()/f'{canonical}.lib'
        if not path.is_file():
            return None
        return {a['name'] for a in read_units(path)[canonical]['atoms'] if a['atomic_number']!=1}
    parent=_NATIVE_PARENT.get(name,name)
    if parent not in _SIDE:
        return None
    return set(('N CA C O '+_SIDE[parent]).split())


def check_pdb(path, output_path=None):
    """Return an English, JSON-compatible report; optionally write a new JSON.

    ``ok`` is false when atom loss, unresolved identity, unsupported chemistry,
    malformed coordinates, or an internal peptide break makes preparation unsafe.
    Missing H/OXT are preparation warnings. Numbering gaps are distinguished from
    geometric breaks; numbering alone cannot establish a missing sequence segment.
    Disulfide distances are candidates for review, never an automatic bond edit.
    Residue ordinals include every ATOM/HETATM residue in source order.
    """
    path=Path(path)
    raw=path.read_bytes()
    lines=raw.decode('utf-8').splitlines()
    issues=[]; residues=[]; current=None; segment=0; model=0; model_tags=[]
    identities={}; atom_count=0; hydrogen_count=0; alternate=[]; declarations=[]

    def issue(code, severity, message, action, **details):
        issues.append(dict(code=code,severity=severity,message=message,action=action,**details))

    for line_number,line in enumerate(lines,1):
        record=line[:6].strip()
        if record=='SSBOND':
            if len(line)<36:
                issue('invalid_disulfide_declaration','error',f'Line {line_number}: truncated SSBOND record.',
                      'Restore the complete disulfide partner identifiers.',line=line_number)
            else:
                declarations.append(dict(line=line_number,endpoints=[
                    dict(chain_id=line[15].strip(),original_id=line[17:21].strip()+line[21].strip()),
                    dict(chain_id=line[29].strip(),original_id=line[31:35].strip()+line[35].strip())]))
        if record=='MODEL':
            model+=1; model_tags.append(line[6:].strip()); current=None; segment+=1
        if record in {'ENDMDL','TER'}:
            segment+=1; current=None
        if record not in {'ATOM','HETATM'}:
            continue
        atom_count+=1
        if len(line)<54:
            issue('malformed_atom_record','error',f'Line {line_number}: truncated atom record.',
                  'Provide fixed-column PDB atom records with complete coordinates.',line=line_number)
            continue
        name=line[12:16].strip(); alt=line[16].strip(); rn=line[17:20].strip()
        chain=line[21].strip(); number=line[22:26].strip(); insertion=line[26].strip()
        try:
            int(number)
        except ValueError:
            issue('invalid_residue_number','error',f'Line {line_number}: invalid residue number {number!r}.',
                  'Correct the PDB residue numbering without changing atom identity.',line=line_number)
        key=(model,segment,chain,number,insertion)
        identity=(model,chain,number,insertion)
        if current is None or current['_key']!=key:
            ordinal=len(residues)+1
            if identity in identities:
                issue('duplicate_residue_identity','error',f'Residue {chain}:{number}{insertion} occurs in separate blocks.',
                      'Assign unique chain and residue identities or correct the residue ordering.',residue_ordinal=ordinal)
            identities[identity]=ordinal
            current=dict(ordinal=ordinal,residue_name=rn,chain_id=chain,original_id=number+insertion,
                         model=model or 1,segment=segment,atoms=[],_key=key)
            residues.append(current)
        elif current['residue_name']!=rn:
            issue('conflicting_residue_names','error',f'Residue ordinal {current["ordinal"]} contains multiple residue names.',
                  'Resolve alternate residue identities explicitly.',residue_ordinal=current['ordinal'])
        try:
            xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
            if not all(math.isfinite(x) for x in xyz):
                raise ValueError
        except ValueError:
            xyz=None
            issue('invalid_coordinates','error',f'Line {line_number}: coordinates are malformed or non-finite.',
                  'Supply finite x, y and z coordinates for this atom.',line=line_number,residue_ordinal=current['ordinal'])
        element=line[76:78].strip().upper() if len(line)>=78 else ''
        is_h=element in {'H','D'} or (not element and name.lstrip('0123456789').startswith(('H','D')))
        hydrogen_count+=is_h
        if not name:
            issue('empty_atom_name','error',f'Line {line_number}: empty atom name.',
                  'Restore an unambiguous PDB atom name.',line=line_number)
        current['atoms'].append(dict(name=name,altloc=alt,element=element,hydrogen=is_h,coordinates_A=xyz,line=line_number))
        if alt:
            alternate.append(line_number)
    if len(model_tags)>1:
        issue('multiple_models','error',f'The PDB contains {len(model_tags)} coordinate models.',
              'Select one explicit MODEL and save it as a separate input PDB.')
    if alternate:
        issue('alternate_locations','error',f'{len(alternate)} atom records have alternate locations.',
              'Choose one coherent residue conformer explicitly; do not combine alternate atoms.',lines=alternate)
    if not atom_count:
        issue('empty_structure','error','No ATOM/HETATM coordinates were found.',
              'Provide a protein PDB with coordinates.')
    if not hydrogen_count:
        issue('missing_hydrogens','warning','No hydrogen coordinates are present.',
              'Add hydrogens during preparation with an explicit pH and residue-state policy.')

    for r in residues:
        names=[a['name'] for a in r['atoms']]
        duplicates=sorted(n for n,count in Counter(names).items() if count>1)
        if duplicates:
            issue('duplicate_atom','error',f'Residue ordinal {r["ordinal"]} has repeated atom names: {", ".join(duplicates)}.',
                  'Resolve duplicate/alternate atom records before preparation.',residue_ordinal=r['ordinal'],atoms=duplicates)
        expected=_heavy_names(r['residue_name'])
        if expected is None:
            issue('unsupported_residue','error',f'Residue ordinal {r["ordinal"]} ({r["residue_name"]}) is unsupported.',
                  'Provide a protein-only input or explicitly parameterize and authorize this ligand, solvent, ion or cap.',residue_ordinal=r['ordinal'])
            continue
        present={a['name'] for a in r['atoms'] if not a['hydrogen']}
        missing=sorted(expected-present)
        if missing:
            backbone=bool({'N','CA','C','O','N1','C1','C2','O1'} & set(missing))
            issue('missing_heavy_atoms','error',f'Residue ordinal {r["ordinal"]} is missing heavy atoms: {", ".join(missing)}.',
                  ('Supply a complete backbone/loop model from an explicitly reviewed reconstruction; no automatic backbone rebuilding is performed.'
                   if backbone else 'Reconstruct and validate the missing side chain in a separate input before starting.'),
                  residue_ordinal=r['ordinal'],missing_atoms=missing)
        extra=sorted(present-expected-{'OXT'})
        if extra:
            issue('unexpected_heavy_atoms','error',f'Residue ordinal {r["ordinal"]} has unexpected heavy atoms: {", ".join(extra)}.',
                  'Check atom naming and the residue chemical state.',residue_ordinal=r['ordinal'],atoms=extra)
        if r['residue_name'] in {'CYX','CYM'} and any(a['name']=='HG' for a in r['atoms']):
            issue('inconsistent_protonation_hydrogens','error',f'Residue ordinal {r["ordinal"]} is {r["residue_name"]} but contains a thiol HG hydrogen.',
                  'Resolve the sulfur oxidation/protonation state before hydrogen addition or bond assignment.',residue_ordinal=r['ordinal'])
        if r['residue_name']=='HIS':
            issue('histidine_tautomer_unspecified','warning',f'Histidine ordinal {r["ordinal"]} has no explicit HID/HIE/HIP label.',
                  'Choose and record a histidine tautomer during hydrogen addition; PROPKA does not choose neutral tautomers.',residue_ordinal=r['ordinal'])
        elif r['residue_name'] in _NATIVE_PARENT:
            issue('explicit_protonation_state','info',f'Residue ordinal {r["ordinal"]} supplies state {r["residue_name"]}.',
                  'Retain this state until a supported, persistent pH audit justifies a change.',residue_ordinal=r['ordinal'])

    def atom_xyz(r, name):
        aliases={'N':'N1','C':'C2'} if _ALIASES.get(r['residue_name'],r['residue_name']) in {'ASC','GLC','ASX','LET'} else {}
        return next((a['coordinates_A'] for a in r['atoms'] if a['name']==aliases.get(name,name)), None)

    for i,r in enumerate(residues):
        after=residues[i+1] if i+1<len(residues) else None
        same_segment=after and r['segment']==after['segment'] and r['model']==after['model'] and r['chain_id']==after['chain_id']
        if same_segment:
            try:
                first=int(r['original_id']); second=int(after['original_id'])
                if second-first!=1:
                    issue('residue_number_gap','warning',f'Residue numbering jumps from {r["original_id"]} to {after["original_id"]} in chain {r["chain_id"] or "(blank)"}.',
                          'Compare with the intended sequence; numbering gaps alone do not prove missing coordinates.',residue_ordinals=[r['ordinal'],after['ordinal']])
            except ValueError:
                pass
            c,n=atom_xyz(r,'C'),atom_xyz(after,'N')
            if c is not None and n is not None:
                distance=math.dist(c,n)
                if not 1.0<=distance<=1.8:
                    issue('peptide_chain_gap','error',f'C–N distance between residue ordinals {r["ordinal"]} and {after["ordinal"]} is {distance:.3f} Å.',
                          'Review the missing segment or bond geometry; declare true chain boundaries with TER and do not join a gap automatically.',
                          residue_ordinals=[r['ordinal'],after['ordinal']],distance_A=distance)
        elif not any(a['name']=='OXT' for a in r['atoms']):
            spec=next((s for s in templates() if s['name']==_ALIASES.get(r['residue_name'],r['residue_name'])),None)
            if not spec or 'terminal' not in spec['modified_groups']:
                issue('missing_terminal_oxygen','warning',f'Chain-end residue ordinal {r["ordinal"]} has no OXT.',
                      'Add only the terminal oxygen during explicit preparation after confirming this is the intended chain end.',residue_ordinal=r['ordinal'])

    sulfurs=[(r,a) for r in residues for a in r['atoms'] if a['name']=='SG' and a['coordinates_A'] is not None]
    candidates=[]
    for i,(r,a) in enumerate(sulfurs):
        for s,b in sulfurs[i+1:]:
            if r['model']!=s['model']: continue
            distance=math.dist(a['coordinates_A'],b['coordinates_A'])
            if distance<=2.5:
                candidates.append(dict(residue_ordinals=[r['ordinal'],s['ordinal']],distance_A=distance))
    if candidates:
        issue('disulfide_review','warning',f'{len(candidates)} SG–SG proximity candidates require bond verification.',
              'Confirm disulfides from the intended structure and graph; distances do not assign oxidation states.')
    for declaration in declarations:
        partners=[]
        for endpoint in declaration['endpoints']:
            matches=[(r,a) for r,a in sulfurs if r['chain_id']==endpoint['chain_id'] and r['original_id']==endpoint['original_id']]
            if len(matches)!=1:break
            partners.append(matches[0])
        if len(partners)!=2 or partners[0][0]['ordinal']==partners[1][0]['ordinal']:
            issue('invalid_disulfide_declaration','error',f'SSBOND line {declaration["line"]} does not identify two unique supplied sulfur atoms.',
                  'Supply both declared disulfide partners and unique residue identities.',line=declaration['line'])
        else:
            distance=math.dist(partners[0][1]['coordinates_A'],partners[1][1]['coordinates_A'])
            declaration['residue_ordinals']=[p[0]['ordinal'] for p in partners];declaration['distance_A']=distance
            if not 1.7<=distance<=2.5:
                issue('disulfide_geometry','error',f'Declared SSBOND line {declaration["line"]} has SG–SG distance {distance:.3f} Å.',
                      'Review and correct the declared bond/coordinates before preparation.',line=declaration['line'],distance_A=distance)
    for r,a in sulfurs:
        if r['residue_name']=='CYX' and not any(r['ordinal'] in c['residue_ordinals'] for c in candidates):
            issue('unpaired_cyx','error',f'CYX ordinal {r["ordinal"]} has no nearby disulfide partner.',
                  'Supply the correct disulfide partner or correct the chemical state.',residue_ordinal=r['ordinal'])
    for r in residues:
        del r['_key']
    report=dict(schema_version=1,source=str(path.resolve()),source_sha256=hashlib.sha256(raw).hexdigest(),
                ok=not any(i['severity']=='error' for i in issues),atom_count=atom_count,residue_count=len(residues),
                model_count=len(model_tags) or (1 if atom_count else 0),hydrogen_count=hydrogen_count,
                issues=issues,messages=[i['message']+' '+i['action'] for i in issues],residues=residues,
                disulfide_candidates=candidates,ssbond_declarations=declarations,
                limitations='Coordinate-only check; missing sequence outside supplied residues, biological assembly, bond orders and force-field validity are not established.')
    if output_path is not None:
        output_path=Path(output_path)
        if output_path.resolve()==path.resolve():
            raise ValueError('Preflight output must not overwrite the input PDB')
        with output_path.open('x') as handle:
            json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    return report
