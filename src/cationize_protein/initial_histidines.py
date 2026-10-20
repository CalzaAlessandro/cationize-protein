"""Initial native histidine charge assignment with PROPKA and OpenMM tautomers."""
from __future__ import annotations

import csv
import io
from pathlib import Path

from .preflight import check_pdb


def _remove_connectivity_atoms(line, deleted):
    if not line.startswith('CONECT') or not deleted:
        return line
    numbers=[int(line[i:i+5]) for i in range(6,len(line.rstrip()),5) if line[i:i+5].strip()]
    if not numbers or numbers[0] in deleted:
        return ''
    kept=[n for n in numbers[1:] if n not in deleted]
    return 'CONECT'+''.join(f'{n:5d}' for n in [numbers[0],*kept])+'\n' if kept else ''


def assign_histidines(pdb_path, output_dir, *, pH=7.4):
    """Return PDB text and per-His decisions; write a persistent PROPKA audit.

    PROPKA determines HIP versus neutral using pKa > pH. Neutral HID/HIE is
    retained when supplied; otherwise OpenMM chooses by H-bond geometry.
    Only histidine labels and ring hydrogens change in the returned PDB.
    Other native protonation changes remain in the later workflow audit.
    """
    from .resources import default_parameters
    from .workflow import _pdb_structure
    from .native_states import apply_native_states
    from .protonation import calculate

    pdb_path,output_dir=Path(pdb_path),Path(output_dir)
    report=check_pdb(pdb_path)
    histidines=[r for r in report['residues'] if r['residue_name'] in {'HIS','HID','HIE','HIP'}]
    original=pdb_path.read_text()
    if not histidines:
        return original,[]
    lines=original.splitlines(keepends=True)
    hist_lines={a['line']:r for r in histidines for a in r['atoms']}
    deleted={int(lines[a['line']-1][6:11]) for r in histidines for a in r['atoms'] if a['name'] in {'HD1','HE2'}}
    # Obtain a neutral tautomer for ambiguous HIS or HIP before the pKa audit.
    # This is a temporary parameterized model; its atomic charges are not
    # inputs to PROPKA, which receives the exact heavy-atom graph.
    neutral=[]
    for number,line in enumerate(lines,1):
        if number in hist_lines:
            if line[12:16].strip() in {'HD1','HE2'}:continue
            line=line[:17]+'HIS'+line[20:]
        elif line.startswith('ANISOU') and int(line[6:11]) in deleted:
            continue
        neutral.append(_remove_connectivity_atoms(line,deleted))
    output_dir.mkdir(parents=True,exist_ok=False)
    model=output_dir/'neutral_model.pdb';model.write_text(''.join(neutral))
    structure,_=_pdb_structure(model,default_parameters(),pH=7.4,allow_automatic_histidines=True)
    prediction=calculate(structure,output_dir/'audit',pH=pH)
    rows=[];targets={}
    for residue in histidines:
        ordinal=residue['ordinal']
        hits=[row for row in prediction['predictions'] if row['residue_ordinal']==ordinal and row['group']=='HIS']
        if len(hits)!=1:raise ValueError(f'PROPKA did not return one pKa for histidine ordinal {ordinal}')
        pka=hits[0]['pKa'];names={a['name'] for a in residue['atoms']}
        supplied=residue['residue_name']
        diagnostic=names & {'HD1','HE2'}
        if supplied in {'HID','HIE'}:
            tautomer=supplied;method='Supplied neutral tautomer'
        elif len(diagnostic)==1:
            tautomer='HID' if 'HD1' in diagnostic else 'HIE';method='Supplied diagnostic ring hydrogen'
        else:
            tautomer=structure.residues[ordinal-1].name;method='OpenMM hydrogen-bond geometry'
        if tautomer not in {'HID','HIE'}:raise ValueError('Initial neutral histidine selection failed')
        state='HIP' if pka>pH else tautomer
        targets[str(ordinal)]=state
        rows.append(dict(residue_ordinal=ordinal,chain_id=residue['chain_id'],residue_id=residue['original_id'],
                         source_state=supplied,pH=pH,pKa=pka,state=state,neutral_tautomer=tautomer,
                         neutral_tautomer_method=method,borderline=abs(pka-pH)<.5))
    final,_,_=apply_native_states(structure,targets,default_parameters())
    from openmm.app import PDBFile
    from openmm import unit
    rendered=io.StringIO();PDBFile.writeFile(final.topology,final.coordinates*unit.angstrom,rendered,keepIds=True)
    generated=[line for line in rendered.getvalue().splitlines() if line.startswith(('ATOM','HETATM'))]
    final_atoms=list(final.atoms)
    if len(generated)!=len(final_atoms):raise ValueError('Histidine serialization changed atom count')
    records={(a.residue.idx+1,a.name):line for a,line in zip(final_atoms,generated)}
    states={int(key):state for key,state in targets.items()};ring={}
    for r in histidines:
        state=states[r['ordinal']]
        ring[r['ordinal']]={'HD1','HE2'} if state=='HIP' else {'HD1'} if state=='HID' else {'HE2'}
    discarded={int(lines[a['line']-1][6:11]) for r in histidines for a in r['atoms']
               if a['name'] in {'HD1','HE2'} and a['name'] not in ring[r['ordinal']]}
    serial=max(int(line[6:11]) for line in lines if line.startswith(('ATOM','HETATM')))
    additions={}
    for r in histidines:
        source_line=lines[r['atoms'][0]['line']-1];present={a['name'] for a in r['atoms']}
        for name in sorted(ring[r['ordinal']]-present):
            serial+=1
            if serial>99999:raise ValueError('Initial histidine H exceeds the PDB serial limit')
            record=records[r['ordinal'],name]
            record=record[:6]+f'{serial:5d}'+record[11:17]+f'{states[r["ordinal"]]:>3s}'+source_line[20:27]+record[27:]+'\n'
            additions.setdefault(r['atoms'][-1]['line'],[]).append(record)
    changed=[]
    identities={(r['chain_id'],r['original_id']):states[r['ordinal']] for r in histidines}
    for number,line in enumerate(lines,1):
        removed=line.startswith(('ATOM','HETATM','ANISOU')) and int(line[6:11]) in discarded
        if not removed:
            if number in hist_lines:
                line=line[:17]+f'{states[hist_lines[number]["ordinal"]]:>3s}'+line[20:]
            elif line.startswith('ANISOU'):
                identity=(line[21].strip(),line[22:26].strip()+line[26].strip())
                if identity in identities:line=line[:17]+f'{identities[identity]:>3s}'+line[20:]
            changed.append(_remove_connectivity_atoms(line,discarded))
        changed.extend(additions.get(number,[]))
    with (output_dir/'histidines.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return ''.join(changed),rows
