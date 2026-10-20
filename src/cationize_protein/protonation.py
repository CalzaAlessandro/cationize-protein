"""Conservative HSA-style PROPKA diagnostics for fixed-ammonium products.

PROPKA is an empirical diagnostic, not a validation of these nonstandard groups.
Two controls reproduce the HSA policy: generic primary amine versus permanent
+1 ion centre; consensus must be at least one pH unit away and repeat on two
consecutive chemical steps. Introduced ammonium groups and native terminal
charges are never changed. The caller applies/rebuilds accepted native states.
"""
from __future__ import annotations

import contextlib
import hashlib
import importlib.metadata
import io
import json
import logging
import math
from pathlib import Path
import string

from .catalog import templates


class ProtonationError(ValueError):
    """Chemistry, mapping or model controls could not be audited safely."""


_ALIASES={'HID':'HIS','HIE':'HIS','HIP':'HIS','CYX':'CYS','CYM':'CYS',
          'ASH':'ASP','GLH':'GLU','LYN':'LYS','TYM':'TYR'}
_HISTIDINES=frozenset(('HID','HIE','HIP'))


def _json(path,record):
    Path(path).write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


def _topology(structure):
    if hasattr(structure,'peptide_bonds') and hasattr(structure,'chains'):
        return structure
    from .topology import inspect_structure
    return inspect_structure(structure)


def prepare_input(structure):
    """Return (heavy-PDB text, exact ordinal/atom/graph mapping) without writes.

    Unique global residue ordinals and synthetic single-character chain IDs avoid
    collisions from repeated PDB residue numbers. Original identities are kept in
    the mapping. Native backbone/unmodified side groups use PROPKA protein rules;
    converted acidic side chains and added fragments use ligand/group rules.
    This avoids treating a converted acid as an incomplete carboxylate and keeps
    native histidine/cysteine groups on C-terminal products visible to PROPKA.
    """
    info=_topology(structure)
    if len(info.residues)>9999 or len(info.atoms)>99999:
        raise ProtonationError('PROPKA PDB index limits exceeded; split input only with a reviewed electrostatic policy')
    if len(info.chains)>62:
        raise ProtonationError('PROPKA supports at most 62 distinct mapped PDB chain identifiers here')
    chain_ids=string.ascii_uppercase+string.ascii_lowercase+string.digits
    atoms={a.index:a for a in info.atoms}
    graph={i:set() for i in atoms}
    for a,b in info.bonds:
        graph[a].add(b);graph[b].add(a)
    lines=[]; mapping=[]; expected=[]; amines=[]; nterms=[]; heavy_keys={}; atom_mapping=[]
    known={s['name']:s for s in templates()}
    for r in info.residues:
        spec=known.get(r.canonical_name)
        if r.modified_groups and spec is None:
            raise ProtonationError(f'Unknown modified template at residue ordinal {r.ordinal}: {r.canonical_name}')
        row=dict(residue_ordinal=r.ordinal,chain_id=r.chain_id,chain_ordinal=r.chain_ordinal,
                 original_id=r.original_id,source_state=r.canonical_name,
                 native_state=r.parent if r.modified_groups else r.canonical_name,
                 context=r.context,modified_groups=list(r.modified_groups),modified_residue=bool(r.modified_groups),
                 mapped_chain_id=chain_ids[r.chain_ordinal-1])
        mapping.append(row)
        native_inverse={actual:native for native,actual in r.native_atom_names.items()}
        backbone={r.native_atom_indices[n] for n in ('N','CA','C','O','OXT') if n in r.native_atom_indices}
        introduced={}
        reused_fragments={'ASL':{'N2','C5','C6','N3'},'NSL':{'N2','C5','C6','N3'},
                          'GLS':{'N2','C6','C7','N3'},'LCL':{'N2','C3','C4','N3'}}
        added_names=(reused_fragments.get(spec['reused_from'],set()) if spec else set())
        if spec and not spec['reused_from']:
            added_names={f'{prefix}{i}' for group in r.modified_groups
                         for prefix in ('S' if group=='side' else 'T',) for i in range(1,5)}
        for group in r.modified_groups:
            name='N3' if spec and spec['reused_from'] else 'S4' if group=='side' else 'T4'
            matching=[i for i in r.atom_indices if atoms[i].name==name]
            if len(matching)!=1 or atoms[matching[0]].atomic_number!=7:
                raise ProtonationError(f'Added {group} ammonium nitrogen is missing/ambiguous at residue ordinal {r.ordinal}: {name}')
            i=matching[0]; hydrogen_neighbors=[j for j in graph[i] if atoms[j].atomic_number==1]
            heavy_neighbors=[j for j in graph[i] if atoms[j].atomic_number!=1]
            if len(hydrogen_neighbors)!=3 or len(heavy_neighbors)!=1 or atoms[heavy_neighbors[0]].atomic_number!=6:
                raise ProtonationError(f'Added amine {r.ordinal}:{name} must retain one carbon bond and three hydrogens (+1 ammonium invariant)')
            introduced[i]=group
            amines.append(dict(residue_ordinal=r.ordinal,atom=name,modified_group=group,
                               fixed_charge_e=1,chain_id=r.chain_id,original_id=r.original_id))
        seen=set()
        for local_i,i in enumerate(r.atom_indices,1):
            a=atoms[i]
            if a.atomic_number==1:
                continue
            original_name=a.name
            # Legacy templates have a complete native->template map; added atoms
            # retain their actual names. New templates already use native names.
            name=native_inverse.get(original_name,original_name)
            try:
                element={6:'C',7:'N',8:'O',16:'S',15:'P',9:'F',17:'CL',35:'BR',53:'I'}[a.atomic_number]
            except KeyError as exc:
                raise ProtonationError('Unsupported heavy element in PROPKA input') from exc
            # PROPKA 3.x infers elements from atom-name columns and ignores the
            # PDB element field. S4/T4 are nitrogens, not sulfur: map them to
            # element-correct, reversible labels before invoking its parser.
            if name.lstrip(string.digits)[:len(element)].upper()!=element:
                if local_i>999:raise ProtonationError('Too many atoms in one residue for PROPKA name mapping')
                name=f'{element}{local_i:03d}'
            if name in seen or len(name)>4:
                raise ProtonationError(f'Ambiguous/unrepresentable PROPKA atom name at residue {r.ordinal}: {name}')
            seen.add(name)
            atom_mapping.append(dict(residue_ordinal=r.ordinal,atom_index=i,source_atom=original_name,propka_atom=name,element=element))
            if i in introduced:
                next(a for a in amines if a['residue_ordinal']==r.ordinal and a['atom']==original_name)['propka_atom']=name
            native_atom=(i in backbone or (original_name not in added_names and 'side' not in r.modified_groups))
            if i in introduced:
                native_atom=False
            record='ATOM  ' if native_atom else 'HETATM'
            resname=_ALIASES.get(r.parent,r.parent) if native_atom else r.canonical_name
            if native_atom and 'side' in r.modified_groups:
                resname='ALA'  # Only the standard backbone remains protein-typed.
            if r.context=='n_terminal' and i==r.native_atom_indices['N']:
                nterms.append(dict(residue_ordinal=r.ordinal,atom=name,source_atom=original_name,held=True))
            try:
                xyz=[float(x) for x in info.coordinates[i]]
                if len(xyz)!=3 or not all(math.isfinite(x) for x in xyz):raise ValueError
                element={6:'C',7:'N',8:'O',16:'S',15:'P',9:'F',17:'CL',35:'BR',53:'I'}[a.atomic_number]
            except (ValueError,KeyError,TypeError) as exc:
                raise ProtonationError('Finite coordinates and supported explicit heavy elements are required') from exc
            if any(x < -999.999 or x>9999.999 for x in xyz):
                raise ProtonationError('Coordinates exceed fixed-column PDB limits')
            serial=len(lines)+1
            name_field=(' '+name.ljust(3)) if len(element)==1 and len(name)<4 else name.ljust(4)
            lines.append(f'{record}{serial:5d} {name_field} {resname:>3s} {row["mapped_chain_id"]}{r.ordinal:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}')
            heavy_keys[i]=(r.ordinal,name)
        if r.ordinal in {chain[-1] for chain in info.chains}:
            lines.append('TER   ')
    for a,b in info.bonds:
        if a in heavy_keys and b in heavy_keys:
            expected.append(sorted([list(heavy_keys[a]),list(heavy_keys[b])]))
    text='\n'.join(lines)+'\nEND\n'
    manifest=dict(residues=mapping,atom_mapping=atom_mapping,introduced_amines=amines,native_n_termini=nterms,
                  expected_heavy_bonds=sorted(expected),expected_heavy_atoms=sorted([list(v) for v in heavy_keys.values()]))
    digest=dict(pdb=text,mapping=mapping,graph=manifest['expected_heavy_bonds'])
    manifest['structure_sha256']=hashlib.sha256(json.dumps(digest,sort_keys=True).encode()).hexdigest()
    return text,manifest


def calculate(structure, output_dir, *, pH=7.4, fixed=False):
    """Run the installed PROPKA Python API; fail on unrecognized added amines.

    The exact inferred heavy-atom graph must match the authoritative source
    topology before pKa calculation. This intentionally blocks bad geometry,
    extra inferred close-contact bonds and missing interchain/disulfide links.
    The sensitivity control keeps a +1 charge and omits amine-specific H-bond/
    titration terms. Neither control uses RESP atomic charges.
    """
    pH=float(pH)
    if not math.isfinite(pH):
        raise ValueError('pH must be finite')
    try:
        from propka.input import read_parameter_file, read_molecule_file
        from propka.lib import loadOptions
        from propka.molecular_container import MolecularContainer
        from propka.parameters import Parameters
    except ImportError as exc:
        raise ImportError('PROPKA audit requires propka>=3.4,<4 in the project Python environment') from exc
    pdb,mapping=prepare_input(structure)
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    inp=out/'input.pdb';inp.write_text(pdb);_json(out/'mapping.json',mapping)
    # Keep the host application's root logging configuration intact.
    logger=logging.getLogger('propka');handler=logging.FileHandler(out/'propka.log',mode='w')
    old_level=logger.level;old_propagate=logger.propagate
    logger.addHandler(handler);logger.setLevel(logging.INFO);logger.propagate=False
    try:
        opts=loadOptions((str(inp),'--pH',str(pH)))
        parameters=read_parameter_file(opts.parameters,Parameters())
        with (out/'stdout.txt').open('w') as stdout,contextlib.redirect_stdout(stdout):
            mol=read_molecule_file(str(inp),MolecularContainer(parameters,opts),stream=io.StringIO(pdb))
            if len(mol.conformation_names)!=1:
                raise ProtonationError('PROPKA input unexpectedly created multiple conformations')
            conf=mol.conformations[mol.conformation_names[0]]
            observed_atoms=sorted([[a.res_num,a.name] for a in conf.atoms if a.element!='H'])
            if observed_atoms!=mapping['expected_heavy_atoms']:
                raise ProtonationError('PROPKA atom mapping differs from the source topology')
            observed=set()
            for a in conf.atoms:
                if a.element=='H':continue
                for b in a.bonded_atoms:
                    if b.element!='H':observed.add(tuple(sorted(((a.res_num,a.name),(b.res_num,b.name)))))
            expected={tuple(tuple(x) for x in pair) for pair in mapping['expected_heavy_bonds']}
            if observed!=expected:
                differences=dict(missing_bonds=sorted(expected-observed),extra_bonds=sorted(observed-expected))
                _json(out/'graph_mismatch.json',differences)
                raise ProtonationError(f'PROPKA inferred heavy graph differs: {len(expected-observed)} missing and {len(observed-expected)} extra bonds; see {out / "graph_mismatch.json"}')
            amines={(a['residue_ordinal'],a['propka_atom']) for a in mapping['introduced_amines']}
            recognized={(g.atom.res_num,g.atom.name) for g in conf.groups if g.type=='N31' and (g.atom.res_num,g.atom.name) in amines}
            if recognized!=amines:
                raise ProtonationError(f'PROPKA failed to recognize introduced generic primary amines: {sorted(amines-recognized)}')
            nterms={(r['residue_ordinal'],r['atom']) for r in mapping['native_n_termini']}
            actual_nterms={(g.atom.res_num,g.atom.name) for g in conf.groups if g.type=='N+'}
            if actual_nterms!=nterms:
                raise ProtonationError('PROPKA native N+ groups do not match the source peptide chain termini')
            expected_native=set()
            for r in mapping['residues']:
                parent=_ALIASES.get(r['native_state'],r['native_state'])
                if parent in {'ASP','GLU','HIS','CYS','TYR','LYS','ARG'} and not (parent in {'ASP','GLU'} and 'side' in r['modified_groups']):
                    expected_native.add((r['residue_ordinal'],parent))
                if r['context']=='c_terminal' and 'terminal' not in r['modified_groups']:
                    expected_native.add((r['residue_ordinal'],'C-'))
            observed_native={(g.atom.res_num,g.residue_type.strip()) for g in conf.groups}
            if not expected_native.issubset(observed_native):
                raise ProtonationError(f'PROPKA failed to recognize native ionizable groups: {sorted(expected_native-observed_native)}')
            if fixed:
                for g in conf.groups:
                    if (g.atom.res_num,g.atom.name) in amines:
                        for a in g.interaction_atoms_for_acids+g.interaction_atoms_for_bases:
                            if a.element=='H':a.group_type=None
                        g.type,g.residue_type='ION','1P'
                        g.titratable=g.model_pka_set=False;g.model_pka=0.;g.charge=1
                        g.interaction_atoms_for_acids=g.interaction_atoms_for_bases=[g.atom]
                        g.atom.group_type='ION'
            mol.calculate_pka()
        by_ordinal={r['residue_ordinal']:r for r in mapping['residues']}
        rows=[]
        original_atoms={(a['residue_ordinal'],a['propka_atom']):a['source_atom'] for a in mapping['atom_mapping']}
        for g in mol.conformations['AVR'].groups:
            if not g.titratable or g.atom.cysteine_bridge:continue
            row=dict(by_ordinal[g.atom.res_num],group=g.residue_type.strip(),atom=g.atom.name,
                     pKa=float(g.pka_value),source_atom=original_atoms[g.atom.res_num,g.atom.name],introduced_amine=(g.atom.res_num,g.atom.name) in amines)
            if not math.isfinite(row['pKa']):raise ProtonationError('PROPKA produced a nonfinite pKa')
            row['protonated_fraction_HH']=1/(1+10**max(-100,min(100,pH-row['pKa'])))
            rows.append(row)
        result=dict(pH=pH,propka_version=importlib.metadata.version('propka'),
                    structure_sha256=mapping['structure_sha256'],mapping=mapping,predictions=rows,
                    exact_heavy_graph_match=True,heavy_bonds=len(expected),recognized_added_amines=len(recognized),
                    fixed_added_amines=fixed,
                    model='Permanent +1 ion-centre sensitivity control; amine H-bond/titration terms omitted' if fixed else 'Standard generic primary-amine N31 model',
                    limitations='Single-structure empirical diagnostics; RESP charges are not PROPKA inputs; modified-group model has not been experimentally calibrated.')
        _json(out/'PROPKA_COMPLETE.json',result)
        return result
    except Exception as exc:
        _json(out/'PROPKA_FAILED.json',dict(error_type=type(exc).__name__,message=str(exc)))
        raise
    finally:
        logger.removeHandler(handler);handler.close();logger.setLevel(old_level);logger.propagate=old_propagate


def assess_predictions(standard, fixed, previous_proposals=None, *, pH=7.4, supported_native_states=None):
    """Return proposals and safely applicable ordinal->state changes.

    ``standard`` and ``fixed`` are prediction lists returned by ``calculate``.
    ``previous_proposals`` accepts the preceding full audit record or its
    ``proposed_native_states`` mapping. The caller must pass the immediately
    preceding chemical-step audit, never an older nonconsecutive proposal.
    Default application capabilities are HID/HIE/HIP only; a backend can opt in
    to other native states explicitly. Modified-product native side changes
    require separately validated fitted templates and are always held here.
    Neutral HIP->HIE uses the HSA protocol's documented default; PROPKA itself
    does not establish which neutral histidine tautomer is preferred.
    """
    pH=float(pH)
    if not math.isfinite(pH):raise ValueError('pH must be finite')
    supported=set(_HISTIDINES if supported_native_states is None else supported_native_states)
    previous=previous_proposals or {}
    previous=previous.get('proposed_native_states',previous)
    previous={str(k):v for k,v in previous.items()}
    def key(r):return r['residue_ordinal'],r['group'].strip(),r['atom']
    if len({key(r) for r in standard})!=len(standard) or len({key(r) for r in fixed})!=len(fixed):
        raise ProtonationError('Duplicate PROPKA group identity in model control')
    standard_native={key(r) for r in standard if not r.get('introduced_amine')}
    fixed_native={key(r) for r in fixed if not r.get('introduced_amine')}
    if standard_native!=fixed_native:
        raise ProtonationError('Native group identities differ between PROPKA model controls')
    lookup={key(r):r for r in fixed};proposals={};applicable={};decisions=[]
    for r in standard:
        group=r['group'].strip();state=r.get('native_state',r['source_state'])
        d={k:r.get(k) for k in ('residue_ordinal','chain_id','chain_ordinal','original_id','source_state')}
        d.update(group=group,atom=r['atom'],pKa_standard=float(r['pKa']),target_state=None)
        if r.get('introduced_amine') or group in {'N+','C-'}:
            d.update(status='held_fixed_group',reason='Introduced amines remain +1; supplied native N/C terminal states are held.');decisions.append(d);continue
        control=lookup.get(key(r))
        if control is None:raise ProtonationError(f'Missing fixed-amine control for native group {key(r)}')
        p1,p2=float(r['pKa']),float(control['pKa'])
        if not math.isfinite(p1) or not math.isfinite(p2):raise ProtonationError('Nonfinite pKa in model control')
        d['pKa_fixed_control']=p2
        neutral=max(p1,p2)<=pH-1.;protonated=min(p1,p2)>=pH+1.
        target=None
        if group=='HIS':target='HIP' if protonated else ('HID' if state=='HID' else 'HIE') if neutral else None
        elif group=='ASP':target='ASH' if protonated else 'ASP' if neutral else None
        elif group=='GLU':target='GLH' if protonated else 'GLU' if neutral else None
        elif group=='LYS':target='LYS' if protonated else 'LYN' if neutral else None
        elif group=='CYS':target='CYS' if protonated else 'CYM' if neutral else None
        elif group=='TYR':target='TYR' if protonated else 'TYM' if neutral else None
        elif group=='ARG':target='ARG' if protonated else 'ARN' if neutral else None
        d['target_state']=target
        if not protonated and not neutral:
            d.update(status='held_borderline_or_disagreement',reason='Both controls must support the same side of pH with a margin of at least 1 pH unit.')
        elif target is None:
            d.update(status='held_unmapped_group',reason='No validated native-state transition policy exists for this PROPKA group.')
        elif target==state:
            d.update(status='retained',reason='Both controls support the supplied state; neutral histidine tautomer is retained.')
        else:
            proposals[str(r['residue_ordinal'])]=target
            reason='Both controls agree at least 1 pH unit from target pH.'
            if group=='HIS' and state=='HIP' and target=='HIE':
                reason+=' HIE is the HSA protocol default neutral tautomer; PROPKA does not select this tautomer.'
            if r.get('modified_residue'):
                d.update(status='held_modified_template',reason=reason+' A native-side state change on an amidated residue needs a matching validated fitted template.')
            elif target not in supported:
                d.update(status='held_unsupported_state',reason=reason+' The active backend has not declared this native protonation state supported.')
            elif previous.get(str(r['residue_ordinal']))!=target:
                d.update(status='proposed_awaiting_persistence',reason=reason+' Retain supplied state until the same proposal recurs on the next chemical step.')
            else:
                applicable[str(r['residue_ordinal'])]=target
                d.update(status='applicable',reason=reason+' The same proposal recurred on two consecutive chemical-step audits.')
        decisions.append(d)
    return dict(proposed_native_states=proposals,applicable_native_states=applicable,decisions=decisions,
                target_pH=pH,margin_pKa_units=1.,required_consecutive_steps=2,
                supported_native_states=sorted(supported),
                policy='Retain all held states; rebuild only applicable states, reparameterize, minimize and audit again. Unsupported proposals are explicit holds, never silent residue retyping.')


def audit(structure, output_dir, previous_proposals=None, pH=7.4, *, supported_native_states=None):
    """Run both models and write one reviewable audit; do not mutate structure.

    ``applicable_native_states`` maps one-based global residue ordinals (strings)
    to Amber native names. Source chains/PDB IDs are diagnostic metadata only.
    A held unsupported proposal sets ``has_unsupported_proposals`` for an
    explicit caller policy (hold/report or fail); it never becomes applicable.
    """
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    standard=calculate(structure,out/'standard',pH=pH,fixed=False)
    fixed=calculate(structure,out/'fixed_amine_control',pH=pH,fixed=True)
    result=assess_predictions(standard['predictions'],fixed['predictions'],previous_proposals,
                              pH=pH,supported_native_states=supported_native_states)
    result.update(structure_sha256=standard['structure_sha256'],mapping=standard['mapping'],
                  introduced_amine_charge_e=1,native_terminal_states_held=True,
                  controls={k:{f:r[f] for f in ('propka_version','exact_heavy_graph_match','recognized_added_amines','model','limitations')}
                            for k,r in [('standard',standard),('fixed_amine_control',fixed)]})
    result['has_unsupported_proposals']=any(d['status'] in {'held_modified_template','held_unsupported_state','held_unmapped_group'} for d in result['decisions'])
    _json(out/'proposals.json',result)
    return result
