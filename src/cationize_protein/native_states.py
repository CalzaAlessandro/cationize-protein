"""Rebuild supported native side-chain hydrogens and reassign the force field.

This is a physical topology change, never residue renaming alone. Native N/C
terminal H atoms, disulfides, and all amidation fragments are held fixed.
Hydrogen placement uses OpenMM's generic H-only relaxation; the caller must
subsequently minimize the rebuilt product with the selected force field.
"""
from __future__ import annotations

import math
import numbers


class UnsupportedNativeState(ValueError):
    """A proposed state cannot be reconstructed with the validated policy."""


_FAMILIES=(frozenset(('HID','HIE','HIP')),frozenset(('ASP','ASH')),
           frozenset(('GLU','GLH')),frozenset(('LYS','LYN')),
           frozenset(('CYS','CYM')),frozenset(('TYR','TYM')))
SUPPORTED_NATIVE_STATES=frozenset().union(*_FAMILIES)
_CHARGE=dict(HID=0,HIE=0,HIP=1,ASP=-1,ASH=0,GLU=-1,GLH=0,
             LYS=1,LYN=0,CYS=0,CYM=-1,TYR=0,TYM=-1)
_BASE=dict(HID='HIS',HIE='HIS',HIP='HIS',ASH='ASP',GLH='GLU',LYN='LYS',
           CYM='CYS',CYX='CYS',TYM='TYR')
_SIDE_PARENTS={**dict.fromkeys(('HID','HIE','HIP'),frozenset(('ND1','NE2'))),
               **dict.fromkeys(('ASP','ASH'),frozenset(('OD1','OD2'))),
               **dict.fromkeys(('GLU','GLH'),frozenset(('OE1','OE2'))),
               **dict.fromkeys(('LYS','LYN'),frozenset(('NZ',))),
               **dict.fromkeys(('CYS','CYM'),frozenset(('SG',))),
               **dict.fromkeys(('TYR','TYM'),frozenset(('OH',)))}
_SIDE_H=dict(HID=[('HD1','ND1')],HIE=[('HE2','NE2')],HIP=[('HD1','ND1'),('HE2','NE2')],
             ASP=[],ASH=[('HD2','OD2')],GLU=[],GLH=[('HE2','OE2')],
             LYS=[('HZ1','NZ'),('HZ2','NZ'),('HZ3','NZ')],LYN=[('HZ2','NZ'),('HZ3','NZ')],
             CYS=[('HG','SG')],CYM=[],TYR=[('HH','OH')],TYM=[])


def _info(structure):
    if hasattr(structure,'peptide_bonds') and hasattr(structure,'chains'):
        return structure
    from .topology import inspect_structure
    return inspect_structure(structure)


def plan_native_states(structure, changes):
    """Validate one-based ordinal->state requests without OpenMM or mutation."""
    info=_info(structure);by_ordinal={r.ordinal:r for r in info.residues};rows=[];delta=0;seen=set()
    for key,target in changes.items():
        try:
            if isinstance(key,bool):raise ValueError
            if isinstance(key,numbers.Integral):ordinal=int(key)
            elif isinstance(key,str) and key.isdecimal():ordinal=int(key)
            else:raise ValueError
            r=by_ordinal[ordinal]
        except (ValueError,KeyError) as exc:
            raise UnsupportedNativeState(f'Invalid one-based residue ordinal {key!r}') from exc
        if ordinal in seen:raise UnsupportedNativeState(f'Duplicate requested residue ordinal {ordinal}')
        seen.add(ordinal)
        if r.modified_groups:
            raise UnsupportedNativeState(f'Protonation change on modified residue ordinal {ordinal} needs a matching validated fitted template; hold this proposal')
        source=r.canonical_name
        if source==target:continue
        family=next((family for family in _FAMILIES if source in family),None)
        if family is None or target not in family:
            raise UnsupportedNativeState(f'Unsupported native state transition at ordinal {ordinal}: {source} -> {target}')
        change=_CHARGE[target]-_CHARGE[source];delta+=change
        rows.append(dict(residue_ordinal=ordinal,source_state=source,target_state=target,
                         chain_id=r.chain_id,original_id=r.original_id,expected_charge_change_e=change))
    return dict(changes=rows,expected_charge_change_e=delta,
                supported_native_states=sorted(SUPPORTED_NATIVE_STATES))



def supported_native_contexts(force_field='ff14SB'):
    """Read actual installed native templates; unavailable contexts are holds."""
    import xml.etree.ElementTree as E
    from .forcefields import data_directory,profile
    root=E.parse(data_directory()/profile(force_field)['xml']).getroot()
    names={r.get('name') for r in root.findall('Residues/Residue')}
    return {state:[context for context,prefix in [('internal',''),('n_terminal','N'),('c_terminal','C')]
                   if prefix+state in names] for state in SUPPORTED_NATIVE_STATES}


def apply_native_states(structure, changes, parameters_dir, force_field='ff14SB'):
    """Return (typed ParmEd structure, OpenMM System, validated change record).

    Only native side-chain H atoms on requested residues are rebuilt. Every other
    atom and coordinate is retained. Explicit (hydrogen,parent) variants prevent
    Modeller's automatic pH choices from changing untargeted native residues or
    introduced ammonium groups. CYM and TYM use an explicit no-SG/no-OH-H list;
    names are restored before force-field assignment, so these remain charged
    thiolate/phenolate states. Modified products are held explicitly.
    """
    import numpy as np
    import parmed
    from openmm import Platform,unit
    from openmm.app import Modeller,Topology,element
    from .parameter_assignment import parameterize
    from .protonation import prepare_input

    info=_info(structure);plan=plan_native_states(info,changes)
    ordered=[ordinal for chain in info.chains for ordinal in chain]
    if ordered!=list(range(1,len(info.residues)+1)):
        raise UnsupportedNativeState('Interleaved peptide-chain residue order is unsupported for native-state rebuilding')
    # Confirm all introduced amines have three H atoms before reconstruction.
    _,before_mapping=prepare_input(info)
    target={row['residue_ordinal']:row['target_state'] for row in plan['changes']}
    contexts=supported_native_contexts(force_field)
    for ordinal,state in target.items():
        context=info.residues[ordinal-1].context
        if context not in contexts[state]:
            raise UnsupportedNativeState(f'{force_field} has no {context} native template for {state}; hold the proposal or validate an explicit compatible template')
    old_atoms=list(structure.atoms);old_residues=list(structure.residues)
    old_xyz=np.asarray(structure.coordinates,dtype=float)
    old_atom_index={id(a):i for i,a in enumerate(old_atoms)}
    old_res_index={id(r):i+1 for i,r in enumerate(old_residues)}
    source_keys={(old_res_index[id(a.residue)],a.name):i for i,a in enumerate(old_atoms)}
    if len(source_keys)!=len(old_atoms):raise UnsupportedNativeState('Duplicate atom identities in source structure')
    source_charge=sum(float(a.charge) for a in old_atoms)
    if not math.isfinite(source_charge):raise UnsupportedNativeState('Source charges are non-finite')
    neighbors={i:[] for i in range(len(old_atoms))}
    for bond in structure.bonds:
        i,j=old_atom_index[id(bond.atom1)],old_atom_index[id(bond.atom2)]
        neighbors[i].append(j);neighbors[j].append(i)

    topology=Topology();new_atoms={};new_residues={};variants=[None]*len(old_residues);delete=[];replaced_hydrogens=set()
    for chain in info.chains:
        ch=topology.addChain(info.residues[chain[0]-1].chain_id)
        for ordinal in chain:
            source_r=old_residues[ordinal-1];state=info.residues[ordinal-1]
            name=_BASE.get(source_r.name,source_r.name) if not state.modified_groups else state.canonical_name
            residue=topology.addResidue(name,ch,str(source_r.number),str(source_r.insertion_code or ''))
            new_residues[ordinal]=residue
            hydrogens=[]
            for atom in source_r.atoms:
                i=old_atom_index[id(atom)]
                new=topology.addAtom(atom.name,element.Element.getByAtomicNumber(atom.atomic_number),residue)
                new_atoms[i]=new
                if atom.atomic_number==1:
                    parents=[j for j in neighbors[i] if old_atoms[j].atomic_number!=1]
                    if len(parents)!=1:raise UnsupportedNativeState('Every source hydrogen must have one heavy-atom parent')
                    parent=old_atoms[parents[0]].name
                    if ordinal in target and parent in _SIDE_PARENTS[target[ordinal]]:
                        delete.append(new);replaced_hydrogens.add((ordinal,atom.name))
                    else:
                        hydrogens.append((atom.name,parent))
            if ordinal in target:
                hydrogens.extend(_SIDE_H[target[ordinal]])
            variants[residue.index]=hydrogens
    for bond in structure.bonds:
        topology.addBond(new_atoms[old_atom_index[id(bond.atom1)]],new_atoms[old_atom_index[id(bond.atom2)]])
    # Build positions in the exact order of the explicitly constructed topology.
    positions=[old_xyz[old_atom_index[id(a)]] for ordinal in (r for chain in info.chains for r in chain)
               for a in old_residues[ordinal-1].atoms]
    modeller=Modeller(topology,np.asarray(positions)*unit.angstrom)
    modeller.delete(delete)
    modeller.addHydrogens(forcefield=None,pH=7.4,variants=variants,platform=Platform.getPlatformByName('Reference'))
    rebuilt=parmed.openmm.load_topology(modeller.topology,xyz=modeller.positions.value_in_unit(unit.angstrom))
    if len(rebuilt.residues)!=len(old_residues):raise UnsupportedNativeState('Hydrogen rebuilding changed residue membership')
    # Modeller follows the source graph chain order. ProteinTopology ordinals are
    # input residue order; ensure this invariant before restoring identities.
    ordered=[ordinal for chain in info.chains for ordinal in chain]
    if ordered!=list(range(1,len(old_residues)+1)):
        raise UnsupportedNativeState('Interleaved peptide-chain residue order is unsupported for native-state rebuilding')
    for ordinal,(old,new) in enumerate(zip(old_residues,rebuilt.residues),1):
        new.name=target.get(ordinal,old.name);new.chain=old.chain;new.number=old.number
        new.insertion_code=old.insertion_code;new.segid=old.segid;new.ter=old.ter
    if structure.box is not None:rebuilt.box=np.array(structure.box,copy=True)

    def keys(s):return {(a.residue.idx+1,a.name):a for a in s.atoms}
    rebuilt_keys=keys(rebuilt)
    # Existing atoms are immobile in Modeller. Restore unit conversion roundoff
    # after checking its magnitude; all heavy coordinates then remain bit exact.
    for key,index in source_keys.items():
        if key not in rebuilt_keys or key in replaced_hydrogens:continue
        atom=rebuilt_keys[key]
        xyz=np.array([atom.xx,atom.xy,atom.xz])
        if not np.allclose(xyz,old_xyz[index],rtol=0,atol=1e-10):
            raise UnsupportedNativeState(f'Hydrogen placement moved an existing atom: {key}')
        atom.xx,atom.xy,atom.xz=old_xyz[index]
    typed,system=parameterize(rebuilt,parameters_dir,force_field=force_field,implicit=False)
    typed_keys=keys(typed)
    def heavy_graph(s):
        return {tuple(sorted(((b.atom1.residue.idx+1,b.atom1.name),(b.atom2.residue.idx+1,b.atom2.name))))
                for b in s.bonds if b.atom1.atomic_number!=1 and b.atom2.atomic_number!=1}
    if heavy_graph(typed)!=heavy_graph(structure):
        raise UnsupportedNativeState('Native-state rebuilding changed the heavy-atom graph or external/disulfide bonds')
    for key,index in source_keys.items():
        old=old_atoms[index]
        if old.atomic_number!=1 and key not in typed_keys:
            raise UnsupportedNativeState(f'Native-state rebuilding lost a heavy atom: {key}')
        if key in typed_keys and key not in replaced_hydrogens:
            atom=typed_keys[key]
            if not np.array_equal(np.array([atom.xx,atom.xy,atom.xz]),old_xyz[index]):
                raise UnsupportedNativeState(f'Parameter assignment changed a retained coordinate: {key}')
    final_charge=sum(float(a.charge) for a in typed.atoms)
    observed_delta=final_charge-source_charge
    if not math.isfinite(final_charge) or abs(observed_delta-plan['expected_charge_change_e'])>2e-4:
        raise UnsupportedNativeState(f'Native-state charge change is {observed_delta:+.6f} e; expected {plan["expected_charge_change_e"]:+d} e')
    _,after_mapping=prepare_input(typed)
    original_amines={(a['residue_ordinal'],a['atom'],a['fixed_charge_e']) for a in before_mapping['introduced_amines']}
    final_amines={(a['residue_ordinal'],a['atom'],a['fixed_charge_e']) for a in after_mapping['introduced_amines']}
    if final_amines!=original_amines:raise UnsupportedNativeState('Introduced +1 ammonium groups changed during native-state rebuilding')
    record=dict(plan,force_field=force_field,source_charge_e=source_charge,final_charge_e=final_charge,
                observed_charge_change_e=observed_delta,heavy_coordinates_exact=True,
                heavy_and_external_graph_exact=True,introduced_ammonium_groups_retained=len(final_amines),
                native_terminal_states_held=True,
                method='Explicit OpenMM hydrogen-parent variants; generic H-only placement, then complete force-field reassignment',
                next_required_step='Minimize with the selected force field and repeat the two-control PROPKA audit.')
    return typed,system,record
