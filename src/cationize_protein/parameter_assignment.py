"""Profile-aware parameter assignment without changing fitted-template provenance."""
from pathlib import Path
import xml.etree.ElementTree as E
import numpy as np
from .catalog import templates
from .forcefields import profile, load_forcefield, data_directory
from .off import read_units


def _graph(structure):
    return {tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in structure.bonds}


def _assigned_atom_types(ff,topology):
    """Use the actual graph match, including equivalent H-name permutations."""
    data=ff._SystemData(topology)
    ff._matchAllResiduesToTemplates(data,topology,{},False)
    return {a.index:data.atomType[a] for a in topology.atoms()}


def _expected_cmap(ff, topology, force_field):
    if not profile(force_field)['cmap']:
        return 0
    root = E.parse(data_directory()/profile(force_field)['xml']).getroot()
    rules = root.findall('CMAPTorsionForce/Torsion')
    atom_types = _assigned_atom_types(ff,topology)
    graph = {a.index: set() for a in topology.atoms()}
    for first, second in topology.bonds():
        graph[first.index].add(second.index); graph[second.index].add(first.index)
    classes = {i: ff._atomTypes[t].atomClass for i, t in atom_types.items()}
    count = 0
    for middle in graph:
        for second in graph[middle]:
            for fourth in graph[middle]-{second}:
                matching = [r for r in rules if all(atom_types[i] == r.get('type'+str(k))
                    for i,k in ((second,2),(middle,3),(fourth,4)))]
                for rule in matching:
                    firsts = [i for i in graph[second]-{middle,fourth} if classes[i] == rule.get('class1')]
                    fifths = [i for i in graph[fourth]-{middle,second} if classes[i] == rule.get('class5')]
                    count += sum(first != fifth for first in firsts for fifth in fifths)
    return count


def parameterize(structure, folder, force_field='ff14SB', implicit=False, names=None):
    """Return (ParmEd structure, OpenMM System), preserving graph and fitted q."""
    profile(force_field)  # Reject unknown profiles before optional runtime imports.
    import parmed
    from openmm import unit
    from openmm.app import NoCutoff
    from .catalog import by_name
    if names is None:
        known={s['name'] for s in templates()}
        names=sorted({residue.name for residue in structure.residues}&known)
    ff = load_forcefield(folder, force_field, implicit, names)
    topology = structure.topology
    if force_field == 'ff14SB':
        from .openmm_backend import parameterize as original_parameterize
        typed, system = original_parameterize(structure, folder, implicit=implicit, names=names)
    else:
        renamed, residue_templates = {}, {}
        if implicit:
            for residue in topology.residues():
                try: spec = by_name(residue.name)
                except StopIteration: continue
                renamed[residue] = residue.name; residue_templates[residue] = residue.name
                residue.name = spec['parent']
        try:
            system = ff.createSystem(topology, nonbondedMethod=NoCutoff, constraints=None,
                rigidWater=False, removeCMMotion=False, residueTemplates=residue_templates)
        finally:
            for residue, name in renamed.items(): residue.name = name
        typed = parmed.openmm.load_topology(topology, system,
            xyz=np.asarray(structure.coordinates)*unit.angstrom)
        from parmed.tools.actions import changeRadii
        changeRadii(typed, 'mbondi2').execute()
    if len(typed.atoms) != len(structure.atoms) or _graph(typed) != _graph(structure):
        raise ValueError('Parameter assignment changed the covalent graph or atom order')
    identity = [(a.residue.name, a.name, a.atomic_number) for a in structure.atoms]
    if [(a.residue.name, a.name, a.atomic_number) for a in typed.atoms] != identity:
        raise ValueError('Parameter assignment changed atom identity or order')
    known = {s['name'] for s in templates()}; errors = []; custom = 0
    for residue in typed.residues:
        if residue.name not in known: continue
        custom += 1
        unit_data = read_units(Path(folder)/(residue.name+'.lib'))[residue.name]
        charges = {a['name']: a['initial_template_charge_e'] for a in unit_data['atoms']}
        if set(charges) != {a.name for a in residue.atoms}:
            raise ValueError('Custom template atom names differ after assignment')
        errors.extend(abs(a.charge-charges[a.name]) for a in residue.atoms)
    error = max(errors, default=0.)
    if error > 1e-7:
        raise ValueError('Parameter assignment changed fitted RESP charges')
    actual_cmap = sum(f.getNumTorsions() for f in system.getForces() if type(f).__name__ == 'CMAPTorsionForce')
    expected_cmap = _expected_cmap(ff, topology, force_field)
    if actual_cmap != expected_cmap:
        raise ValueError(f'CMAP count mismatch: {actual_cmap} != {expected_cmap}')
    typed._cationize_system = system
    template_types=_assigned_atom_types(ff,topology)
    typed._cationize_atom_types=[template_types[i] for i in range(len(typed.atoms))]
    typed._cationize_assignment = dict(force_field=force_field, implicit_GBn2=bool(implicit),
        custom_residues=custom, max_custom_charge_error_e=error, cmap_terms=actual_cmap,
        expected_cmap_terms=expected_cmap, custom_charge_reference='ff14SB RESP',
        transfer='native parent-context types; fitted custom charges retained')
    return typed, system
