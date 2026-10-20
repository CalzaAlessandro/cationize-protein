"""Protein force-field profiles and explicit transfer of custom RESP products.

OpenMM is an engine, not a force field. Product RESP charges are retained across
profiles; alternate profiles transfer native protein parameters by analogy.
ff19SB backbone type identity is inherited from the parent residue so CMAP is
not silently lost when the product has a new three-letter residue name.
"""
from pathlib import Path
import xml.etree.ElementTree as E

_PROFILES = {
    'ff14SB': dict(xml='amber14/protein.ff14SB.xml', water='TIP3P', cmap=False),
    'ff19SB': dict(xml='amber19/protein.ff19SB.xml', water='OPC', cmap=True),
    'ff99SBildn': dict(xml='amber99sbildn.xml', water='TIP3P', cmap=False),
}


def profile(name):
    if name not in _PROFILES:
        raise ValueError('Unknown protein force field: '+str(name))
    return dict(name=name, **_PROFILES[name])


def data_directory():
    import openmm.app
    return Path(openmm.app.__file__).parent/'data'


def force_field_status(name):
    spec = profile(name)
    try:
        exists = (data_directory()/spec['xml']).is_file()
    except ImportError:
        exists = False
    return dict(**spec, available=exists, parameter_family='Amber protein + fixed custom RESP charges',
                validation_scope='ff14SB HSA/capped-model checks; alternate fields require profile-specific transfer/export validation')


def custom_xml(folder, force_field, names=None, *, data_dir=None):
    from .catalog import by_name
    from .off import read_units
    from .topology import _REUSED_MAPS
    native = E.parse(Path(data_dir or data_directory())/profile(force_field)['xml']).getroot()
    types = {t.get('name'): t for t in native.find('AtomTypes')}
    templates = {r.get('name'): {a.get('name'): a.get('type') for a in r.findall('Atom')}
                 for r in native.find('Residues')}
    nb = native.find('NonbondedForce')
    per_residue_charge = any(a.tag == 'UseAttributeFromResidue' and a.get('name') == 'charge' for a in nb)
    root = E.Element('ForceField')
    custom_types = E.SubElement(root, 'AtomTypes')
    residues = E.SubElement(root, 'Residues')
    output_nb = None
    if not per_residue_charge:
        output_nb = E.SubElement(root, 'NonbondedForce', **nb.attrib)

    def normalized(value):
        return value.removeprefix('protein-')

    def basic_type(amber_type, atomic_number):
        aliases = {'2C':'CT', '3C':'CT', 'CX':'CT', 'C8':'CT', 'CO':'C'}
        symbols = {1:'H',6:'C',7:'N',8:'O',16:'S'}
        for wanted in (amber_type, aliases.get(amber_type, amber_type)):
            matches = [name for name, row in types.items()
                       if normalized(row.get('class','')) == wanted
                       and row.get('element') == symbols[atomic_number]]
            preferred = 'protein-'+wanted
            if preferred in matches:
                return preferred
            if matches:
                return matches[0]
        raise ValueError(f'{force_field}: no supported protein type for {amber_type}')

    def nonbonded(type_name):
        cls = types[type_name].get('class')
        for row in nb:
            if row.tag == 'Atom' and (row.get('type') == type_name or row.get('class') == cls):
                return row
        raise ValueError('Missing native nonbonded definition for '+type_name)

    found = []
    for library in sorted(Path(folder).glob('*.lib')):
        name = library.stem
        if names is not None and name not in names:
            continue
        spec = by_name(name)
        source_name = ('N' if spec['context']=='n_terminal' else 'C' if spec['context']=='c_terminal' else '')+spec['parent']
        if source_name not in templates:
            raise ValueError(f'{force_field}: native parent template {source_name} is unavailable')
        parent = templates[source_name]
        mapping = _REUSED_MAPS.get(spec.get('reused_from'), {})
        inverse = {new: old for old, new in mapping.items()}
        unit = read_units(library)[name]
        residue = E.SubElement(residues, 'Residue', name=name)
        for atom in unit['atoms']:
            role = inverse.get(atom['name'], atom['name'])
            # Keep parent-specific N/CA/C types, including ff19SB CMAP identities.
            carbonyl_changed = ('side' in spec['modified_groups'] and role in ('CG','CD','OD1','OE1'))
            type_name = parent.get(role) if atom['atomic_number'] != 1 and not carbonyl_changed else None
            type_name = type_name or basic_type(atom['amber_type'], atom['atomic_number'])
            if not per_residue_charge:
                unique = f'cation-{name}-{atom["name"]}'
                attributes = dict(types[type_name].attrib)
                attributes['name'] = unique
                E.SubElement(custom_types, 'Type', **attributes)
                parameters = nonbonded(type_name)
                E.SubElement(output_nb, 'Atom', type=unique, charge=str(atom['initial_template_charge_e']),
                             sigma=parameters.get('sigma'), epsilon=parameters.get('epsilon'))
                type_name = unique
            E.SubElement(residue, 'Atom', name=atom['name'], type=type_name,
                         charge=str(atom['initial_template_charge_e']))
        for i,j in unit['bonds']:
            E.SubElement(residue, 'Bond', atomName1=unit['atoms'][i]['name'], atomName2=unit['atoms'][j]['name'])
        if spec['context'] != 'n_terminal':
            E.SubElement(residue, 'ExternalBond', atomName=mapping.get('N','N'))
        if spec['context'] != 'c_terminal':
            E.SubElement(residue, 'ExternalBond', atomName=mapping.get('C','C'))
        if spec['parent'] == 'CYX':
            E.SubElement(residue, 'ExternalBond', atomName='SG')
        found.append(name)
    if names is not None and set(found) != set(names):
        raise ValueError('One or more requested custom templates are unavailable')
    return E.tostring(root, encoding='unicode')


def load_forcefield(folder, force_field='ff14SB', implicit=False, names=None):
    import io
    from openmm.app import ForceField
    if force_field == 'ff14SB':
        from .openmm_backend import forcefield
        return forcefield(folder, implicit, names)
    files = [profile(force_field)['xml']]
    if implicit:
        files.append('implicit/gbn2.xml')
    return ForceField(*files, io.StringIO(custom_xml(folder, force_field, names)))
