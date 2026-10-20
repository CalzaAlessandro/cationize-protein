"""Independent exported-template checks using OpenMM, with no Amber executable."""
import argparse
import json
from pathlib import Path
import numpy as np
from .qm import digest


def capped_structure(mapping, coordinates, charges):
    import parmed as pmd
    from parmed.periodic_table import Element, Mass
    atoms = mapping['atoms']
    xyz, q = np.asarray(coordinates, float), np.asarray(charges, float)
    if xyz.shape != (len(atoms), 3) or q.shape != (len(atoms),):
        raise ValueError('Capped model coordinate/charge dimensions differ')
    if [a['index'] for a in atoms] != list(range(1, len(atoms)+1)):
        raise ValueError('Capped model atom indices must be consecutive')
    if not np.isfinite(xyz).all() or not np.isfinite(q).all():
        raise ValueError('Capped model contains non-finite data')
    order = sorted(range(len(atoms)), key=lambda i: (atoms[i]['residue_index'], i))
    structure, lookup = pmd.Structure(), {}
    for i in order:
        a = atoms[i]
        z = a['atomic_number']
        atom = pmd.Atom(name=a['name'], type=a['amber_type'], atomic_number=z,
                        mass=Mass[Element[z]], charge=float(q[i]))
        atom.xx, atom.xy, atom.xz = xyz[i]
        structure.add_atom(atom, a['residue'], a['residue_index'], chain='A')
        lookup[i+1] = atom
    for (a, b), bond_order in zip(mapping['bonds'], mapping['bond_orders']):
        structure.bonds.append(pmd.Bond(lookup[a], lookup[b], order=bond_order))
    return structure, [i+1 for i in order]


def validate(project, name):
    from .openmm_backend import parameterize
    from .resp import mk_coordinates
    from .qm import completed
    project = Path(project).resolve()
    model = project/'work/models'/name
    fit_root = project/'work/fits'/name
    fit_marker = fit_root/'FIT_COMPLETE.json'
    fit = json.loads(fit_marker.read_text())
    library = project/'parameters'/(name+'.lib')
    if fit['library_sha256'] != digest(library):
        raise ValueError('Exported library hash differs from fitted template')
    if fit['model_sha256'] != digest(model/'atom_mapping.json'):
        raise ValueError('Model changed after fitting')
    attempts = sorted(fit_root.glob('attempt_*/stage2.charges'))
    matching = [p for p in attempts if (p.parent/(name+'.lib')).exists()
                and digest(p.parent/(name+'.lib')) == fit['library_sha256']]
    if len(matching) != 1:
        raise ValueError('Cannot identify unique exported fit attempt')
    mapping = json.loads((model/'atom_mapping.json').read_text())
    natoms = len(mapping['atoms'])
    charges = np.asarray([float(v) for v in matching[0].read_text().split()]).reshape(2, natoms)[0]
    evidence = []
    for label in ('alpha', 'extended'):
        qm = completed(model/label, project/'work/qm'/name/label)
        if qm is None:
            raise ValueError('Validated Gaussian result missing')
        xyz = mk_coordinates((qm/'gaussian.log').read_text(), natoms)
        structure, order = capped_structure(mapping, xyz, charges)
        # Both ordinary ff14SB and GBn2 must accept the exported template.
        for implicit in (False, True):
            typed, system = parameterize(structure, project/'parameters', implicit=implicit, names=[name])
            expected = charges[np.asarray(order)-1]
            error = float(np.max(np.abs(np.asarray([a.charge for a in typed.atoms])-expected)))
            if error > 1e-7 or len(typed.bonds) != len(mapping['bonds']):
                raise ValueError('Exported template changed fitted charges or covalent graph')
            from openmm import Context, VerletIntegrator, Platform, unit
            integrator = VerletIntegrator(.001*unit.picoseconds)
            context = Context(system, integrator, Platform.getPlatformByName('Reference'))
            context.setPositions(typed.coordinates*unit.angstrom)
            state = context.getState(getEnergy=True, getForces=True)
            energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
            forces = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole/unit.nanometer))
            del context, integrator
            if not np.isfinite(energy) or not np.isfinite(forces).all():
                raise ValueError('Exported template gives non-finite OpenMM energy/forces')
            evidence.append(dict(conformer=label, implicit_GBn2=implicit,
                                 charge_max_error_e=error, potential_energy_kJ_mol=energy,
                                 n_atoms=len(typed.atoms), n_bonds=len(typed.bonds)))
    result = dict(name=name, status='RESP_capped_Amber_and_OpenMM_validated',
                  fit_marker_sha256=digest(fit_marker), library_sha256=digest(library),
                  validation_code_sha256=digest(Path(__file__)),
                  backend_code_sha256=digest(Path(__file__).with_name('openmm_backend.py')), evidence=evidence)
    destination = fit_root/'VALIDATED_TEMPLATE.json'
    if destination.exists():
        old = json.loads(destination.read_text())
        if old != result:
            raise ValueError('Existing exported-template validation differs')
        return old
    with destination.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project', type=Path)
    parser.add_argument('name', help='Three-letter custom residue name')
    args = parser.parse_args()
    print(json.dumps(validate(args.project, args.name), indent=2))
