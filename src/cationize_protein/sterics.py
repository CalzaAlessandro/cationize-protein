"""Deterministic local bond rotations constrained by an explicit covalent graph.

All indices are zero based; coordinates and LJ half-Rmin radii are in Å. This
module ranks a geometric repulsive overlap heuristic, never a physical energy.
The caller declares chemically rotatable bonds and protects the backbone; no
bond orders, aromaticity, or chemical rotatability are guessed from distances.
"""
from __future__ import annotations

from dataclasses import dataclass
import numbers
from typing import Iterable, Sequence

import numpy as np
from scipy.spatial.distance import cdist


@dataclass(frozen=True)
class RotatableBond:
    """Oriented bond: base stays fixed and the tip-side component can rotate."""

    base: int
    tip: int


@dataclass(frozen=True)
class Rotation:
    base: int
    tip: int
    moving_atoms: tuple[int, ...]


def _index(value, atom_count: int) -> int:
    if not isinstance(value, numbers.Integral) or not 0 <= value < atom_count:
        raise ValueError(f'Atom index {value!r} outside a graph of {atom_count} atoms')
    return int(value)


def _graph(bonds, atom_count: int):
    graph = [set() for _ in range(atom_count)]
    for bond in bonds:
        if len(bond) != 2:
            raise ValueError('Every bond must contain two atom indices')
        i, j = (_index(value, atom_count) for value in bond)
        if i == j:
            raise ValueError('Self-bonds are invalid')
        graph[i].add(j)
        graph[j].add(i)
    return graph


def derive_rotations(
    bonds: Iterable[tuple[int, int]], rotatable_bonds: Sequence[RotatableBond],
    residue_ids: Sequence[str], *, protected_atoms: Iterable[int] = (),
    disulfide_bonds: Iterable[tuple[int, int]] = (),
) -> tuple[Rotation, ...]:
    """Derive distal components, refusing ring/protected/interresidue movement.

    ``residue_ids`` must uniquely identify complete residues across chains.
    Caller-declared disulfide endpoints may never move. An unlisted disulfide
    connecting different residues is also rejected by the interresidue check.
    Any protected atom in a distal component causes rejection, even if that
    atom lies on the rotation axis. Supply every backbone atom as protected.
    """
    count = len(residue_ids)
    if not count or any(not isinstance(value, str) or not value for value in residue_ids):
        raise ValueError('residue_ids must contain a nonempty ID for every atom')
    graph = _graph(bonds, count)
    protected = {_index(value, count) for value in protected_atoms}
    for pair in disulfide_bonds:
        if len(pair) != 2:
            raise ValueError('Every disulfide bond must have two endpoints')
        i, j = (_index(value, count) for value in pair)
        if j not in graph[i]:
            raise ValueError('Declared disulfide is not present in the bond graph')
        protected.update((i, j))
    axes, seen = [], set()
    for bond in rotatable_bonds:
        if not isinstance(bond, RotatableBond):
            raise ValueError('rotatable_bonds must contain RotatableBond objects')
        base, tip = _index(bond.base, count), _index(bond.tip, count)
        if tip not in graph[base]:
            raise ValueError('Rotatable bond is absent from the covalent graph')
        edge = tuple(sorted((base, tip)))
        if edge in seen:
            raise ValueError('A rotatable bond may only be declared once')
        seen.add(edge)
        distal, stack = {tip}, [tip]
        while stack:
            node = stack.pop()
            for neighbor in graph[node]:
                if {node, neighbor} == {base, tip} or neighbor in distal:
                    continue
                distal.add(neighbor)
                stack.append(neighbor)
        if base in distal:
            raise ValueError('Rotation would move a ring')
        if residue_ids[base] != residue_ids[tip] or any(residue_ids[index] != residue_ids[tip] for index in distal):
            raise ValueError('Rotation would move an interresidue covalent component')
        if distal & protected:
            raise ValueError('Rotation would move a protected atom or disulfide')
        axes.append(Rotation(base, tip, tuple(sorted(distal))))
    if not axes:
        raise ValueError('At least one rotatable bond is required')
    return tuple(axes)


def _arrays(coordinates, lj_radii_A):
    xyz = np.asarray(coordinates, dtype=float)
    radii = np.asarray(lj_radii_A, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError('Coordinates must be a finite nonempty N by 3 array')
    if radii.shape != (len(xyz),) or not np.isfinite(radii).all() or np.any(radii < 0):
        raise ValueError('LJ radii must be finite nonnegative values for every atom')
    return xyz, radii


def _exclusions(graph, moving):
    excluded = np.zeros((len(moving), len(graph)), dtype=bool)
    for row, index in enumerate(moving):
        seen, front = {int(index)}, {int(index)}
        for _ in range(3):
            front = {neighbor for node in front for neighbor in graph[node]} - seen
            seen.update(front)
        excluded[row, list(seen)] = True
    return excluded


def _score(coordinates, moving, threshold, excluded):
    if not len(moving):
        return 0.0
    distances = cdist(coordinates[moving], coordinates)
    with np.errstate(over='ignore', invalid='ignore'):
        ratio = np.maximum(threshold / np.maximum(distances, 0.25), 1.0)
        penalty = ratio**12 - 1.0
        penalty[excluded] = 0.0
        result = float(penalty.sum())
    if not np.isfinite(result):
        raise ValueError('Overlap score overflow; inspect LJ radii and coordinates')
    return result


def overlap_score(coordinates, bonds, lj_radii_A, moving_atoms, *, scale: float = 0.75) -> float:
    """Sum repulsive overlap of moving atoms against the entire supplied system.

    The threshold is ``scale * (LJ radius i + LJ radius j)``. Distances below
    0.25 Å are clamped; self and pairs separated by at most three covalent bonds
    are excluded. Moving-moving interactions are counted from both endpoints,
    preserving the prior deterministic heuristic. This is not a force-field
    energy and does not test torsional preferences or electrostatic suitability.
    """
    xyz, radii = _arrays(coordinates, lj_radii_A)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('scale must be finite and positive')
    graph = _graph(bonds, len(xyz))
    moving = np.array(sorted({_index(value, len(xyz)) for value in moving_atoms}), dtype=int)
    threshold = scale * (radii[moving, None] + radii[None, :])
    return _score(xyz, moving, threshold, _exclusions(graph, moving))


def _rotate(coordinates, rotation: Rotation, angle: int):
    if angle == 0:
        return coordinates.copy()
    axis = coordinates[rotation.tip] - coordinates[rotation.base]
    length = np.linalg.norm(axis)
    if not np.isfinite(length) or length <= 1e-12:
        raise ValueError('A rotatable bond has a degenerate geometric axis')
    axis = axis / length
    moving = np.array(rotation.moving_atoms, dtype=int)
    vectors = coordinates[moving] - coordinates[rotation.base]
    theta = np.deg2rad(angle)
    result = coordinates.copy()
    result[moving] = (coordinates[rotation.base] + vectors * np.cos(theta)
                      + np.cross(axis, vectors) * np.sin(theta)
                      + np.outer(vectors @ axis, axis) * (1 - np.cos(theta)))
    # Axis endpoints are mathematically invariant; retain exact supplied values.
    result[rotation.tip] = coordinates[rotation.tip]
    return result


def choose_rotamer(
    coordinates, bonds, lj_radii_A, rotatable_bonds: Sequence[RotatableBond], *,
    residue_ids: Sequence[str], protected_atoms: Iterable[int] = (),
    disulfide_bonds: Iterable[tuple[int, int]] = (), beam_width: int = 8,
    scale: float = 0.75,
):
    """Return copied coordinates and a deterministic 60° beam-search report.

    The explicit bond order determines the sequence of rotations. Stable angle
    tuple ordering resolves tied scores. Every covalent length is verified and
    protected coordinates remain unchanged. The original pose is retained when
    necessary to guarantee a nonincreasing score. Perform physical minimization
    separately; a low geometric overlap score is not structural validation.
    """
    xyz, radii = _arrays(coordinates, lj_radii_A)
    if len(residue_ids) != len(xyz):
        raise ValueError('residue_ids must identify every coordinate')
    if not isinstance(beam_width, numbers.Integral) or beam_width < 1:
        raise ValueError('beam_width must be a positive integer')
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('scale must be finite and positive')
    bonds = tuple(tuple(bond) for bond in bonds)
    protected = tuple(protected_atoms)
    disulfides = tuple(tuple(pair) for pair in disulfide_bonds)
    axes = derive_rotations(bonds, rotatable_bonds, residue_ids,
                            protected_atoms=protected, disulfide_bonds=disulfides)
    graph = _graph(bonds, len(xyz))
    moving = np.array(sorted({index for axis in axes for index in axis.moving_atoms}), dtype=int)
    excluded = _exclusions(graph, moving)
    threshold = scale * (radii[moving, None] + radii[None, :])
    for rotation in axes:
        if np.linalg.norm(xyz[rotation.tip] - xyz[rotation.base]) <= 1e-12:
            raise ValueError('A rotatable bond has a degenerate geometric axis')
    initial = _score(xyz, moving, threshold, excluded)
    beam = [(initial, (), xyz.copy())]
    for rotation in axes:
        candidates = []
        for _, angles, old in beam:
            for angle in (0, 60, 120, 180, 240, 300):
                new = _rotate(old, rotation, angle)
                candidates.append((_score(new, moving, threshold, excluded), angles + (angle,), new))
        beam = sorted(candidates, key=lambda item: (item[0], item[1]))[:beam_width]
    final, angles, result = min(beam, key=lambda item: (item[0], item[1]))
    if final > initial:
        final, angles, result = initial, (0,) * len(axes), xyz.copy()
    for i, j in bonds:
        if not np.isclose(np.linalg.norm(result[i] - result[j]), np.linalg.norm(xyz[i] - xyz[j]), rtol=1e-10, atol=1e-8):
            raise RuntimeError('A rotation changed a covalent bond length')
    fixed = set(protected) | {index for pair in disulfides for index in pair}
    if fixed and not np.array_equal(result[sorted(fixed)], xyz[sorted(fixed)]):
        raise RuntimeError('A rotation moved protected atoms')
    report = dict(
        initial_steric_score=initial, final_steric_score=final,
        rotatable_bonds=[[rotation.base, rotation.tip] for rotation in axes],
        rotations_degrees=list(angles), moving_atoms=moving.tolist(),
        maximum_displacement_A=float(np.max(np.linalg.norm(result - xyz, axis=1))),
        criterion='Geometric overlap using 0.75 times LJ half-Rmin sums by default; excludes pairs within three covalent bonds; not a physical energy',
        overlap_scale=float(scale), beam_width=int(beam_width),
    )
    return result, report
