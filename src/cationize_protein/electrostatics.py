"""Portable APBS sampling using caller-supplied topology charges and radii.

Distances are Å, charges are elementary charges, and APBS potentials are kBT/e.
No protonation, charges, radii, or residue identities are inferred here. APBS is
an external runtime dependency supplied explicitly through :class:`APBSConfig`.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Mapping, Sequence

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree


class APBSError(RuntimeError):
    """APBS execution or convergence validation failed; logs remain available."""


@dataclass(frozen=True)
class APBSConfig:
    """Explicit solver and PB settings; a symmetric monovalent salt is used."""

    executable: str | Path
    temperature_K: float = 298.0
    solute_dielectric: float = 2.0
    solvent_dielectric: float = 78.5
    salt_M: float = 0.15
    ion_radius_A: float = 2.0
    spacing_A: float = 0.75
    padding_A: float = 20.0
    coarse_factor: float = 1.6
    max_grid_points: int = 40_000_000
    timeout_seconds: float = 3600.0

    def __post_init__(self):
        if not str(self.executable).strip():
            raise ValueError('An APBS executable is required')
        for key in ('temperature_K', 'solute_dielectric', 'solvent_dielectric',
                    'ion_radius_A', 'spacing_A', 'padding_A', 'timeout_seconds'):
            value = getattr(self, key)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{key} must be finite and positive')
        if not np.isfinite(self.salt_M) or self.salt_M < 0:
            raise ValueError('salt_M must be finite and nonnegative')
        if not np.isfinite(self.coarse_factor) or self.coarse_factor < 1:
            raise ValueError('coarse_factor must be at least one')
        if not isinstance(self.max_grid_points, int) or self.max_grid_points < 8:
            raise ValueError('max_grid_points must be an integer of at least eight')


@dataclass(frozen=True)
class CandidateGrid:
    """Union of local spheres, per-candidate point indices, and ion exclusions."""

    points: np.ndarray
    indices: dict[str, np.ndarray]
    blocked: np.ndarray


@dataclass(frozen=True)
class CandidateRank:
    candidate_id: str
    score_kBT_per_e: float
    point_index: int
    accessible_grid_points: int


def _coordinates(values, name: str, *, allow_empty: bool = False) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if allow_empty and array.size == 0:
        return np.empty((0, 3), dtype=float)
    if array.ndim != 2 or array.shape[1] != 3 or not np.isfinite(array).all():
        raise ValueError(f'{name} must be a finite N by 3 array')
    if not allow_empty and not len(array):
        raise ValueError(f'{name} cannot be empty')
    return array


def _vector(values, size: int, name: str, *, nonnegative: bool = False) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (size,) or not np.isfinite(array).all():
        raise ValueError(f'{name} must contain {size} finite values')
    if nonnegative and np.any(array < 0):
        raise ValueError(f'{name} cannot contain negative values')
    return array


def build_candidate_grid(
    candidate_centers: Mapping[str, Sequence[Sequence[float]]],
    coordinates, lj_radii_A, ion_radius_A: float, *, spacing_A: float = 1.0,
    sphere_radius_A: float = 5.0, max_points: int = 5_000_000,
) -> CandidateGrid:
    """Sample the union of local spheres on a global lattice, blocking ion clashes.

    A site may have several oxygen centers. IDs are opaque strings (for example
    ``chain:residue:group``); separate groups sharing a residue remain separate.
    A point is blocked when its distance to any atom is no greater than the
    supplied LJ radius plus the ion radius. No PB radius is substituted for LJ.
    """
    xyz = _coordinates(coordinates, 'coordinates', allow_empty=True)
    radii = _vector(lj_radii_A, len(xyz), 'lj_radii_A', nonnegative=True)
    if not np.isfinite(ion_radius_A) or ion_radius_A < 0:
        raise ValueError('ion_radius_A must be finite and nonnegative')
    if any(not np.isfinite(x) or x <= 0 for x in (spacing_A, sphere_radius_A)):
        raise ValueError('spacing and sphere radius must be finite and positive')
    if not isinstance(max_points, int) or max_points < 1:
        raise ValueError('max_points must be a positive integer')
    if any(not isinstance(key, str) or not key for key in candidate_centers):
        raise ValueError('Candidate IDs must be nonempty strings')
    pieces, keys = [], []
    total = 0
    for candidate_id in sorted(candidate_centers):
        centers = np.asarray(candidate_centers[candidate_id], dtype=float)
        if centers.shape == (3,):
            centers = centers.reshape(1, 3)
        centers = _coordinates(centers, 'candidate centers')
        spheres = []
        for center in centers:
            lower_float = np.ceil((center - sphere_radius_A) / spacing_A)
            upper_float = np.floor((center + sphere_radius_A) / spacing_A)
            # Reject extreme lattice indices before an integer cast can wrap.
            if np.max(np.abs(np.r_[lower_float, upper_float])) > 2**52:
                raise ValueError('Grid indices exceed exact integer range')
            lower, upper = lower_float.astype(np.int64), upper_float.astype(np.int64)
            extent = upper - lower + 1
            if int(extent[0]) * int(extent[1]) * int(extent[2]) > max_points:
                raise ValueError('Candidate sphere exceeds max_points')
            lattice = np.stack(np.meshgrid(*(np.arange(a, b + 1) for a, b in zip(lower, upper)), indexing='ij'), axis=-1).reshape(-1, 3)
            keep = np.linalg.norm(lattice * spacing_A - center, axis=1) <= sphere_radius_A + 1e-10
            spheres.append(lattice[keep])
        local = np.unique(np.concatenate(spheres), axis=0)
        total += len(local)
        if total > max_points:
            raise ValueError('Candidate grid exceeds max_points')
        pieces.append(local)
        keys.append(candidate_id)
    if not pieces:
        return CandidateGrid(np.empty((0, 3)), {}, np.empty(0, dtype=bool))
    lattice, inverse = np.unique(np.concatenate(pieces), axis=0, return_inverse=True)
    points = lattice * spacing_A
    blocked = np.zeros(len(points), dtype=bool)
    if len(points):
        tree = cKDTree(points)
        for hits in tree.query_ball_point(xyz, radii + ion_radius_A):
            blocked[hits] = True
    indices, offset = {}, 0
    for key, local in zip(keys, pieces):
        indices[key] = inverse[offset:offset + len(local)]
        offset += len(local)
    return CandidateGrid(points, indices, blocked)


def rank_candidates(scores, indices: Mapping[str, Sequence[int]], blocked):
    """Return accessible local minima ordered by potential, then stable site ID.

    Tied points select the smallest global point index. Unavailable IDs are
    returned separately; a fully blocked sphere never acquires a score.
    """
    scores = np.asarray(scores, dtype=float)
    mask = np.asarray(blocked, dtype=bool)
    if scores.ndim != 1 or mask.shape != scores.shape or not np.isfinite(scores).all():
        raise ValueError('scores and blocked must be matching finite vectors')
    if any(not isinstance(key, str) or not key for key in indices):
        raise ValueError('Candidate IDs must be nonempty strings')
    rows, unavailable = [], []
    for candidate_id in sorted(indices):
        raw = np.asarray(indices[candidate_id])
        if raw.ndim != 1 or (raw.size and not np.issubdtype(raw.dtype, np.integer)):
            raise ValueError('Candidate indices must be integer vectors')
        points = np.unique(raw.astype(np.int64))
        if np.any(points < 0) or np.any(points >= len(scores)):
            raise ValueError('Candidate index outside score vector')
        valid = points[~mask[points]]
        if not len(valid):
            unavailable.append(candidate_id)
            continue
        best = int(valid[np.argmin(scores[valid])])
        rows.append(CandidateRank(candidate_id, float(scores[best]), best, len(valid)))
    return sorted(rows, key=lambda row: (row.score_kBT_per_e, row.candidate_id)), unavailable


@dataclass(frozen=True)
class DXGrid:
    """Scalar OpenDX grid with mutually orthogonal, possibly rotated axes."""

    values: np.ndarray
    origin: np.ndarray
    deltas: np.ndarray

    def __post_init__(self):
        values = np.asarray(self.values, dtype=float)
        origin = np.asarray(self.origin, dtype=float)
        deltas = np.asarray(self.deltas, dtype=float)
        if values.ndim != 3 or any(n < 1 for n in values.shape) or not np.isfinite(values).all():
            raise ValueError('DX values must be a nonempty finite three dimensional array')
        if origin.shape != (3,) or deltas.shape != (3, 3) or not np.isfinite(origin).all() or not np.isfinite(deltas).all():
            raise ValueError('Invalid DX origin or deltas')
        lengths = np.linalg.norm(deltas, axis=1)
        if np.any(lengths == 0) or not np.isfinite(lengths).all():
            raise ValueError('DX axes must have nonzero finite lengths')
        normalized = deltas / lengths[:, None]
        if not np.allclose(normalized @ normalized.T, np.eye(3), atol=1e-8, rtol=0):
            raise ValueError('Only orthogonal DX grids are supported')
        object.__setattr__(self, 'values', values)
        object.__setattr__(self, 'origin', origin)
        object.__setattr__(self, 'deltas', deltas)

    def sample(self, points) -> np.ndarray:
        """Trilinear interpolation in lattice coordinates; outside points fail."""
        points = _coordinates(points, 'sample points', allow_empty=True)
        if not len(points):
            return np.empty(0)
        lattice = np.linalg.solve(self.deltas.T, (points - self.origin).T).T
        maximum = np.asarray(self.values.shape) - 1
        if np.any(lattice < -1e-9) or np.any(lattice > maximum + 1e-9):
            raise ValueError('Sample point outside DX grid')
        lattice = np.clip(lattice, 0, maximum)
        axes = tuple(np.arange(size, dtype=float) for size in self.values.shape)
        result = RegularGridInterpolator(axes, self.values, bounds_error=True)(lattice)
        if not np.isfinite(result).all():
            raise ValueError('Nonfinite interpolated potential')
        return result


def read_dx(path: str | Path, *, max_items: int = 40_000_000) -> DXGrid:
    """Read bounded ASCII scalar OpenDX data, rejecting malformed/skew fields."""
    if not isinstance(max_items, int) or max_items < 1:
        raise ValueError('max_items must be a positive integer')
    counts, origin, deltas, items = None, None, [], None
    values = []
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            if 'class gridpositions counts' in line:
                match = re.search(r'counts\s+(\d+)\s+(\d+)\s+(\d+)\s*$', line)
                if not match:
                    raise ValueError('Malformed DX counts')
                counts = tuple(int(value) for value in match.groups())
                if min(counts) < 1 or counts[0] * counts[1] * counts[2] > max_items:
                    raise ValueError('DX grid exceeds max_items or has empty axes')
            elif line.startswith('origin '):
                try:
                    origin = np.array([float(value) for value in line.split()[1:]])
                except ValueError as exc:
                    raise ValueError('Malformed DX origin') from exc
            elif line.startswith('delta '):
                try:
                    deltas.append([float(value) for value in line.split()[1:]])
                except ValueError as exc:
                    raise ValueError('Malformed DX delta') from exc
            elif 'data follows' in line:
                match = re.search(r'\bitems\s+(\d+)\s+data follows\s*(.*)$', line)
                if not match or not re.search(r'\brank\s+0\b', line) or 'binary' in line:
                    raise ValueError('Expected ASCII scalar DX data')
                items = int(match.group(1))
                if counts is None or items != counts[0] * counts[1] * counts[2] or items > max_items:
                    raise ValueError('DX item count differs from grid dimensions')
                pending = match.group(2).split()
                try:
                    values.extend(float(value) for value in pending)
                    while len(values) < items:
                        data_line = next(handle).split('#', 1)[0]
                        values.extend(float(value) for value in data_line.split())
                except (ValueError, StopIteration) as exc:
                    raise ValueError('Truncated or invalid DX values') from exc
                if len(values) != items:
                    raise ValueError('DX contains excess values')
                break
    if counts is None or origin is None or len(deltas) != 3 or items is None:
        raise ValueError('Incomplete DX header')
    return DXGrid(np.asarray(values).reshape(counts), origin, np.asarray(deltas))


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _sha_array(array) -> str:
    array = np.ascontiguousarray(array, dtype='<f8')
    return hashlib.sha256(str(array.shape).encode() + array.tobytes()).hexdigest()


def _labels(values, size, default, name):
    labels = [default] * size if values is None else list(values)
    if len(labels) != size or any(not re.fullmatch(r'[A-Za-z0-9_+*\-]+', str(value)) for value in labels):
        raise ValueError(f'{name} must contain one whitespace-free token per atom')
    return labels


def _check_apbs_log(log: str):
    failures = (
        r'\bfailed\s+to\s+converg', r'\bnot\s+converg',
        r'\bconvergence\s+(?:failure|failed|not\s+achieved)',
        r'\berror\s+tolerance\s+not\b', r'\b(?:fatal|error)\s*[:!]',
        r'^\s*(?:[A-Za-z]\w*\s*:\s*)?(?:error(?!\s+tolerance\b)|fatal)\b',
        r'\bunconverged\b', r'\bsegmentation\s+fault\b',
        r'\bmaximum\s+(?:number\s+of\s+)?iterations?\s+(?:reached|exceeded)',
        r'\b(?:VASSERT|assertion\s+failed)\b',
    )
    for pattern in failures:
        if re.search(pattern, log, flags=re.IGNORECASE | re.MULTILINE):
            raise APBSError('APBS reported a solver error or nonconvergence; inspect apbs.log')
    if not re.search(r'\bconverged\b|thanks\s+for\s+using\s+APBS|total\s+execution\s+time', log, re.IGNORECASE):
        raise APBSError('APBS did not report normal completion; inspect apbs.log')


def run_apbs(
    coordinates, charges_e, pb_radii_A, sample_points, output_dir: str | Path,
    config: APBSConfig, *, atom_names=None, residue_names=None,
    residue_numbers=None, provenance: Mapping | None = None, keep_map: bool = False,
):
    """Solve the current supplied charge distribution and interpolate at points.

    The output directory must be new. PQR, APBS input, and logs are preserved;
    ``PB_COMPLETE.json`` is created only after process, log, DX, and interpolation
    checks pass. The DX map is retained only when requested. Caller provenance
    is nested under ``inputs`` and cannot overwrite computed checksums/settings.
    A normal APBS completion marker plus absence of known solver failures is
    required; this is a log audit, not an independent numerical convergence test.
    """
    xyz = _coordinates(coordinates, 'coordinates')
    charges = _vector(charges_e, len(xyz), 'charges_e')
    radii = _vector(pb_radii_A, len(xyz), 'pb_radii_A', nonnegative=True)
    points = _coordinates(sample_points, 'sample_points', allow_empty=True)
    names = _labels(atom_names, len(xyz), 'X', 'atom_names')
    residues = _labels(residue_names, len(xyz), 'MOL', 'residue_names')
    numbers = np.ones(len(xyz), dtype=int) if residue_numbers is None else np.asarray(residue_numbers)
    if numbers.shape != (len(xyz),) or not np.issubdtype(numbers.dtype, np.integer):
        raise ValueError('residue_numbers must contain one integer per atom')
    inputs = dict(provenance or {})
    json.dumps(inputs, allow_nan=False)
    length = np.ptp(xyz, axis=0) + 2 * config.padding_A
    center = xyz.min(axis=0) + np.ptp(xyz, axis=0) / 2
    dimensions_float = 32 * np.ceil(length / (32 * config.spacing_A)) + 1
    if not np.isfinite(dimensions_float).all() or np.prod(dimensions_float) > config.max_grid_points:
        raise ValueError('APBS dimensions exceed max_grid_points')
    dimensions = dimensions_float.astype(int)
    if len(points) and np.any(np.abs(points - center) > length / 2 + 1e-9):
        raise ValueError('Sample points lie outside the configured fine grid')
    executable = shutil.which(str(config.executable))
    if executable is None:
        raise APBSError(f'APBS executable unavailable: {config.executable}')
    executable = Path(executable).resolve()
    solver_sha = _sha_file(executable)
    requested_out = Path(output_dir)
    if requested_out.exists() or requested_out.is_symlink():
        raise FileExistsError(f'APBS output directory already exists: {requested_out}')
    out = requested_out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    pqr = out / 'molecule.pqr'
    with pqr.open('x') as handle:
        for index, (name, residue, number, coord, charge, radius) in enumerate(zip(names, residues, numbers, xyz, charges, radii), 1):
            handle.write(f'ATOM {index} {name} {residue} {number} {coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f} {charge:.12f} {radius:.10f}\n')
    vector = lambda values: ' '.join(format(value, '.12g') for value in values)
    text = f'''read
 mol pqr molecule.pqr
end
elec name protein
 mg-auto
 dime {vector(dimensions)}
 cglen {vector(length * config.coarse_factor)}
 fglen {vector(length)}
 cgcent {vector(center)}
 fgcent {vector(center)}
 mol 1
 npbe
 bcfl sdh
 pdie {config.solute_dielectric:.12g}
 sdie {config.solvent_dielectric:.12g}
 ion charge 1 conc {config.salt_M:.12g} radius {config.ion_radius_A:.12g}
 ion charge -1 conc {config.salt_M:.12g} radius {config.ion_radius_A:.12g}
 chgm spl2
 srfm smol
 srad 1.4
 swin 0.3
 sdens 10
 temp {config.temperature_K:.12g}
 calcenergy no
 calcforce no
 write pot dx potential
end
quit
'''
    (out / 'apbs.in').write_text(text)
    started = time.monotonic()
    from .dependencies import external_program_environment
    try:
        with (out / 'apbs.log').open('x') as log_handle:
            result = subprocess.run([str(executable), 'apbs.in'], cwd=out,
                                    stdout=log_handle, stderr=subprocess.STDOUT,
                                    timeout=config.timeout_seconds, check=False,
                                    env=external_program_environment(executable))
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise APBSError('APBS execution failed; inspect apbs.log') from exc
    if result.returncode != 0:
        raise APBSError(f'APBS exited with code {result.returncode}; inspect apbs.log')
    _check_apbs_log((out / 'apbs.log').read_text(errors='replace'))
    map_path = out / 'potential.dx'
    if not map_path.is_file():
        raise APBSError('APBS did not produce potential.dx')
    try:
        field = read_dx(map_path, max_items=config.max_grid_points)
        scores = field.sample(points)
    except ValueError as exc:
        raise APBSError(f'Invalid APBS potential field: {exc}') from exc
    settings = asdict(config)
    settings['executable'] = str(config.executable)
    metadata = dict(
        solver_executable=str(executable), solver_sha256=solver_sha,
        inputs=inputs, config=settings, potential_units='kBT/e',
        coordinates_sha256=_sha_array(xyz), charges_sha256=_sha_array(charges),
        pb_radii_sha256=_sha_array(radii), sample_points_sha256=_sha_array(points),
        pqr_sha256=_sha_file(pqr), apbs_input_sha256=_sha_file(out / 'apbs.in'),
        apbs_log_sha256=_sha_file(out / 'apbs.log'), DX_sha256=_sha_file(map_path),
        grid_points=dimensions.tolist(), grid_spacing_A=(length / (dimensions - 1)).tolist(),
        map_grid_points=list(field.values.shape), map_origin_A=field.origin.tolist(),
        map_deltas_A=field.deltas.tolist(),
        fine_lengths_A=length.tolist(), protein_charge_e=float(charges.sum()),
        wall_seconds=time.monotonic() - started, map_retained=bool(keep_map),
        convergence_check='Normal completion marker and no recognized failure in APBS log',
    )
    if not keep_map:
        map_path.unlink()
    marker = out / '.PB_COMPLETE.json.tmp'
    marker.write_text(json.dumps(metadata, indent=2, allow_nan=False) + '\n')
    marker.replace(out / 'PB_COMPLETE.json')
    return scores, metadata
