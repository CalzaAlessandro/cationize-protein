"""Byte-preserving residue-name migration; never edit active MD files in place.

The mapping is supplied by the caller. Names must have the same length so the
reverse mapping can restore even unusual fixed-column spacing byte for byte.
Only ATOM/HETATM PDB names and Amber RESIDUE_LABEL fields are eligible.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
from collections import Counter
from pathlib import Path
from collections.abc import Mapping


def _mapping_bytes(mapping: Mapping[str, str], width: int) -> dict[bytes, bytes]:
    result = {}
    for old, new in mapping.items():
        if not isinstance(old, str) or not isinstance(new, str):
            raise ValueError("Residue names must be strings")
        if not re.fullmatch(r"[A-Za-z0-9]{1," + str(width) + r"}", old):
            raise ValueError(f"Invalid source residue name: {old!r}")
        if not re.fullmatch(r"[A-Za-z0-9]{1," + str(width) + r"}", new):
            raise ValueError(f"Invalid target residue name: {new!r}")
        if len(old) != len(new):
            raise ValueError("Name lengths must match for exact reversible migration")
        result[old.encode("ascii")] = new.encode("ascii")
    if len(set(result.values())) != len(result):
        raise ValueError("Residue mapping must be injective")
    return result


def _check_collisions(labels: set[bytes], mapping: dict[bytes, bytes]) -> None:
    collisions = labels & (set(mapping.values()) - set(mapping))
    if collisions:
        raise ValueError("Residue name collision prevents exact inverse: " +
                         ", ".join(sorted(x.decode("ascii") for x in collisions)))


def _replace_field(field: bytes, mapping: dict[bytes, bytes]) -> bytes:
    label = field.strip()
    return field.replace(label, mapping[label], 1) if label in mapping else field


def rename_pdb(data: bytes, mapping: Mapping[str, str]) -> bytes:
    """Rename columns 18-20 of ATOM/HETATM records, preserving every other byte."""
    names = _mapping_bytes(mapping, 3)
    lines = data.splitlines(keepends=True)
    atoms = [i for i, line in enumerate(lines)
             if line[:6] in (b"ATOM  ", b"HETATM")]
    if any(len(lines[i].rstrip(b"\r\n")) < 20 for i in atoms):
        raise ValueError("Truncated PDB atom record")
    _check_collisions({lines[i][17:20].strip() for i in atoms}, names)
    for i in atoms:
        lines[i] = lines[i][:17] + _replace_field(lines[i][17:20], names) + lines[i][20:]
    return b"".join(lines)


def _residue_section(data: bytes) -> tuple[list[bytes], int, int, int]:
    lines = data.splitlines(keepends=True)
    indices = [i for i, line in enumerate(lines)
               if line.rstrip(b"\r\n ") == b"%FLAG RESIDUE_LABEL"]
    if len(indices) != 1:
        raise ValueError("Expected exactly one RESIDUE_LABEL section")
    format_index = indices[0] + 1
    while format_index < len(lines) and lines[format_index].startswith(b"%COMMENT"):
        format_index += 1
    match = (re.fullmatch(rb"%FORMAT\(\s*\d+[aA](\d+)\s*\)\s*", lines[format_index])
             if format_index < len(lines) else None)
    if match is None or int(match[1]) != 4:
        raise ValueError("RESIDUE_LABEL must use the Amber a4 format")
    end = format_index + 1
    while end < len(lines) and not lines[end].startswith(b"%FLAG "):
        end += 1
    return lines, format_index + 1, end, 4


def rename_parm7(data: bytes, mapping: Mapping[str, str]) -> bytes:
    """Rename only the a4 fields in the RESIDUE_LABEL section of an Amber parm7."""
    names = _mapping_bytes(mapping, 4)
    lines, start, end, width = _residue_section(data)
    fields = []
    for i in range(start, end):
        if lines[i].startswith(b"%COMMENT"):
            continue
        body = lines[i].rstrip(b"\r\n")
        if len(body) % width:
            raise ValueError("Incomplete a4 field in RESIDUE_LABEL")
        fields.extend(body[j:j + width] for j in range(0, len(body), width))
    _check_collisions({field.strip() for field in fields}, names)
    for i in range(start, end):
        if lines[i].startswith(b"%COMMENT"):
            continue
        body = lines[i].rstrip(b"\r\n")
        ending = lines[i][len(body):]
        lines[i] = b"".join(_replace_field(body[j:j + width], names)
                             for j in range(0, len(body), width)) + ending
    return b"".join(lines)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parm7_invariants(data: bytes) -> tuple[bytes, dict[str, str], dict[str, int]]:
    lines, start, end, width = _residue_section(data)
    stable = b"".join(lines[:start] + lines[end:])
    labels = Counter()
    for line in lines[start:end]:
        if not line.startswith(b"%COMMENT"):
            body = line.rstrip(b"\r\n")
            labels.update(body[i:i + width].strip().decode("ascii")
                          for i in range(0, len(body), width) if body[i:i + width].strip())
    flags = {}
    section_name = None
    chunks = []
    for line in lines:
        if line.startswith(b"%FLAG "):
            if section_name is not None:
                flags[section_name] = _sha256(b"".join(chunks))
            section_name = line[6:].strip().decode("ascii")
            chunks = []
        chunks.append(line)
    if section_name is not None:
        flags[section_name] = _sha256(b"".join(chunks))
    flags.pop("RESIDUE_LABEL", None)
    return stable, flags, dict(sorted(labels.items()))


def _pdb_invariants(data: bytes) -> tuple[bytes, bytes, bytes, dict[str, int]]:
    atoms = [line for line in data.splitlines(keepends=True)
             if line[:6] in (b"ATOM  ", b"HETATM")]
    stable = b"".join(line[:17] + b"   " + line[20:]
                       if line[:6] in (b"ATOM  ", b"HETATM") else line
                       for line in data.splitlines(keepends=True))
    coordinates = b"".join(line[30:54] for line in atoms)
    atom_order = b"".join(line[:17] + line[20:30] for line in atoms)
    counts = Counter(line[17:20].strip().decode("ascii") for line in atoms)
    return stable, coordinates, atom_order, dict(sorted(counts.items()))


def _validated_copy(data: bytes, mapping: Mapping[str, str], kind: str) -> tuple[bytes, dict]:
    transform = rename_pdb if kind == "pdb" else rename_parm7
    renamed = transform(data, mapping)
    inverse = {new: old for old, new in mapping.items()}
    if transform(renamed, inverse) != data:
        raise ValueError("Inverse transformation does not restore original bytes")
    validation = {"exact_inverse": True}
    if kind == "pdb":
        before = _pdb_invariants(data)
        after = _pdb_invariants(renamed)
        if before[:3] != after[:3]:
            raise ValueError("PDB coordinates, atom order or non-label bytes changed")
        validation.update(all_other_bytes_identical=True, coordinates_identical=True,
                          atom_order_identical=True, non_label_sha256=_sha256(before[0]),
                          coordinates_sha256=_sha256(before[1]), atom_order_sha256=_sha256(before[2]))
        old_labels, new_labels = before[3], after[3]
        label_count_unit = "atom_records"
    else:
        before = _parm7_invariants(data)
        after = _parm7_invariants(renamed)
        if before[:2] != after[:2]:
            raise ValueError("Amber sections outside RESIDUE_LABEL changed")
        validation.update(all_other_bytes_identical=True, non_label_sha256=_sha256(before[0]),
                          unchanged_section_sha256=before[1])
        old_labels, new_labels = before[2], after[2]
        label_count_unit = "residues"
    return renamed, {"source_sha256": _sha256(data), "renamed_sha256": _sha256(renamed),
                     "source_bytes": len(data), "renamed_bytes": len(renamed),
                     "old_label_counts": old_labels, "new_label_counts": new_labels,
                     "label_count_unit": label_count_unit, "validation": validation}


def export_series(source_root: str | Path, output_root: str | Path,
                  mapping: Mapping[str, str],
                  levels: tuple[str, ...] = ("0", "25", "50", "75", "100")) -> dict:
    """Create a new export tree atomically, refusing any existing destination.

    Each source level needs protein.parm7, minimized.pdb and minimized.rst7.
    system.parm7 is included when present. Restart coordinates are referenced by
    relative symlink, never copied or altered. Consumers must treat the shared
    coordinates as read-only. No trajectories, job scripts or MD inputs migrate.
    """
    _mapping_bytes(mapping, 3)
    source_root = Path(source_root).resolve()
    requested_output = Path(output_root).absolute()
    if os.path.lexists(requested_output):
        raise FileExistsError(f"Refusing to overwrite existing export: {requested_output}")
    output_root = requested_output.parent.resolve() / requested_output.name
    levels = tuple(levels)
    if not levels or len(set(levels)) != len(levels) or any(
            not isinstance(level, str) or level in ("", ".", "..") or
            Path(level).name != level or "/" in level or "\\" in level for level in levels):
        raise ValueError("Levels must be unique, nonempty directory names without traversal")
    if os.path.lexists(output_root):
        raise FileExistsError(f"Refusing to overwrite existing export: {output_root}")
    if not source_root.is_dir():
        raise FileNotFoundError(source_root)
    output_root.parent.mkdir(parents=True, exist_ok=True)
    report = {"schema_version": 1, "mapping": dict(mapping),
              "source_root": str(source_root), "export_root": str(output_root),
              "original_inputs_unchanged": True, "cases": {}}
    source_hashes = {}
    with tempfile.TemporaryDirectory(dir=output_root.parent, prefix=".residue-rename-") as tmp:
        staging = Path(tmp) / "exports"
        staging.mkdir()
        for level in levels:
            source = source_root / level
            staged = staging / level
            staged.mkdir()
            final = output_root / level
            files = {}
            for name in ("protein.parm7", "minimized.pdb", "system.parm7"):
                original = source / name
                if name == "system.parm7" and not original.exists():
                    continue
                data = original.read_bytes()
                renamed, record = _validated_copy(data, mapping, "pdb" if name.endswith(".pdb") else "parm7")
                source_hashes[original] = record["source_sha256"]
                (staged / name).write_bytes(renamed)
                if _file_sha256(staged / name) != record["renamed_sha256"]:
                    raise ValueError(f"Export verification failed: {name}")
                record.update(source=str(original), destination=str(final / name))
                files[name] = record
            coordinates = source / "minimized.rst7"
            coordinates_hash = _file_sha256(coordinates)
            source_hashes[coordinates] = coordinates_hash
            relative = os.path.relpath(coordinates, final)
            (staged / "minimized.rst7").symlink_to(relative)
            case = {"files": files, "coordinates": {"source": str(coordinates),
                    "sha256": coordinates_hash, "relative_symlink": relative,
                    "treat_shared_source_as_read_only": True,
                    "contents_copied": False}}
            report["cases"][level] = case
            (staged / "rename_manifest.json").write_text(json.dumps(
                {"schema_version": 1, "mapping": dict(mapping), **case}, indent=2, sort_keys=True) + "\n")
        for path, expected_hash in source_hashes.items():
            if _file_sha256(path) != expected_hash:
                raise ValueError(f"Source changed during migration; no exports published: {path}")
        (staging / "migration_manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        (staging / "README.txt").write_text(
            "Independent residue-name exports\n\n"
            "Only PDB ATOM/HETATM residue-name fields and Amber RESIDUE_LABEL fields differ.\n"
            "Source structures, force-field values, MD inputs and trajectories are unchanged.\n"
            "minimized.rst7 is a relative symlink to the original coordinates: treat it as read-only.\n"
            "If system.parm7 is present, pair it with the original solvated restart/trajectory,\n"
            "not with the dry-protein minimized.rst7. No solvated coordinates were copied.\n"
            "Manifests record source paths, old/new label counts, hashes, exact inverse checks,\n"
            "all untouched Amber-section hashes and PDB coordinate/atom-order checks.\n")
        if os.path.lexists(output_root):
            raise FileExistsError(f"Export destination appeared during migration: {output_root}")
        staging.rename(output_root)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--map", action="append", required=True, metavar="OLD=NEW",
                        help="Equal-length residue names; repeat for each change")
    parser.add_argument("--levels", nargs="+", default=["0", "25", "50", "75", "100"])
    args = parser.parse_args(argv)
    mapping = {}
    for item in args.map:
        if item.count("=") != 1:
            parser.error("Each mapping must be OLD=NEW")
        old, new = item.split("=")
        if old in mapping:
            parser.error(f"Duplicate source residue in mapping: {old}")
        mapping[old] = new
    report = export_series(args.source, args.output, mapping, levels=tuple(args.levels))
    print(json.dumps({"export_root": report["export_root"], "levels": list(report["cases"]),
                      "original_inputs_unchanged": report["original_inputs_unchanged"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
