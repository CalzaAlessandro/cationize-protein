#!/usr/bin/env python3
"""Run a minimal, self-contained cationize-protein preflight check and progressive cationization.

Uses OpenMM and analytical Coulomb electrostatics; does not require APBS,
Amber, or GROMACS installations.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run minimal cationize-protein quickstart check and progressive run."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("work_quickstart"),
        help="Target output directory (default: work_quickstart)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pdb_path = Path("examples/pdb_checks/reference_peptide.pdb").resolve()
    output_dir = args.output_dir.resolve()

    if not pdb_path.exists():
        print(f"Error: input PDB not found: {pdb_path}", file=sys.stderr)
        return 1

    # Locate cationize-protein executable
    cationize_cmd = [sys.executable, "-c",
                     "from cationize_protein.cli import main; raise SystemExit(main())"]

    print(f"[Quickstart] Preflight checking {pdb_path.name}...")
    check_cmd = [*cationize_cmd, "check", str(pdb_path)]
    res_check = subprocess.run(check_cmd, capture_output=True, text=True)
    if res_check.returncode != 0:
        print(f"PDB check failed:\n{res_check.stderr}", file=sys.stderr)
        return res_check.returncode
    print(f"[Quickstart] {res_check.stdout.strip()}")

    print(f"[Quickstart] Running progressive cationization with five checkpoints...")
    run_cmd = [
        *cationize_cmd,
        "run",
        "--pdb",
        str(pdb_path),
        "--output",
        str(output_dir),
        "--engine",
        "openmm",
        "--solver",
        "coulomb",
        "--force-field",
        "ff14SB",
        "--formats",
        "amber",
        "openmm",
        "gromacs",
    ]
    res_run = subprocess.run(run_cmd, capture_output=True, text=True)
    if res_run.returncode != 0:
        print(f"Cationization run failed:\n{res_run.stderr}", file=sys.stderr)
        return res_run.returncode

    # Verify stage checkpoints
    expected_stages = ["0", "25", "50", "75", "100"]
    for stage in expected_stages:
        stage_dir = output_dir / stage
        assert stage_dir.is_dir(), f"Missing checkpoint directory for stage {stage}"

    # Verify atom counts
    pdb_0 = output_dir / "0" / "protein.pdb"
    pdb_100 = output_dir / "100" / "protein.pdb"
    atoms_0 = count_pdb_atoms(pdb_0)
    atoms_100 = count_pdb_atoms(pdb_100)

    assert atoms_0 == 46, f"Expected 46 baseline atoms, got {atoms_0}"
    assert atoms_100 == 79, f"Expected 79 stage 100 atoms, got {atoms_100}"
    print(f"[Quickstart] Baseline atoms: {atoms_0} -> Final atoms: {atoms_100}")

    # Verify exports in stage 100
    exp_100 = output_dir / "100" / "exports"
    expected_exports = [
        "system.xml",
        "system.parm7",
        "system.rst7",
        "system.top",
        "system.gro",
        "protein.itp",
        "protein.lib",
        "parameters.frcmod",
        "forcefield.itp",
    ]
    for exp_file in expected_exports:
        assert (exp_100 / exp_file).exists(), f"Missing exported topology: {exp_file}"

    print(f"[Quickstart] Verified topologies in {exp_100.name}/: OpenMM, Amber, GROMACS")
    print("[Quickstart] Success: complete run verified in", output_dir)
    return 0


def count_pdb_atoms(path: Path) -> int:
    with open(path, "r", encoding="utf-8") as f:
        return sum(1 for line in f if line.startswith(("ATOM  ", "HETATM")))



if __name__ == "__main__":
    raise SystemExit(main())
