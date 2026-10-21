"""Standard-library-only dependency discovery and private environment bootstrap.

Discovery describes capabilities; it never substitutes a different physical model.
Only the installer mutates files, and only in a dedicated environment/private tools.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import shlex
import stat
import subprocess
import sys
import tempfile
from urllib.parse import urlparse
from urllib.request import Request, urlopen
import zipfile


class DependencyError(RuntimeError):
    """A dependency cannot safely be installed or used as requested."""


PYTHON_PACKAGES = {
    "numpy": ("numpy", "1.26"),
    "scipy": ("scipy", "1.11"),
    "ParmEd": ("parmed", "4.3"),
    "openmm": ("openmm", "8.6"),
    "propka": ("propka", "3.5"),
    "pdbfixer": ("pdbfixer", "1.12"),
}
AMBER_PROGRAMS = ("pmemd.cuda", "pmemd", "sander", "tleap", "pbsa", "antechamber", "resp")
GROMACS_PROGRAMS = ("gmx", "gmx_mpi")
ENV_MARKER = ".cationize-protein-environment.json"
FORCE_FIELDS = {
    "ff14SB": {
        "supported": True,
        "description": "Amber ff14SB; custom RESP template validation is checked separately",
    },
    "ff19SB": {
        "supported": False,
        "description": "Amber ff19SB; requires parent/custom-template compatibility and CMAP validation",
    },
    "ff99SBildn": {
        "supported": False,
        "description": "Amber ff99SB-ILDN; requires parent/custom-template compatibility validation",
    },
}


def validate_force_field(name: str) -> str:
    """Validate a selection; numerical compatibility is reported separately."""
    if name in FORCE_FIELDS:
        return name
    raise DependencyError(
        f"Force field {name!r} is not supported for this custom RESP library. "
        "Choose ff14SB, ff19SB or ff99SBildn; each still requires compatible, validated custom templates."
    )


def _python_probe(python: str | Path) -> dict:
    # Import in the interpreter being checked, not in the bootstrap interpreter.
    script = r'''
import contextlib, importlib, importlib.metadata as md, io, json, sys
required = json.loads(sys.argv[1])
out = {'executable':sys.executable, 'version':list(sys.version_info[:3]),
       'isolated':sys.prefix != sys.base_prefix, 'prefix':sys.prefix, 'packages':{}}
for distribution, spec in required.items():
    module, minimum = spec
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            importlib.import_module(module)
        out['packages'][distribution] = {'importable':True, 'version':md.version(distribution), 'minimum':minimum}
    except Exception as exc:
        out['packages'][distribution] = {'importable':False, 'version':None, 'minimum':minimum, 'error':type(exc).__name__}
print(json.dumps(out))
'''
    try:
        result = subprocess.run(
            [str(python), "-I", "-B", "-c", script, json.dumps(PYTHON_PACKAGES)],
            capture_output=True, text=True, timeout=90,
        )
        if result.returncode:
            raise DependencyError("The selected Python could not complete dependency checks.")
        data = json.loads(result.stdout)
    except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
        raise DependencyError("The selected Python is unavailable or failed dependency checks.") from exc
    data["usable"] = (3, 10) <= tuple(data["version"]) < (3, 15)
    for entry in data["packages"].values():
        found = tuple(int(x) for x in re.findall(r"\d+", entry["version"] or "")[:3])
        minimum = tuple(int(x) for x in entry["minimum"].split("."))
        entry["satisfies"] = entry["importable"] and found >= minimum
    return data


def find_executable(name: str, root: Path, environ: dict | None = None) -> str | None:
    """Prefer an existing Amber installation or a project-private APBS."""
    environ = os.environ if environ is None else environ
    candidates = []
    if name in AMBER_PROGRAMS and environ.get("AMBERHOME"):
        candidates.append(Path(environ["AMBERHOME"]).expanduser() / "bin" / name)
    if name == "apbs":
        candidates.append(Path(root) / ".tools" / "apbs" / "bin" / name)
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.absolute())
    found = shutil.which(name, path=environ.get("PATH", ""))
    return str(Path(found).absolute()) if found else None


def external_program_environment(executable: str | Path, environ: dict | None = None) -> dict:
    """Find APBS companion libraries in the child process, without shell changes."""
    env = dict(os.environ if environ is None else environ)
    executable = Path(executable).absolute()
    library = executable.parent.parent / "lib"
    if executable.name == "apbs" and library.is_dir():
        for name in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
            env[name] = str(library) + (os.pathsep + env[name] if env.get(name) else "")
        env["PATH"] = str(executable.parent) + (os.pathsep + env["PATH"] if env.get("PATH") else "")
    return env


def _force_field_report(root: Path, python: str | Path) -> dict:
    """Ask the force-field registry using the target environment's OpenMM."""
    script = """
import json,sys
sys.path.insert(0, sys.argv[1])
from cationize_protein.forcefields import force_field_status
print(json.dumps({name:force_field_status(name) for name in json.loads(sys.argv[2])}, default=str))
"""
    try:
        result = subprocess.run([str(python), "-I", "-B", "-c", script, str(root / "src"), json.dumps(list(FORCE_FIELDS))],
                                capture_output=True, text=True, timeout=60)
        if result.returncode:
            raise ValueError("registry probe failed")
        statuses = json.loads(result.stdout)
        return {name: {**FORCE_FIELDS[name], **status, "supported": bool(status.get("available"))}
                for name, status in statuses.items() if name in FORCE_FIELDS}
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return {name: {**entry, "supported": False, "error": "Force-field registry could not be checked in the selected Python"}
                for name, entry in FORCE_FIELDS.items()}


def inspect_dependencies(root: Path, python: str | Path | None = None,
                         environ: dict | None = None, engine: str = "openmm",
                         solver: str = "coulomb", force_field: str = "ff14SB") -> dict:
    """Return a JSON-ready read-only report for one explicit engine/model choice.

    Executable detection is not a numerical validation or a licence determination.
    This function imports third-party modules only in a subprocess of ``python``.
    """
    root = Path(root).absolute()
    programs = {name: find_executable(name, root, environ) for name in (*AMBER_PROGRAMS, *GROMACS_PROGRAMS, "apbs")}
    try:
        py = _python_probe(python or sys.executable)
    except DependencyError as exc:
        py = {"executable": str(python or sys.executable), "usable": False,
              "packages": {}, "error": str(exc)}
    runtime = bool(py["usable"] and all(p.get("satisfies") for p in py["packages"].values())
                   and len(py["packages"]) == len(PYTHON_PACKAGES))
    amber = any(programs[name] for name in ("pmemd", "sander"))
    capabilities = {
        "runtime": {"available": runtime},
        "openmm_minimization": {"available": bool(py["packages"].get("openmm", {}).get("satisfies"))},
        "amber_minimization": {"available": amber, "note": "Executable detected; simulation and licence not tested"},
        "amber_topology": {"available": bool(programs["tleap"])},
        "gromacs": {"available": any(programs[name] for name in GROMACS_PROGRAMS),
                    "note": "Optional external engine; not needed to export GROMACS files"},
        "apbs": {"available": bool(programs["apbs"]), "model": "Poisson-Boltzmann; executable presence only"},
        "pbsa": {"available": bool(programs["pbsa"]), "model": "Amber PBSA; executable presence only"},
        "coulomb": {"available": True, "model": "Screened Coulomb approximation; not a Poisson-Boltzmann solver"},
    }
    engine_available = (capabilities["openmm_minimization"]["available"] if engine == "openmm" else
                        amber if engine == "amber" else bool(programs.get(engine)))
    force_fields = _force_field_report(root, python or sys.executable)
    supported_ff = bool(force_fields.get(force_field, {}).get("supported"))
    from .optional_modelers import inspect_modelers
    modelers = inspect_modelers(python or sys.executable, environ)
    selected = {"engine": engine, "solver": solver, "force_field": force_field,
                "engine_available": engine_available,
                "solver_available": bool(capabilities.get(solver, {}).get("available")),
                "force_field_supported": supported_ff}
    return {"schema_version": 1, "project_root": str(root), "python": py,
            "executables": programs, "capabilities": capabilities, "force_fields": force_fields,
            "optional_modelers": modelers,
            "selected": selected,
            "ready": runtime and engine_available and selected["solver_available"] and supported_ff,
            "scope": "Dependencies only; chemistry, custom templates, solver output and exports require separate validation"}


def format_report(report: dict) -> str:
    lines = [f"Project: {report['project_root']}", f"Python: {report['python']['executable']}"]
    if report["python"].get("error"):
        lines.append(report["python"]["error"])
    for name, entry in report["python"]["packages"].items():
        status = entry["version"] if entry.get("satisfies") else "missing, broken or below minimum"
        lines.append(f"  {name}: {status} (requires >= {entry['minimum']})")
    lines.append("External programs (detected only; computational validation is separate):")
    lines.extend(f"  {name}: {path or 'not found'}" for name, path in report["executables"].items())
    lines.append("GROMACS: " + ("detected" if report['capabilities']['gromacs']['available'] else "not found (optional)"))
    lines.append('Optional structure modeling (not required for PDBFixer preparation):')
    for entry in report.get('optional_modelers', {}).values():
        lines.append(f'  {entry["name"]}: {entry["reason"]}')
        lines.append(f'    Python: {entry["python"]}; installation: {entry["url"]}')
    lines.append("Python runtime: " + ("complete" if report['capabilities']['runtime']['available'] else "incomplete"))
    choice = report["selected"]
    lines.append(f"Selected: engine={choice['engine']}, potential={choice['solver']}, force field={choice['force_field']}")
    if choice["solver"] == "coulomb":
        lines.append("Screened Coulomb is an explicit approximation, not a Poisson-Boltzmann solver.")
    if not choice["engine_available"]:
        lines.append("Selected minimization engine is unavailable; no engine was substituted.")
    if not choice["solver_available"]:
        lines.append("Selected potential solver is unavailable; no physical model was substituted.")
    if not choice["force_field_supported"]:
        lines.append("The selected protein force-field definition is unavailable or could not be checked.")
    scope = report["force_fields"].get(choice["force_field"], {}).get("validation_scope")
    if scope:
        lines.append(f"Force-field validation scope: {scope}")
    lines.append("Selected workflow dependencies are available." if report["ready"] else
                 "Some selected workflow dependencies are unavailable; install them before requesting that workflow.")
    lines.append("Amber, GROMACS, APBS, MODELLER and ESMFold2 can be installed later; discovery is repeated when requested.")
    lines.append(report["scope"])
    return "\n".join(lines)


def _environment_python(target: Path) -> Path:
    return target / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def _validate_target(root: Path, target: Path) -> Path:
    root = root.resolve()
    target = Path(os.path.abspath(target.expanduser()))
    # macOS aliases /var to /private/var. Reject symlinks *inside* the
    # checkout before canonicalizing, while allowing its external path aliases.
    for part in (target, *target.parents):
        if part.resolve() == root:
            break
        if part.is_symlink():
            raise DependencyError("Environment path must not contain symlinks.")
    target = target.resolve()
    if target == root or target in root.parents:
        raise DependencyError("Environment path must not be the project directory or its parent.")
    # Limit private environments to this checkout. This also prevents a typo from
    # turning an unrelated environment elsewhere into an install destination.
    try:
        relative = target.relative_to(root)
    except ValueError as exc:
        raise DependencyError("Environment path must be inside this project checkout.") from exc
    if relative.parts[0] in {"src", "tests", "parameters", "config", "docs", "tools", ".tools", ".git"}:
        raise DependencyError("Environment path overlaps project source or data.")
    return target


def _create_venv(python: str, target: Path) -> None:
    result = subprocess.run([python, "-m", "venv", str(target)], capture_output=True, text=True)
    if result.returncode:
        raise DependencyError("Could not create the private environment. Install Python 3.12 with venv/ensurepip support and rerun; no system Python was changed.")


def ensure_environment(root: Path, target: Path, python: str) -> Path:
    """Create or reuse only an installer-owned isolated environment."""
    target = _validate_target(root, target)
    marker = target / ENV_MARKER
    if target.exists():
        if not target.is_dir() or not marker.is_file() or marker.is_symlink():
            raise DependencyError("Existing environment directory is not owned by this installer; choose a fresh path.")
        try:
            content = json.loads(marker.read_text())
        except (OSError, ValueError) as exc:
            raise DependencyError("Environment ownership marker is invalid.") from exc
        if content != {"owner": "cationize-protein", "schema_version": 1}:
            raise DependencyError("Environment ownership marker is invalid.")
        probe = _python_probe(_environment_python(target))
        if not probe["usable"] or not probe.get("isolated") or Path(probe.get("prefix", "")).resolve() != target.resolve():
            raise DependencyError("Existing environment does not match its isolated Python; it was left untouched.")
        return _environment_python(target)
    try:
        target.mkdir(parents=True, exist_ok=False)
        marker.write_text(json.dumps({"owner": "cationize-protein", "schema_version": 1}) + "\n")
        _create_venv(python, target)
        probe = _python_probe(_environment_python(target))
        if not probe["usable"] or not probe.get("isolated"):
            raise DependencyError("New environment is not an isolated supported Python (3.10-3.14).")
    except (OSError, subprocess.SubprocessError) as exc:
        raise DependencyError("Private environment creation failed; the partial directory was retained for inspection.") from exc
    return _environment_python(target)


def install_python_package(root: Path, python: Path) -> None:
    """Use the private interpreter's pip, keeping pip output out of credential logs."""
    env = dict(os.environ)
    # Explicit per-user pip options must not redirect writes outside the venv.
    for name in list(env):
        if name.startswith("PIP_") or name in {"PYTHONHOME", "PYTHONPATH"}:
            env.pop(name)
    env["PIP_CONFIG_FILE"] = os.devnull
    env["PIP_REQUIRE_VIRTUALENV"] = "true"
    env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
    try:
        result = subprocess.run(
            [str(python), "-m", "pip", "install", "--no-input", "--no-cache-dir", "--disable-pip-version-check", str(root.resolve())],
            env=env, capture_output=True, text=True,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise DependencyError("Private pip could not run. No global install was attempted.") from exc
    if result.returncode:
        # pip may echo authenticated index URLs, so never include its raw output.
        raise DependencyError("Private pip installation failed. Check network access, Python-version wheel support and the project dependencies; no global install was attempted.")


def extract_apbs_archive(archive: Path, destination: Path, bundle: str) -> None:
    """Copy a release bundle without extractall, links or unsafe archive paths."""
    if destination.exists() or destination.is_symlink():
        raise DependencyError("APBS extraction destination must be new.")
    with zipfile.ZipFile(archive) as source:
        entries = []
        seen = set()
        expanded = 0
        for item in source.infolist():
            name = item.filename
            parts = PurePosixPath(name).parts
            mode = item.external_attr >> 16
            raw_parts = name.rstrip("/").split("/")
            if (not parts or name.startswith("/") or "\\" in name or ":" in name or
                    any(p in {"..", ".", ""} for p in raw_parts) or
                    any(p in {"..", ".", ""} for p in parts) or parts[0] != bundle or
                    stat.S_ISLNK(mode) or (stat.S_IFMT(mode) not in (0, stat.S_IFREG, stat.S_IFDIR)) or
                    name in seen):
                raise DependencyError("APBS archive contains an unsafe path, link or file type.")
            seen.add(name)
            expanded += item.file_size
            if expanded > 1024 * 1024 * 1024 or len(seen) > 100000:
                raise DependencyError("APBS archive exceeds safe extraction limits.")
            if len(parts) == 1 and not item.is_dir():
                raise DependencyError("APBS archive has an invalid bundle root.")
            entries.append((item, parts[1:]))
        if not any(parts == ("bin", "apbs") and not item.is_dir() for item, parts in entries):
            raise DependencyError("APBS archive does not contain bin/apbs.")
        # Validation is complete before the first write.
        destination.mkdir(parents=True)
        for item, parts in entries:
            path = destination.joinpath(*parts)
            if item.is_dir():
                path.mkdir(parents=True, exist_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                with source.open(item) as incoming, path.open("xb") as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
                path.chmod(0o755 if item.external_attr >> 16 & 0o111 else 0o644)


def install_apbs(root: Path) -> str:
    """Install one SHA256-pinned official APBS release into .tools/apbs."""
    machine = platform.machine().lower()
    machine = "x86_64" if machine in {"amd64", "x64"} else machine
    key = f"{platform.system()}-{machine}"
    try:
        manifest = json.loads((root / "config" / "install-manifest.json").read_text())
        asset = manifest["apbs"]["assets"].get(key)
    except (OSError, ValueError, KeyError) as exc:
        raise DependencyError("APBS install manifest is missing or invalid.") from exc
    if not asset:
        raise DependencyError(f"No pinned official APBS binary for {key}. Supply an existing compatible apbs on PATH or build APBS locally using the official instructions; no solver was substituted.")
    url = asset["url"]
    parsed = urlparse(url)
    if (parsed.scheme != "https" or parsed.netloc != "github.com" or
            not parsed.path.startswith("/Electrostatics/apbs/releases/download/") or
            not re.fullmatch(r"[a-f0-9]{64}", asset["sha256"])):
        raise DependencyError("APBS asset must use an official HTTPS release URL and a pinned SHA256.")
    tools = root / ".tools"
    destination = tools / "apbs"
    if tools.is_symlink() or (tools.exists() and not tools.is_dir()) or destination.exists() or destination.is_symlink():
        raise DependencyError("Private APBS destination already exists or is a symlink; it was left untouched.")
    tools.mkdir(exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(prefix=".apbs-install-", dir=tools) as temporary:
            temporary = Path(temporary)
            archive = temporary / "release.zip"
            digest = hashlib.sha256()
            size = 0
            with urlopen(Request(url, headers={"User-Agent": "cationize-protein-installer"}), timeout=60) as incoming, archive.open("xb") as outgoing:
                # GitHub redirects release objects to its own release CDN.
                final = urlparse(incoming.geturl())
                if final.scheme != "https" or final.hostname not in {"github.com", "release-assets.githubusercontent.com", "objects.githubusercontent.com"}:
                    raise DependencyError("APBS download redirected outside the official release hosts.")
                while chunk := incoming.read(1024 * 1024):
                    size += len(chunk)
                    if size > 200 * 1024 * 1024:
                        raise DependencyError("APBS download exceeds the allowed archive size.")
                    digest.update(chunk)
                    outgoing.write(chunk)
            if digest.hexdigest() != asset["sha256"]:
                raise DependencyError("APBS download SHA256 mismatch; no executable was installed.")
            staged = temporary / "unpacked"
            extract_apbs_archive(archive, staged, asset["bundle"])
            staged.rename(destination)
    except DependencyError:
        raise
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise DependencyError("Official APBS download or private extraction failed; no global changes were made.") from exc
    return str(destination / "bin" / "apbs")


def _select_python(explicit: str | None) -> str:
    choices = ([explicit] if explicit else
               [shutil.which(name) for name in ("python3.12", "python3.13", "python3.11", "python3.14", "python3.10")] + [sys.executable])
    seen = set()
    for candidate in choices:
        if not candidate:
            continue
        resolved = shutil.which(candidate) if not Path(candidate).is_absolute() else candidate
        if not resolved or resolved in seen:
            continue
        seen.add(resolved)
        try:
            version = subprocess.run([str(resolved), "-I", "-c", "import sys,json; print(json.dumps(list(sys.version_info[:3])))"],
                                     capture_output=True, text=True, timeout=15)
            if version.returncode == 0 and (3, 10) <= tuple(json.loads(version.stdout)) < (3, 15):
                return str(Path(resolved).absolute())
        except (OSError, ValueError, subprocess.TimeoutExpired):
            continue
    raise DependencyError("A working Python 3.10-3.14 with venv is required for the current OpenMM wheels; Python 3.12 is recommended. Install a compatible Python separately or select it with --python. The installer does not replace system Python.")


def installer_main(argv: list[str] | None = None, *, root: Path) -> int:
    parser = argparse.ArgumentParser(description="Install the complete Python runtime in a private venv; report optional Amber, GROMACS, APBS, MODELLER and ESMFold2 without downloading them.")
    parser.add_argument("--check", action="store_true", help="Read-only checks; do not create or install anything")
    parser.add_argument("--venv", type=Path, help="Dedicated environment inside this checkout (default: .venv)")
    parser.add_argument("--python", help="Python used to create the venv (prefers python3.12)")
    parser.add_argument("--engine", choices=("openmm", "amber", "pmemd", "sander"), help="Optional saved minimizer preference; choose later at runtime (default: openmm)")
    parser.add_argument("--solver", choices=("apbs", "pbsa", "coulomb"), help="Saved potential preference (default: coulomb, included screened approximation); APBS/PBSA are separate installations")
    parser.add_argument("--force-field", choices=tuple(FORCE_FIELDS), help="Optional saved protein force-field preference; choose later at runtime (default: ff14SB)")
    parser.add_argument("--json", action="store_true", help="Emit a machine-readable dependency report")
    args = parser.parse_args(argv)
    root = Path(root).resolve()
    target = args.venv or root / ".venv"
    if not target.is_absolute():
        target = root / target
    try:
        target = _validate_target(root, target)
        settings = target / "cationize-protein-settings.json"
        if settings.is_symlink() or (settings.exists() and not settings.is_file()):
            raise DependencyError("Saved installation settings must be a regular file, not a symlink.")
        try:
            stored = json.loads(settings.read_text()) if settings.is_file() else {}
        except (OSError, ValueError) as exc:
            raise DependencyError("Saved installation settings are invalid; inspect cationize-protein-settings.json.") from exc
        if not isinstance(stored, dict) or stored.get("schema_version", 1) != 1:
            raise DependencyError("Saved installation settings have an unsupported format.")
        force_field = validate_force_field(args.force_field or stored.get("force_field", "ff14SB"))
        engine = args.engine or stored.get("engine", "openmm")
        solver = args.solver or stored.get("solver", "coulomb")
        if engine not in {"openmm", "amber", "pmemd", "sander"} or solver not in {"apbs", "pbsa", "coulomb"}:
            raise DependencyError("Saved installation settings contain an unknown engine or solver.")
        if args.check:
            python = _environment_python(target) if target.exists() else _select_python(args.python)
            report = inspect_dependencies(root, python, engine=engine, solver=solver, force_field=force_field)
            report["environment"] = {"target": str(target), "exists": target.exists()}
            if not target.exists():
                report["ready"] = False
                report["environment"]["note"] = "Private environment is absent; Python packages shown belong to the inspected bootstrap interpreter. Run install.py to create the environment."
        else:
            if not args.json:
                print(f"Checking the private Python environment: {target}", flush=True)
            python = ensure_environment(root, target, _select_python(args.python))
            if not args.json:
                print("Installing missing Python requirements in the private environment...", flush=True)
            install_python_package(root, python)
            report = inspect_dependencies(root, python, engine=engine, solver=solver, force_field=force_field)
            report["environment"] = {"target": str(target), "exists": True}
            report['installation_complete'] = bool(report['capabilities']['runtime']['available'])
            (target / "cationize-protein-settings.json").write_text(json.dumps({
                "schema_version": 1, "engine": engine, "solver": solver, "force_field": force_field,
                "project_root": str(root), "executables": report["executables"],
            }, indent=2) + "\n")
        if args.json:
            print(json.dumps(report, indent=2))
        else:
            print(format_report(report))
            if report.get("installation_issue"):
                print(report["installation_issue"])
            if report["environment"].get("note"):
                print(report["environment"]["note"])
            if not args.check:
                print("Python package installation complete." if report['installation_complete'] else
                      "Python package installation incomplete; inspect the dependency report.")
                if os.name=='nt':
                    print(f"Activate in PowerShell: & '{python.parent / 'Activate.ps1'}'")
                else:
                    activation=os.path.relpath(python.parent/'activate',Path.cwd())
                    print('Activate: source '+shlex.quote(activation))
                print('Then run: cationize-protein --help')
                if solver == "apbs" and report["executables"]["apbs"]:
                    print(f"APBS executable: {report['executables']['apbs']}")
        complete=report['ready'] if args.check else report['installation_complete']
        return 0 if complete else 1
    except DependencyError as exc:
        if args.json:
            print(json.dumps({"ready": False, "error": str(exc)}))
        else:
            print(f"Installation check: {exc}", file=sys.stderr)
        return 1
