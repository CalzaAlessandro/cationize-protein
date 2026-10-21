"""Read-only discovery of optional modeling software in separate Pythons."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

SOFTWARE = {
    'modeller': {'name': 'MODELLER', 'environment': 'CATIONIZE_MODELLER_PYTHON',
                 'url': 'https://salilab.org/modeller/'},
    'esmfold': {'name': 'ESMFold2', 'environment': 'CATIONIZE_ESMFOLD_PYTHON',
                'url': 'https://huggingface.co/biohub/ESMFold2'},
}

_PROBE = r'''
import importlib.metadata as md, importlib.util as util, json, pathlib, sys
method = sys.argv[1]
def package(name):
    try:
        spec = util.find_spec(name)
        return spec
    except (ImportError, ValueError):
        return None
def version(name):
    try: return md.version(name)
    except md.PackageNotFoundError: return None
if method == 'modeller':
    spec = package('modeller')
    available = bool(spec and any((pathlib.Path(p)/'automodel').is_dir()
                    for p in (spec.submodule_search_locations or [])))
    result = dict(available=available, version=version('modeller'),
        reason='code detected; licence is checked only when modeling starts' if available else
               'not installed in this Python; install the official MODELLER distribution and configure its licence')
else:
    spec = package('transformers')
    folding = bool(spec and any((pathlib.Path(p)/'models/esmfold2/modeling_esmfold2.py').is_file()
                   for p in (spec.submodule_search_locations or [])))
    available = bool(package('torch') and folding and package('accelerate'))
    result = dict(available=available, version=version('transformers'),
        reason='Python code detected; model weights and hardware are checked when modeling starts' if available else
               'not installed in this Python; PyTorch, Accelerate and Transformers with EsmFold2Model are needed')
print(json.dumps(result))
'''


def inspect_modelers(python=None, environ=None):
    """Find code without importing modelers, testing licences or loading weights.

    Explicit external interpreters never fall back to the package environment.
    Otherwise, other Python interpreters on PATH are checked after the active one.
    The executable path is not resolved through symlinks, to retain venv identity.
    """
    env = os.environ if environ is None else environ
    report = {}
    for method, spec in SOFTWARE.items():
        configured = env.get(spec['environment'])
        candidates = [configured or str(python or sys.executable)]
        if not configured:
            candidates += [found for name in ('python3.15','python3.14','python3.13','python3.12',
                                              'python3.11','python3.10','python3','python')
                           if (found := shutil.which(name, path=env.get('PATH', '')))]
        checked = set()
        for executable in candidates:
            found = shutil.which(os.path.expanduser(executable), path=env.get('PATH', ''))
            if found:
                found = os.path.abspath(found)
            if (found or executable) in checked:
                continue
            checked.add(found or executable)
            entry = dict(name=spec['name'], url=spec['url'], python=found or executable,
                         available=False, version=None, optional=True,
                         license_validated=None, weights_available=None)
            if not found:
                entry['reason'] = f'Optional {spec["name"]} Python interpreter was not found: {executable}'
            else:
                try:
                    result = subprocess.run([found, '-I', '-B', '-c', _PROBE, method],
                        capture_output=True, text=True, timeout=20, env=dict(env))
                    data = json.loads(result.stdout)
                    if result.returncode or not isinstance(data, dict) or not isinstance(data.get('available'), bool):
                        raise ValueError
                    entry.update(available=data['available'], version=data.get('version'),
                                 reason=data.get('reason', 'Python code detected; calculation not tested'))
                except (OSError, ValueError, subprocess.TimeoutExpired):
                    entry['reason'] = f'Optional {spec["name"]} interpreter could not complete its read-only probe'
            if method not in report or entry['available']:
                report[method] = entry
            if entry['available']:
                break
        report[method]['interpreters_checked'] = len(checked)
    return report
