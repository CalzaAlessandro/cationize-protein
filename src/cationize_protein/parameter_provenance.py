"""Verify the actual source version used for an immutable historical fit."""
from pathlib import Path
import re
from .qm import digest

def known_resp_source(project, sha):
    if not isinstance(sha, str) or not re.fullmatch('[0-9a-f]{64}', sha):
        return False
    candidates = [Path(__file__).with_name('resp.py'),
                  Path(project)/'parameters/provenance'/f'resp_{sha}.py']
    return any(p.is_file() and not p.is_symlink() and digest(p) == sha for p in candidates)
