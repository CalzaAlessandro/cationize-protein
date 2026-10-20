"""Run development fits one after another, without a batch scheduler."""
import logging
from pathlib import Path
from .catalog import templates

LOG = logging.getLogger(__name__)


def run_serial(project, names=None):
    from .parameter_jobs import worker
    from .parameter_batch import collect
    project = Path(project).resolve()
    available = {s['name'] for s in templates() if not s['reused_from']}
    names = sorted(available) if names is None else list(names)
    if len(names) != len(set(names)) or not set(names).issubset(available) or not names:
        raise ValueError('Select unique new-fit template names; the four reused HSA fits are preserved')
    done = []
    try:
        for index, name in enumerate(names, 1):
            LOG.info('Parameterization %d/%d: %s', index, len(names), name)
            worker(project, name)
            done.append(name)
    finally:
        collect(project)
    return dict(completed=done, execution='sequential; no scheduler')
