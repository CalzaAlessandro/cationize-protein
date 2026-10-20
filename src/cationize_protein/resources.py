"""Find the custom residue library in a checkout or an installed wheel."""
from pathlib import Path
import sysconfig

def default_parameters():
    source = Path(__file__).resolve().parents[2]/'parameters'
    if (source/'manifest.json').is_file():
        return source
    return Path(sysconfig.get_path('data'))/'share/cationize-protein/parameters'
