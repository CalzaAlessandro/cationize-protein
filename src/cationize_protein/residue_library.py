"""Copy hash-verified residue OFF libraries into a portable run output."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile


def export_residue_library(source, destination):
    source,destination=Path(source),Path(destination)
    manifest=json.loads((source/'manifest.json').read_text())
    allowed={'validated','reused_validated_HSA_template','RESP_capped_Amber_and_OpenMM_validated'}
    accepted=[]
    for row in manifest.get('templates',[]):
        if row.get('status') not in allowed:continue
        name=row.get('name','')
        if len(name)!=3 or not name.isalnum():raise ValueError('Invalid residue template name')
        path=source/(name+'.lib')
        if not path.is_file() or path.is_symlink():raise ValueError(f'Validated template {name} is missing')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row.get('library_sha256'):
            raise ValueError(f'Validated template {name} hash differs from its manifest')
        accepted.append(name)
    data={path.name:path.read_bytes() for path in source.iterdir() if path.is_file() and
          (path.suffix in {'.json','.xml','.csv'} or path.name in {name+'.lib' for name in accepted})}
    text='''# Individual cationized residue templates

Each NAME.lib is one residue in Amber OFF format. The .lib extension is a
filename convention; a duplicate .off file is not required. Original template
bytes and fitted charges are preserved. manifest.json records chemistry,
terminal context, validation status and source hashes.

Example: ASC.lib is the internal Asp product; ASX.lib is its N-terminal
counterpart. Use the template matching the actual residue and chain context.

These templates use the source Amber type names and need compatible Amber
force-field parameters. They are distinct from ../25/exports/protein.lib,
which contains the complete protein unit and uses its own printed type map.
Do not combine a raw residue library with unrelated remapped atom types.

The ready-to-use protein files and full numerical parameters are in each
percentage's exports directory. For ff19SB, retain the CMAP information in
the numerical topology or XML/JSON bundle; frcmod alone does not encode CMAP.

The bundled cationized_residues.xml is the ff14SB residue extension. Other
profiles are assigned by cationize-protein from their native parent contexts.
'''
    data['README.md']=text.encode()
    record=dict(schema_version=1,template_count=len(accepted),templates=sorted(accepted),
                files_sha256={name:hashlib.sha256(value).hexdigest() for name,value in sorted(data.items())})
    if destination.exists() or destination.is_symlink():
        if destination.is_symlink() or not destination.is_dir():raise ValueError('Unsafe existing residue library')
        marker=destination/'LIBRARY_EXPORT.json'
        if marker.is_symlink() or not marker.is_file() or json.loads(marker.read_text())!=record:
            raise ValueError('Existing residue library manifest differs; no files were overwritten')
        for name,digest in record['files_sha256'].items():
            path=destination/name
            if path.is_symlink() or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
                raise ValueError(f'Existing residue library file {name} changed or has a hash mismatch')
        return record
    with tempfile.TemporaryDirectory(prefix='.residue-library-',dir=destination.parent) as tmp:
        staged=Path(tmp)/'library';staged.mkdir()
        for name,value in data.items():(staged/name).write_bytes(value)
        (staged/'LIBRARY_EXPORT.json').write_text(json.dumps(record,indent=2)+'\n')
        os.rename(staged,destination)
    return record
