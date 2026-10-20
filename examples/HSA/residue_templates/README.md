# Individual cationized residue templates

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
