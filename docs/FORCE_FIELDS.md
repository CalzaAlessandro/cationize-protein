# Force fields, physical exports and complete parameters

The protein field is selected independently of the minimizer and potential solver.

| Profile | OpenMM native source | Product parameters | CMAP |
|---|---|---|---|
| ff14SB | amber14/protein.ff14SB.xml | Original RESP charges and Amber protein-type bonded/LJ analogy | No |
| ff19SB | amber19/protein.ff19SB.xml | Same fitted RESP; selected native parent-context types | Parent-specific backbone map retained |
| ff99SBildn | amber99sbildn.xml | Same fitted RESP; explicit native-type/class transfer | No |

OpenMM 8.6+ supplies all three profiles. ff19SB custom backbone types inherit
the native parent/context type, including the zeroed ordinary backbone torsions
and corresponding CMAP selection. Its CMAP grid/index survives Amber/GROMACS
conversion and is tested numerically. The original RESP fit is not repeated
for each profile. See the [ff19SB authors' implementation notes](https://github.com/csimmerling/ff19SB_201907/blob/master/doc/forcefield_files)
and [paper](https://pubmed.ncbi.nlm.nih.gov/31714766/).

CMT, the terminal product of a native CYM thiolate, currently lacks the required
CCYM parent in ff19SB/ff99SBildn. Selecting that actual product fails explicitly.
An unrelated unavailable product does not prevent use of the other 31 templates.
The default ff14SB library supports all 32. Native protonation states are checked
against actual installed context templates; e.g. the ff14SB distribution has
internal ASH/GLH/LYN/CYM but does not provide every corresponding terminal state.

## Complete parameter bundle

Every checkpoint export contains the full physical `system.xml`,
`PARAMETERS.json` and atom-type/charge CSV, plus all ordinary parameters in
`parameters.frcmod` and a complete `protein.lib`. These are self-contained for
the system that was generated. They do not load a different Amber data tree
implicitly. Amber OFF libraries are used; PREP/PREPC is not a runtime requirement.

`PARAMETERS.json` prints atomic charges, masses, LJ values, every bond/angle,
proper/improper torsion, 1-4/exclusion entry, CMAP indices and complete grids.
`rebuild.py` reconstructs Amber files from the physical XML and full-precision
JSON with OpenMM+ParmEd. `rebuild.in` permits direct LEaP reconstruction for
systems without CMAP. With CMAP it stops and points to the Python rebuild;
a conventional frcmod alone must not silently lose the map.

Amber parm7 contains numeric parameters and does not require a matching leaprc.
Executable support still matters: ff19SB requires an engine that implements its
CMAP terms. The authors originally required patched Amber 18/AmberTools 19;
Amber 20+ is recommended for these maps. No file export can provide missing
functionality in an older executable. Existing fits and HSA production settings
are preserved separately.

## Validation of formats

Physical exports exclude GB and minimization restraints. Checks compare the
atom graph/order, q, LJ sigma/epsilon, 1-4/exclusion graph and coefficients, CMAP
count, vacuum energy and atomic forces after reloading with ParmEd/OpenMM at
identical source coordinates. A zero-epsilon particle's sigma is physically
irrelevant and may change convention; epsilon, charge and actual energy/forces
must still agree. Coordinate file precision differs by format.

Amber exports parm7/rst7, OpenMM exports physical XML/JSON, and GROMACS exports
top/gro with separate forcefield.itp/protein.itp. A PDB is always included.
GROMACS atom types use a CP_ namespace and all includes are relative. A
successful ParmEd/OpenMM round trip is recorded separately from an actual
external-engine single-point test. No MD is launched by these conversions.

OpenMM minimization uses GBn2 and one CPU thread by default. Amber minimization
uses installed serial sander/pmemd, GBn2/mbondi3 and explicit atom-index restraint
groups. GPU and MPI executables are not used by this serial backend. Both minimizers use unscreened GBn2 (0M salt). Subsequent
PB ranking uses its separately recorded intrinsic radii/dielectric/salt settings.
The two minimizers are separate implementations and need not reach identical
minima. The selected engine/field/model is recorded in immutable run provenance.
