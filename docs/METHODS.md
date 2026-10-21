# Parameterization and electrostatic ranking

## Chemistry and charge states

The target reaction replaces one carboxylate oxygen by an ethylenediamine
amide, yielding `–C(=O)–NH–CH2–CH2–NH3+`. The introduced ammonium is held
protonated. Carboxylates on Asp/Glu side chains and biological C termini
are distinct candidate groups. N-terminal, C-terminal, histidine-tautomer,
thiolate and disulfide contexts have separate templates where needed.

The 32-state catalog covers the 20 standard amino acids at the C terminus,
with HIE as the default histidine tautomer, plus HID/HIP and CYX/CYM
variants; internal and N-terminal Asp/Glu side-chain products; and
C-terminal acidic residues with side-only or both carboxyls modified.
Unmodified ASH/GLH states retain their acid hydrogen and remain chemical
carboxyl candidates. State transitions on already amidated products need an
explicit matching fitted template. Unsupported post-translational modifications
are rejected.

## Capped quantum models and RESP

Native geometry, atom types and cap charges are read from the installed
Amber definitions. Internal residues use ACE/NME caps, N-terminal residues
use NME, and C-terminal residues use ACE. For the disulfide variant an
explicit neutral ACE–CYX–NME partner preserves the S–S bond. Caps and this
partner retain their standard charges during fitting. Existing four HSA
fits are preserved exactly under the renamed identifiers.

For each new template, two initial backbone conformations are optimized
with Gaussian at HF/6-31G(d), using a frozen backbone torsion and frozen
heteroatom–H bond lengths to preserve the selected protonation state.
These constraints are part of this implementation's model; the resulting
geometries are not fully unconstrained gas-phase minima. Merz–Kollman ESP
samples feed two-stage, two-conformer RESP. Equivalent atom charges are
shared across conformers. Element identities in RESP inputs are checked
against explicit atomic numbers, independently of the exported atom names.

The reacted carbonyl uses the ff14SB C/O amide types. Other inherited
types are preserved, with protein atom types assigned to the introduced
amide, methylenes and ammonium. Bond, angle, torsion and Lennard-Jones
terms are transferred by chemical analogy from the protein force field.
No GAFF parameter assignment or new QM torsion fitting is performed.

Acceptance requires normal Gaussian termination, completed optimization
and ESP generation; intact bonds/protonation and inherited stereochemistry;
fixed cap charges; equivalent conformer charge sets; the expected residue
charge within 3 × 10⁻⁵ e; finite ESP relative RMS below the documented 0.30
screen; complete ff14SB capped-model parameters in tleap; and finite
OpenMM energies/forces with exact fitted charge retention in both ff14SB
and GBn2. The ESP threshold screens failed fits; parameter quality also
requires independent validation. Model, log, library and source hashes
connect all accepted outputs to the exact calculation.

These checks establish internal consistency and capped-model
compatibility. Independent peptides, multiple proteins, torsional
energetics and experimental comparisons remain necessary before a claim
of broad transferability.

## Progressive site ranking and steric relaxation

The default screened-Coulomb solver evaluates the current atomic charges in
a uniform solvent dielectric with Debye–Hückel screening. Its defaults are
298 K, dielectric 78.5 and 0.15 M monovalent salt. This approximation omits
the protein dielectric boundary and reaction field and runs entirely in Python.
APBS and Amber PBSA provide optional Poisson–Boltzmann calculations; APBS uses
the nonlinear equation with solute/solvent dielectric constants 2/78.5.
Solver identity and settings are recorded for every run, and an existing run
retains its recorded solver on resumption. The models can rank sites differently.
Candidate sampling
uses a 1 Å lattice in the union of 5 Å spheres around its two carboxylate
oxygens. Points overlapping atomic LJ exclusion radii plus the explicitly
selected positive-probe radius are removed. The most negative accessible
local potential ranks first; deterministic identifiers break exact ties.

Substitution preserves the remaining protein graph, peptide bonds,
disulfides and inherited coordinates. Deterministic bond rotations address
steric overlap while protecting backbone atoms, rings and disulfides.
Restrained OpenMM GBn2 minimization then relaxes the structure with ff14SB
and the custom RESP templates by default. Installed serial Amber and
ff19SB/ff99SBildn profiles are explicit alternatives. Each field is recomputed
from the updated state. This ranking is an electrostatic model for choosing sites, not a
prediction of chemical reaction rates.

## Native protonation and full parameters

PROPKA audits the baseline and each substitution at pH7.4 using generic-amine
and fixed+1 controls. Native changes require consensus at least one pH unit
from the target and persistence over two chemical steps. Compatible changes
rebuild side-chain hydrogens, reassign charges/parameters, minimize again and
repeat the diagnostic. Introduced ammonium and native terminal states are held;
missing native contexts and modified-product states are reported explicitly.

Every structure export includes all charges, masses, bonded/LJ/1-4 and CMAP
parameters in numeric JSON/XML, a full-system OFF library, ordinary Amber
frcmod and a type mapping. See [FORCE_FIELDS.md](FORCE_FIELDS.md) for the
profile transfer, complete rebuild and engine compatibility checks.

## Software references

- [OpenMM force fields and implicit-solvent models](https://docs.openmm.org/latest/userguide/application/02_running_sims.html).
- [APBS documentation](https://apbs.readthedocs.io/en/latest/index.html) and
  [continuum-solvation background](https://apbs.readthedocs.io/en/latest/background.html).

Gaussian/RESP development is separate from the installed Python runtime.
The repository does not redistribute Gaussian or standard Amber library
files.
