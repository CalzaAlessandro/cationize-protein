# HSA cationization example

This folder contains the existing **585-residue HSA series** at five
cationization levels. The structures were minimized with the APBS-guided
protocol and renamed using the current residue names. They were produced in the HSA campaign before this CLI release.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/hsa_workflow_dark.gif">
  <img src="../../docs/images/hsa_workflow.gif" width="100%" alt="Four molecular illustrations in a closed cycle beside the saved minimized HSA structures from 0 to 100 percent cationization">
</picture>

## Files

```text
HSA/
├── 1n5u_native_ready.pdb
├── 0/
├── 25/
│   ├── minimized.pdb
│   ├── protein.parm7
│   └── minimized.rst7
├── 50/
├── 75/
├── 100/
├── residue_templates/
├── series_summary.csv
├── modification_sites.csv
└── PROVENANCE.json
```

Each percentage has the same three structure files. The PDB is ready for
inspection. Pair `protein.parm7` with `minimized.rst7` from the **same folder**
to load the parameterized protein. These files describe dry, minimized
structures; they contain no production trajectories.

`residue_templates/` contains the individual custom `.lib` files in Amber OFF
format. `modification_sites.csv` lists every modified residue and identifies
sites added since the preceding checkpoint. `PROVENANCE.json` records source
artifacts, hashes and validation checks without linking to local source paths.

## Series summary

| Cationization | Modified residues | Actual fraction of 99 groups | Protein charge, e | Disulfides |
|---|---:|---:|---:|---:|
| 0% | 0 | 0.00% | −14 | 17 |
| 25% | 25 | 25.25% | +36 | 17 |
| 50% | 50 | 50.51% | +85 | 17 |
| 75% | 74 | 74.75% | +133 | 17 |
| 100% | 99 | 100.00% | +183 | 17 |

Charges are rounded sums of topology partial charges and include the native
protonation states used in the campaign. The topology/PDB atom order,
coordinate agreement, cumulative modification sets and 17 disulfides were
checked at all five levels. Modified-residue charges agree with the distributed
OFF libraries within Amber topology printing precision.

## Residue names

| Name | Original group | Context |
|---|---|---|
| ASC | Asp side-chain carboxyl | Internal residue |
| GLC | Glu side-chain carboxyl | Internal residue |
| ASX | Asp side-chain carboxyl | N-terminal residue |
| LET | C-terminal carboxyl | C-terminal Leu in this HSA model |

The native N-terminal amine is retained in ASX. LET modifies the terminal
carboxyl group; it does not modify the Leu side chain. The complete package
also includes templates for other terminal residue identities.

## Start a new calculation

With the environment activated and from the repository root:

```sh
cationize-protein prepare examples/HSA/1n5u_native_ready.pdb \
  --output hsa_prepared --ph 7.4
cationize-protein run --pdb hsa_prepared/protein.pdb --output hsa_new_series \
  --solver apbs --force-field ff14SB --formats amber openmm gromacs
```

This creates a separate series using the current implementation and initial
histidine assignment. APBS must be installed for this command. It does not
modify the example files or launch MD.

## Figure

All structures are fitted on the same 585 Cα atoms and displayed with one
camera and scale. Protein cartoons are colored by secondary structure;
modified residues appear as heavy-atom sticks. Native backbone aliases are
used only in visualization copies so VMD can recognize the modified residues.
The original example PDB names remain unchanged.

The figure shows five minimized structures at increasing cationization levels.
The four cycle insets illustrate separate operations. The electrostatic and
protonation examples use native HSA. The amidation and minimization examples
use an Asp–Asp–Leu peptide with the distributed ASC parameters.

### Electrostatic ranking

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_potential_dark.png">
  <img src="../../docs/images/cycle_step_potential.png" width="650" alt="Native HSA surface colored by screened-Coulomb electrostatic potential, from minus 2 to plus 2 kBT per elementary charge">
</picture>

The surface is colored by a screened-Coulomb calculation using the native
topology charges, a solvent dielectric of 78.5, 0.15 M monovalent salt and
298 K. Red indicates negative potential and blue positive potential. The
display uses a 3 Å grid and a fixed range of −2 to +2 kBT/e. This illustrates
the default solver; the saved HSA cationization series was guided by APBS.

### Amidation and rotation

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_amidation_dark.png">
  <img src="../../docs/images/cycle_step_amidation.png" width="650" alt="Asp side-chain fragment before modification and the corresponding ASC fragment with its aminoethylamide group">
</picture>

The Asp side-chain carboxyl is converted into an aminoethylamide. The views
show the same protein-bound site before modification and after local
relaxation. The script samples rotations of the added group to reduce steric
overlap before minimization. Hydrogen atoms are omitted in this view.

### Energy minimization

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_minimization_dark.png">
  <img src="../../docs/images/cycle_step_minimization.png" width="650" alt="Overlay of the actual modified peptide side-chain coordinates before and after restrained OpenMM energy minimization">
</picture>

The overlay compares the modified side chain immediately before minimization
(orange) and after convergence (green), fitted on the local backbone. The
calculation used ff14SB, the custom RESP parameters and OpenMM GBn2 with
backbone restraints. Both geometries were evaluated with the same energy
function; the energy decreased by 35.4 kJ/mol.

### PROPKA protonation

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_protonation_dark.png">
  <img src="../../docs/images/cycle_step_protonation.png" width="650" alt="Native HSA His3 in the neutral HIE state, with a PROPKA pKa of 5.78 at pH 7.4">
</picture>

PROPKA predicts a pKa of 5.78 for His3 in the native HSA example. At pH 7.4,
the neutral HIE state is retained. The introduced terminal amines remain
fixed at +1 under the cationization protocol.

Rendering, source hashes and calculation details are recorded in
[HSA_PREVIEW.json](../../docs/images/HSA_PREVIEW.json).
