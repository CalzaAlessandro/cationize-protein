# HSA cationization example

This folder contains the existing **585-residue HSA series** at five
cationization levels. The structures were minimized with the APBS-guided
protocol and renamed using the current residue names. They were produced in the HSA campaign before this CLI release.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/hsa_workflow_dark.gif">
  <img src="../../docs/images/hsa_workflow.gif" width="100%" alt="A text-only cycle at left and a molecular slideshow at right, showing the saved HSA structure and step images separately for each percentage">
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

The left panel contains the four text labels and highlights the active phase.
The right panel shows each percentage in turn: the saved HSA structure,
its electrostatic surface, a chemical-step example, the corresponding local
minimization and the final PROPKA check. Native HSA has no amidation page.

The chemical examples are the last sites added before each percentage
checkpoint. Their coordinates come from the original HSA campaign, including
the saved starting placement and local minimization. These local steps
preceded the final solvated minimization of the complete checkpoint.

| Checkpoint | Chemical step | Site modified | Product |
|---|---:|---|---|
| 25% | 25/99 | Leu585, C-terminal carboxyl | LET |
| 50% | 50/99 | Glu376, side-chain carboxyl | GLC |
| 75% | 74/99 | Glu442, side-chain carboxyl | GLC |
| 100% | 99/99 | Asp249, side-chain carboxyl | ASC |

### Electrostatic ranking

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_potential_dark.png">
  <img src="../../docs/images/cycle_step_potential.png" width="650" alt="Native HSA surface colored by screened-Coulomb electrostatic potential, from minus 5 to plus 5 kBT per elementary charge">
</picture>

The surface is colored by a screened-Coulomb calculation using the native
topology charges, a solvent dielectric of 78.5, 0.15 M monovalent salt and
298 K. Red indicates negative potential and blue positive potential. The
display uses the same 3 Å grid and fixed range of −5 to +5 kBT/e for all five
percentages. Each field is calculated from its own checkpoint coordinates
and topology charges. This illustrates the default solver; the saved HSA
cationization series was guided by APBS. Native HSA is shown above; the
animation includes all five surfaces.

### Amidation and rotation

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_amidation_dark.png">
  <img src="../../docs/images/cycle_step_amidation.png" width="650" alt="HSA Asp249 before modification and after conversion to ASC in chemical step 99, reaching the 100 percent checkpoint">
</picture>

The example shows Asp249 before modification and after conversion to ASC
at step 99. The views use the saved HSA coordinates, fitted on the local
fragment. Rotations of the added group reduce steric overlap before
minimization. Hydrogen atoms are omitted. The animation uses the respective
sites listed in the table for the other checkpoints.

### Energy minimization

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_minimization_dark.png">
  <img src="../../docs/images/cycle_step_minimization.png" width="650" alt="HSA ASC249 immediately before and after its actual restrained Amber local minimization at step 99">
</picture>

The overlay compares the modified side chain immediately before minimization
(orange) and after convergence (green), fitted on the local backbone. The
archived HSA calculation used `pmemd.cuda`, ff14SB, the custom RESP parameters,
implicit solvent (`igb=5`, 0.15 M salt) and backbone restraints of
10 kcal/(mol Å²). At step 99, the starting geometry is the saved rotamer
correction used before the successful minimization. Displacements are shown
at their actual scale. Native HSA is displayed as its minimized reference.

### PROPKA protonation

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../docs/images/cycle_step_protonation_dark.png">
  <img src="../../docs/images/cycle_step_protonation.png" width="650" alt="His3 in 100 percent cationized HSA, with its saved neutral HIE geometry and original PROPKA pKa of 5.93 at pH 7.4">
</picture>

Each frame uses the His3 geometry and archived PROPKA result for that
percentage. Predicted pKa values are 5.78, 6.02, 5.96, 5.53 and 5.93 at
0%, 25%, 50%, 75% and 100%, respectively. The neutral HIE state is retained
at pH 7.4. Introduced terminal amines remain fixed at +1. The enlarged
example above shows the 100% checkpoint.

Rendering, source hashes and calculation details are recorded in
[HSA_PREVIEW.json](../../docs/images/HSA_PREVIEW.json).
