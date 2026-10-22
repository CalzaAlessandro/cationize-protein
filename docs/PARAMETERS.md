# Residue parameters

The library describes carboxyl amidation with ethylenediamine. The product is
`–C(=O)–NH–CH₂–CH₂–NH₃⁺`, with the introduced ammonium fixed at +1.
The complete residue charge also includes its other ionizable groups and
terminal charges. Separate templates describe internal, N-terminal and
C-terminal contexts.

## What is included

The [manifest](../parameters/manifest.json) lists 32 chemical states: Asp/Glu
side-chain products, their terminal variants, and modified C termini for the
standard amino acids. Histidine tautomers, protonated histidine, disulfide-linked
cysteine and thiolate cysteine have separate entries. C-terminal Asp and Glu
can contain modifications on the side chain, the terminus, or both groups.

Four accepted HSA templates are retained under their current names: `ASC`,
`GLC`, `ASX` and `LET`. Their fitted charges are unchanged. The remaining
28 entries were fitted with the capped-model procedure below. The
[numerical validation table](../parameters/validation_summary.csv) gives
the charge and fit status of every entry.

## How the new templates were built

### 1. Build a capped amino-acid model

Native atoms, geometry and atom types come from the installed Amber residue
definitions. The required carboxyl group is converted into the amide product.
The cap arrangement preserves the chemical context of the residue:

| Residue context | Quantum model |
|---|---|
| Internal residue | ACE–modified residue–NME |
| N-terminal residue | Modified residue–NME |
| C-terminal residue | ACE–modified residue |
| Disulfide-linked cysteine | An explicit disulfide partner with ACE/NME caps |

ACE and NME represent the neighboring peptide groups. Their standard charges
remain fixed during fitting. The disulfide partner also retains its native
charges. The fitted product therefore joins a protein with the intended
terminal state, peptide connections and disulfide chemistry.

### 2. Optimize geometry and calculate the electrostatic potential

Two starting backbone conformations, α-like and extended, are optimized with
[Gaussian](https://gaussian.com/) at **HF/6-31G(d)**. A backbone torsion and
heteroatom–hydrogen bond lengths are held fixed to retain the chosen
conformation and protonation state. Merz–Kollman electrostatic-potential samples
are calculated for the resulting geometries.

### 3. Fit the partial charges with RESP

The two conformers enter a joint, two-stage **RESP** fit using AmberTools.
Corresponding atoms share their charges between conformers; chemically
equivalent atoms are constrained together. Cap and disulfide-partner charges
remain fixed. The modified residue must recover its specified integer net
charge within numerical rounding.

The RESP charges describe the modified residue. During cationization, PROPKA
selects among supported native protonation states; it does not refit these
partial charges or titrate the introduced ammonium.

### 4. Assign protein force-field terms

The reacted carbonyl receives Amber amide C/O types. The new amide, methylene
and ammonium atoms use protein atom types. Bond, angle, torsion and
Lennard-Jones parameters are transferred by chemical analogy from the protein
force field. **No new torsional-energy fit was performed.**

The default profile is ff14SB. The ff19SB and ff99SBildn profiles retain the
same custom RESP charges and use their corresponding native parent-residue
contexts. ff19SB also retains its backbone CMAP terms. Profile assignment and
current coverage are described in [Force fields](FORCE_FIELDS.md).

### 5. Check the exported residue

Acceptance checks cover Gaussian optimization and ESP completion, bond
connectivity, stereochemistry, protonation, cap charges, equivalent atoms,
integer residue charge and the fit residuals. The capped model must have
complete Amber parameters and finite OpenMM energies and forces in vacuum and
GBn2. Export checks compare the atom graph, charges and physical parameters
after format conversion.

These checks establish consistency of the fitted models and their exports.
Additional peptide simulations, torsion benchmarks and experimental comparisons
would be needed to assess transferability to a new chemical environment.
Detailed criteria and completed calculations are in [Methods](METHODS.md) and
[Validation](VALIDATION.md).

## Where to find the parameters

| Location | Contents |
|---|---|
| `parameters/NAME.lib` | One Amber OFF residue template with atom names, types, charges and connectivity |
| `parameters/NAME.json` | Atom mapping and fitted model data for newly fitted entries |
| `parameters/manifest.json` | Parent residue, terminal context, modified groups, charge, status and file hashes |
| `parameters/validation_summary.csv` | Numerical charge and fit checks |
| `parameters/cationized_residues.xml` | The ff14SB OpenMM residue extension |
| `OUTPUT/residue_templates/` | Verified copies of the library used for a cationization run |
| `OUTPUT/PERCENT/exports/` | Complete parameters for that generated protein: OFF/frcmod, XML/JSON and selected engine files |

The single-residue OFF files use the original Amber atom-type names. The
complete protein export has its own type mapping and printed numerical
parameters. For a ligand complex, follow the [assembly tutorial](COMPLEXES.md).

## Extend or refit the library

Work in a separate development directory and retain the accepted library.
Gaussian and AmberTools (`resp`, `respgen`, `espgen`, `antechamber` and `tleap`)
are needed for this development procedure. RDKit prepares the capped models.
These tools are separate from the dependencies used to run cationization with
the bundled templates.

With the project environment activated, install the Python parameter-development
extra from the repository directory:

```sh
python -m pip install '.[parameterization]'
```

This installs RDKit. Install Gaussian and AmberTools separately, set `AMBERHOME`
to the Amber installation and make their executables available on `PATH`.

For an existing catalog state, the development sequence is:

1. Prepare its capped models and atom mapping with
   [`models.prepare`](../src/cationize_protein/models.py), using the installed
   Amber definitions. Create a development manifest containing the selected
   chemical states. This supplies `work/models/NAME/alpha/`,
   `work/models/NAME/extended/` and `atom_mapping.json`.
2. Run the prepared fits sequentially. For example, for an N-terminal Glu
   side-chain product:

   ```sh
   cationize-protein parameterize --project parameter-development --names GLX
   ```

   The project must already contain its model inputs and manifest. This command
   runs the Gaussian calculations, RESP fitting and exported-template checks.
   The current worker discovers Gaussian as `g16` on `PATH`; it does not submit
   Slurm jobs.
3. Review the geometries, ESP residuals and charge/energy checks. Retain the
   Gaussian inputs/logs, fitted charges, atom mapping and acceptance records.
4. Add the accepted entry to a **copy of the complete library**, updating its
   manifest, hashes, XML and validation tables together. Check it in a short
   peptide and verify the requested engine exports before using it in a protein.
5. Select the complete custom library when running cationization:

   ```sh
   cationize-protein run --pdb prepared/protein.pdb --output custom-series --parameters custom-parameters
   ```

Adding a chemical product beyond the present catalog also requires extending
the [chemical-state catalog](../src/cationize_protein/catalog.py),
[capped-model builder](../src/cationize_protein/models.py), atom typing and
the corresponding structure/protonation transitions. A new three-letter name
alone does not define those operations. Terminal states and disulfide contexts
need their own mapped templates. New force-field families would require an
independent parameterization and validation procedure.
