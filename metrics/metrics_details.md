# Metrics Details: Formulas and Script Calculations

This document lists all metrics computed by the FoldKit analyser and metrics scripts, with:

**Related (not covered below):** `trim_models.py` only harmonises residue ranges (standalone script; shared with `trim_superimposeLSQ.py` trimming); it does not compute packing or similarity metrics. See **README.md** (File management). **Post-processing CSV tools** `interface_molecule_report_csv.py` and `contact_molecule_report_csv.py` do not compute new metrics; they parse analyser text reports into tables (Section 2.10, Section 3.5). Example filenames in this repo use `model_01`, `results.txt` (interface / packing merged output), `contact_results.txt` (`contact_analyser.py -o`), `set_a` / `set_b` (`--sets`), and `./out`; see **README.md** (opening “Example naming” paragraph).

- **Mathematical definition** and **script calculation** for each metric.
- **Reader-friendly formula**: a short, plain-language description of what is being computed and how.
- **Libraries and functions**: BioPython, NumPy, SciPy (optional), pandas (where used), **biotite** (optional in Dali workflows), as summarised in the table below.
- **Derivation / references** where applicable.

---

## Libraries used (summary)


| Script                  | Libraries                                                                                                                        | Main functions used                                                                                                                                                                                                                                                                                                                                                                                    |
| ----------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| packing_metrics.py      | **BioPython** (Bio.PDB.PDBParser, PDBIO, PDBExceptions), **NumPy**, **SciPy** (optional)                                         | `PDBParser.get_structure()`, manual CRYST1 read; `structure` → model → chain → residue → atom; `atom.coord`, `atom.element`, `atom.name`; `np.radians`, `np.sqrt`, `np.cos`; `scipy.spatial.distance.pdist()` or fallback loop with `np.linalg.norm`                                                                                                                                                   |
| interface_analyser_asu_charge.py / interface_analyser_lattice_charge.py | **BioPython** (PDBParser, NeighborSearch), **Bio.PDB.SASA** (ShrakeRupley, optional), **NumPy** | `parser.get_structure()`, `structure.get_chains()`, `chain.get_atoms()`; `NeighborSearch(atom_list).search_all(threshold, 'R')` → interface residues; contact residues = interface residues filtered by same criteria (distance + type) via `_residue_passes_contact_criteria()`; `ShrakeRupley().compute()` for BSA; `np.linalg.norm`, `np.mean`, etc. Includes shape complementarity and **charge-tag** complementarity metrics. Text report → CSV: Section 2.10. |
| interface_analyser_asu_ec.py / interface_analyser_lattice_ec.py / electrostatic_complementarity.py | **BioPython** (PDBParser, NeighborSearch), **NumPy** | McCoy-style electrostatic complementarity (EC): sample facing surface points and compute electrostatic potentials; per-interface Pearson correlation \(r\) and lattice-weighted summaries (Fisher z). Definitions: `EC_details.md`; summary integration: Section 2.5.4 and Section 2.7.1. |
| contact_analyser.py     | **BioPython** (PDBParser, PDBExceptions), **NumPy**                                                                              | Same structure traversal; `chain.get_atoms()`, `atom.coord`; `np.linalg.norm`; manual contact loop and type classification; optional full ASU contact list in text / `*_asu_contacts.txt`. Text report → CSV: Section 3.5.                                                                                                                                                                                                                                                                             |
| dali_score.py           | **BioPython** (PDBParser), **NumPy**, **biotite** (optional)                                                                     | `PDBParser.get_structure()`, structure traversal for CA atoms; `np.linalg.norm`, `np.exp` for distance matrices and Dali formula; `biotite.structure.superimpose_structural_homologs()` for automatic residue equivalences when available                                                                                                                                                              |


---

## 1. packing_metrics.py

**Function selection / usage**

- These metrics are computed in `packing_metrics.py`, typically called by `crystal_packing_analyser.py` for each structure.
- At the CLI level, you select whether to run packing metrics by:
  - Running `packing_metrics.py` directly on one or more PDB/CIF files, or
  - Enabling the “packing” stage in `crystal_packing_analyser.py` (it runs packing by default unless explicitly skipped in that orchestrator).
- Within `packing_metrics.py`, the main driver function `calculate_metrics()` calls:
  - `_calculate_unit_cell_volume()` → Section 1.1
  - `_analyse_molecular_content()` → Section 1.2
  - `_calculate_matthews_metrics()` → Section 1.3
  - `_calculate_packing_density()` and `_analyse_void_spaces()` → Sections 1.4–1.5

### 1.1 Unit cell volume

**Reader-friendly formula**

*Unit cell volume* is the volume of the crystallographic box. For a general box with side lengths *a*, *b*, *c* and angles α, β, γ between them, the volume is: *(a × b × c)* times a geometric factor that depends only on the three angles. The script reads *a*, *b*, *c*, α, β, γ from the PDB `CRYST1` line, converts angles to radians, then applies this formula.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Larger unit cell volume | More physical space per unit cell (often due to larger macromolecule(s), more solvent, higher Z, or looser packing). |
| Smaller unit cell volume | Less space per unit cell (often tighter packing and/or smaller contents). |
| Unusual values / defaults | If `CRYST1` is missing, defaults are used (a=b=c=1 Å, α=β=γ=90°), so the reported volume is not meaningful. |

**Mathematical definition**

For a general parallelepiped with edge lengths a, b, c (Å) and angles alpha, beta, gamma (degrees):

```
V = a*b*c*sqrt(1 + 2*cos(alpha)*cos(beta)*cos(gamma)
               - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2)
```

**Implementation (simplified):** Read a, b, c (Å) and angles α, β, γ (degrees) from the PDB CRYST1 line. Convert angles to radians. Volume = a × b × c × sqrt(1 + 2×cos(α)×cos(β)×cos(γ) − cos²(α) − cos²(β) − cos²(γ)).

**Script calculation**

- Unit cell parameters are read from the PDB `CRYST1` line (defaults: a=b=c=1, α=β=γ=90°, space_group P1 if missing).
- Angles are converted to radians, then:
  - `volume = a * b * c * np.sqrt(1 + 2*cos(α)*cos(β)*cos(γ) - cos(α)**2 - cos(β)**2 - cos(γ)**2)`
- Implemented in `_calculate_unit_cell_volume()` (and similarly in `contact_analyser._calculate_cell_volume()`).

**Libraries / functions used**

- **Input:** Unit cell is read manually by parsing the file for a line starting with `CRYST1` and slicing fields (no BioPython for CRYST1). Structure is loaded with **Bio.PDB.PDBParser** → `parser.get_structure('crystal', pdb_file)` for use in later steps.
- **NumPy:** `np.radians()` to convert angles, `np.cos()`, `np.sqrt()` for the volume expression.

---

### 1.2 Molecular content

**Reader-friendly formula**  
*Molecular content* is a set of counts and one sum: how many atoms, residues, and chains; molecular weight (sum of atomic masses over all atoms); and how many atoms per element. The script walks the structure model → chain → residue → atom and uses each atom’s element (or name) to look up mass and to update counts.

**What it means (how to interpret)**

| Metric | Value pattern | Typical interpretation |
| --- | --- | --- |
| `total_atoms`, `residue_count`, `chain_count` | Higher | Larger system and/or more molecules/chains in the model. |
| `molecular_weight` | Higher | Heavier contents (more atoms and/or heavier elements); depends on what’s present in the structure file. |
| `element_composition` | Unexpected elements / unusual proportions | Can indicate ligands/metals/ions, or incorrect/missing element assignment in the input. |

**Implementation (simplified):** Loop over all models → chains → residues → atoms. For each atom: increment total_atoms; add atomic_weights[element] to molecular_weight (use 12.011 for unknown); increment element_composition[element]. Count chains and residues. Result: total_atoms, molecular_weight, element_composition, chain_count, residue_count.


| Metric                  | Mathematical definition                         | Script calculation                                                      |
| ----------------------- | ----------------------------------------------- | ----------------------------------------------------------------------- |
| **total_atoms**         | Count of all atoms in the structure             | Loop over all models → chains → residues → atoms; increment counter.    |
| **molecular_weight**    | M = sum_over_atoms( m_element ) (Da) | Sum of `atomic_weights[element]` per atom (default 12.011 for unknown). |
| **element_composition** | Count per element                               | Dict of element → count over all atoms.                                 |
| **chain_count**         | Number of chains                                | Count of chains across all models.                                      |
| **residue_count**       | Number of residues                              | Count of residues across all models.                                    |


- Implemented in `_analyse_molecular_content()`.
- Atomic weights (Da): C 12.011, N 14.007, O 15.999, S 32.06, P 30.974, H 1.008, etc.

**Libraries / functions used**

- **BioPython:** `parser.get_structure()` returns a `Structure`; iteration `for model in structure`, `for chain in model`, `for residue in chain`, `for atom in residue`. Per atom: `atom.element` (string), `atom.name` (fallback if element missing). No PDBIO/Select used for this metric.
- **Script:** Built-in dict `atomic_weights` keyed by element symbol; default carbon mass for unknown.

---

### 1.3 Matthews coefficient and solvent content

**Reader-friendly formula**  
*Matthews coefficient* is unit-cell volume per dalton of protein (Å³/Da). *Protein volume* is estimated as mass divided by a fixed density (≈0.81 Da/Å³). *Solvent volume* is unit-cell volume minus protein volume. *Solvent content* is that solvent volume as a percentage of unit-cell volume.

---

### 1.6 Lattice-wide packing metrics: lattice_packing_analyser.py

**Purpose.** `lattice_packing_analyser.py` computes packing-style metrics for a *symmetry-expanded coordinate model* (a supercell / lattice fragment) where many copies are present in one Cartesian frame (typically P1). Unlike `packing_metrics.py`, which treats the input as a single crystal unit cell, this script is intended for multi-copy PDBs such as `supercell.pdb`.

**Important limitation (boundary effects).** The reported values depend on how the supercell was constructed (number of copies, spacing, and how much “empty” space exists at the boundary). For meaningful comparisons, use the **same expansion protocol** across structures.

#### 1.6.1 Volume definition (required for densities)

The script supports two volume sources:

- **CRYST1 volume** (`--volume-source cryst1`): parse the PDB `CRYST1` record and compute a unit-cell volume \(V\) using the standard parallelepiped formula (same as Section 1.1). This is appropriate when your supercell PDB includes a meaningful `CRYST1` record.
- **Bounding-box volume** (`--volume-source bbox`): compute the axis-aligned coordinate bounding box over all atoms, optionally padded by `--bbox-pad` Å. Let \(\Delta x, \Delta y, \Delta z\) be the box side lengths; then:

```
V_bbox = Δx * Δy * Δz
```

Default (`--volume-source auto`): use CRYST1 if present, otherwise bounding box.

#### 1.6.2 Atom density and mass density (whole supercell)

Let \(N_{\mathrm{atoms}}\) be the number of atoms in the file and \(M\) the total mass estimated by summing atomic weights by element (Da). For the chosen volume \(V\) (Å³):

```
atom_density = N_atoms / V                (atoms/Å^3)
mass_density = M / V                      (Da/Å^3)
```

**Interpretation.**
- Higher densities indicate a more “filled” coordinate volume, but the absolute number is sensitive to the chosen volume definition (CRYST1 vs bounding box) and supercell boundaries.

#### 1.6.3 Packing density (fraction and percent)

The script estimates a total “occupied” volume by summing per-atom tabulated volumes \(v_i\) (Å³) by element:

```
V_atoms = sum_i(v_i)
packing_density_fraction = V_atoms / V
packing_density_percent = 100 * packing_density_fraction
```

**Interpretation.**
- This is a geometric proxy for packing fraction. Values near 0 indicate mostly empty volume (or a very large bounding box); values closer to 1 indicate a tightly filled volume. In protein crystals, this number is strongly method-dependent because \(v_i\) are approximate and the coordinate volume is user-defined in supercell mode.

#### 1.6.4 Matthews-like ratio and solvent content (heuristic)

The script reports a Matthews-like ratio for the whole supercell:

```
lattice_matthews_a3_per_da = V / M
```

It also reports a heuristic “solvent content” based on the same constant used in `packing_metrics.py` (protein density ≈ 0.81 Da/Å³):

```
V_protein = M / 0.81
V_solvent = V - V_protein
estimated_solvent_content_fraction = V_solvent / V
estimated_solvent_content_percent = 100 * estimated_solvent_content_fraction
```

**Interpretation and caveats.**
- These are convenient summary ratios for comparing similarly-built supercells.
- They are **not** equivalent to crystallographic solvent content unless the CRYST1 unit cell is meaningful for the model and the contents correspond to a true crystallographic asymmetric unit/packing.


**What it means (how to interpret)**

| Metric | Value pattern | Typical interpretation |
| --- | --- | --- |
| Matthews coefficient (Å³/Da) | Higher | More space per mass → typically higher solvent fraction / looser packing. |
| Matthews coefficient (Å³/Da) | Lower | Less space per mass → typically tighter packing / lower solvent fraction (may be unrealistic if extremely low). |
| Solvent content (%) | Higher | More solvent in the unit cell; often easier diffusion but potentially more disorder. |
| Solvent content (%) | Lower | Less solvent; often tighter crystals (but very low values can be suspect). |

**Mathematical definition**

**Matthews coefficient** (Å³/Da):

```
Vm = Vcell / M
```

**Estimated protein volume** (Å³):

```
Vprotein = M / rho,   rho ≈ 0.81 Da/Å^3
```

**Solvent volume** (Å³):

```
Vsolvent = Vcell - Vprotein
```

**Solvent content** (%):

```
solvent_percent = (Vsolvent / Vcell) * 100
```

**Implementation (simplified):** Matthews = unit_cell_volume / molecular_weight (0 if molecular_weight ≤ 0). Protein_volume = molecular_weight / 0.81 (Å³). Solvent_volume = unit_cell_volume − protein_volume. Solvent_content = (solvent_volume / unit_cell_volume) × 100, clamped to [0, 100].

- **Script:** In `_calculate_matthews_metrics()`.

**Libraries / functions used**

- **Script only:** Arithmetic using `unit_cell_volume` and `molecular_weight` from previous steps; no extra BioPython or NumPy calls in this function.

**Derivation / references**

- The **Matthews coefficient** is defined as the ratio of unit-cell volume to protein mass:
  - B. W. Matthews, “Solvent content of protein crystals”, *J. Mol. Biol.* **33**, 491–497 (1968). DOI: [10.1016/0022-2836(68)90205-2](https://doi.org/10.1016/0022-2836(68)90205-2)
- The typical protein partial specific volume (≈0.73–0.75 mL/g) and average crystal solvent content (40–60%) imply a bulk density on the order of **1.35 g/cm³** (≈ **0.81 Da/Å³**), which motivates the protein-volume estimate `V_protein` used in the script.
- The solvent fraction follows from a volume balance: `V_cell = V_protein + V_solvent`, so `solvent_fraction = V_solvent / V_cell`.

---

### 1.4 Packing density and void fraction

**Reader-friendly formula**  
*Total atomic volume* is the sum of a fixed volume per atom type (from van der Waals–style tables). *Packing density* is that sum divided by unit-cell volume (what fraction of the cell is “filled” by atoms). *Packing efficiency* is the same as packing density expressed as a percentage. *Void fraction* is one minus packing density (the unfilled fraction).

**What it means (how to interpret)**

| Metric | Value pattern | Typical interpretation |
| --- | --- | --- |
| Packing density | Higher | Unit cell is more “filled” by atom volumes → tighter packing (lower voids). |
| Packing density | Lower | More empty space relative to atom volumes → looser packing / more void. |
| Void fraction | Higher | More unfilled fraction of the unit cell (inverse of packing density). |
| Void fraction | Lower | Less unfilled space (tighter packing). |

**Mathematical definition**

**Total atomic volume** (Å³):

```
V_atoms = sum_over_atoms( V_element )
```

(equivalently, sum_over_atoms(...); the subscript here is descriptive, not a variable name.)

**Packing density** (dimensionless):

```
packing_density = V_atoms / V_cell
```

**Packing efficiency** (%):

```
packing_efficiency_percent = packing_density * 100
```

**Void fraction:**

```
void_fraction = 1 - packing_density
```

**Implementation (simplified):** For each atom, look up atomic_volumes[element] (default 16.44 Å³ for C). Total_atomic_volume = sum over all atoms. Packing_density = total_atomic_volume / unit_cell_volume (0 if volume ≤ 0). Packing_efficiency_percent = packing_density × 100. Void_fraction = 1 − packing_density.

**Libraries / functions used**

- **BioPython:** Same structure traversal (model → chain → residue → atom); `atom.element`, `atom.name` for element lookup.
- **Script:** Dict `atomic_volumes` by element; default carbon volume.
- **NumPy:** Used only for results (no np calls in this block beyond possible type handling).

**Derivation / references**

- **Atomic volumes** are approximated from empirical van der Waals volumes, e.g.:
  - A. Bondi, “van der Waals Volumes and Radii”, *J. Phys. Chem.* **68**, 441–451 (1964). DOI: [10.1021/j100785a001](https://doi.org/10.1021/j100785a001); and later tabulations used in protein packing analyses.
- **Packing density** is the fractional occupancy of the unit cell by atom-based spheres:
  - `V_atoms = sum of V_element` over atoms, so `phi = V_atoms / V_cell`.
  - `phi` close to 0.7–0.8 is typical for tightly packed proteins; higher `void_fraction` indicates looser packing.

---

### 1.5 Void / interatomic metrics (_analyse_void_spaces)

**Reader-friendly formula**  
Each atom is assigned a *radius* from its element’s volume assuming a sphere (V = ⁴⁄₃πr³ → r = (3V/(4π))^(1/3)). *Average interatomic distance* is the mean of all pairwise distances between atoms; *minimum* is the smallest such distance. *Atom density* is number of atoms divided by unit-cell volume. The script collects all atom coordinates, then either uses SciPy’s `pdist` for pairwise distances or a double loop with `np.linalg.norm`.

**What it means (how to interpret)**

| Metric | Value pattern | Typical interpretation |
| --- | --- | --- |
| Average interatomic distance | Lower | More tightly packed atoms on average (often correlates with tighter packing). |
| Average interatomic distance | Higher | Looser average packing (more space between atoms). |
| Minimum interatomic distance | Very low | Potential steric clashes / overlaps (or artefacts/alternate conformations); interpret cautiously. |
| Atom density (atoms/Å³) | Higher | More atoms per volume → tighter packing and/or larger contents relative to cell volume. |
| Atom density (atoms/Å³) | Lower | Fewer atoms per volume → looser packing and/or more solvent space. |

**Mathematical definition**

**Atomic radius from volume** (sphere model):

```
V = (4/3) * pi * r^3
r = (3*V/(4*pi))^(1/3)
```

**Average interatomic distance** (Å):

```
d_bar = (1/|P|) * sum_over_pairs_(i,j in P)( ||x_i - x_j|| )
```

where P is the set of all unordered atom pairs.

**Minimum interatomic distance** (Å):

```
d_min = min_over_pairs_(i,j)( ||x_i - x_j|| )
```

**Atom density** (atoms/Å³):

```
rho_atom = N_atoms / V_cell
```

**Implementation (simplified):** Per atom, radius = (3 × atomic_volumes[element] / (4π))^(1/3). Collect all atom coordinates; compute all pairwise distances (SciPy pdist or double loop). Average interatomic distance = mean(distances); minimum = min(distances). Atom density = number of atoms / unit_cell_volume (0 if no atoms or volume ≤ 0).

**Libraries / functions used**

- **BioPython:** Structure iteration; `atom.coord` (NumPy-like array), `atom.element`, `atom.name`.
- **NumPy:** `np.array(coords)`, `np.pi`; `np.mean(distances)`, `np.min(distances)`; fallback: `np.linalg.norm(coords[i] - coords[j])` in a double loop.
- **SciPy (optional):** `scipy.spatial.distance.pdist(coords)` to compute all pairwise distances in one call; script falls back to the loop if SciPy is not available.

**Derivation / references**

- The atomic radius from volume assumes a **spherical atom model**: `V = (4/3) * pi * r^3` → `r = (3*V/(4*pi))^(1/3)`.
- Average and minimum interatomic distances are standard descriptors of packing tightness and are analogous to nearest-neighbour statistics in condensed-matter physics.
- Atom-density measures mirror typical crystallographic discussions of packing fraction and nearest-neighbour spacing in protein crystals.

---

## 2. Interface analysis (charge-tag and EC)

**Function selection / usage**

- Interface metrics are computed via the ASU and lattice entrypoints:
  - `interface_analyser_asu_charge.py` / `interface_analyser_lattice_charge.py` (charge-tag metrics)
  - `interface_analyser_asu_ec.py` / `interface_analyser_lattice_ec.py` (McCoy EC)

  All variants share the core implementation in `interface_analyser_base.py`.
- At the CLI level, you select these metrics by:
  - Running one of the ASU or lattice interface entrypoints directly on a PDB/CIF (specifying chains or letting the script enumerate interfaces), or
  - Allowing `crystal_packing_analyser.py` to run the interface stage on each structure.
- Within the interface analysers, the main pipeline:
  - **Interface residues:** `NeighborSearch(atom_list).search_all(contact_distance, 'R')` → residue pairs within threshold, different chains.
  - **Contact residues:** Those interface residues filtered by applying the same criteria (distance + type-specific limits) as the contact list: keep residue iff it has ≥1 atom that forms a pair with the other chain satisfying those criteria.
  - Calls `_find_contacts()` and `_filter_contacts_by_limits()` → Sections 2.1–2.4, 2.6.
  - **BSA:** `_calculate_buried_surface_area(chain1, chain2, model)` using `Bio.PDB.SASA.ShrakeRupley` when available (Section 2.2); else fallback estimate.
  - Calls `_calculate_contact_area()` and `_calculate_interface_complementarity()` → Sections 2.3–2.5, 2.7.
- **Contact distance:** Default 5.0 Å; atom pairs with distance ≤ this are candidates; only those passing type-specific limits (below) are counted.
- **Multi-copy / symmetry-expanded models:** When the input file contains many chains representing repeated molecules in a lattice fragment, the same pipeline reports all pairwise interfaces present in the file and optional focal-copy metrics (`reference_chain_id` / `--reference-chain`); see Section 2.7.2.

**Type-specific distance limits (contacts outside these are ignored):**


| Type          | Limit (Å) |
| ------------- | --------- |
| H-bonds       | ≤ 3.5     |
| Electrostatic | ≤ 5.0     |
| Hydrophobic   | 3.5–4.5   |
| Van der Waals | ≤ 5.0     |


- Contacts are classified by atom types (N/O/S–N/O/S → H-bond; C–C → hydrophobic; N/O/S/P → electrostatic; else van der Waals), then only kept if distance is within the limit for that type. All downstream metrics (BSA, contact area, complementarity) use only these filtered contacts.

**Libraries / functions used (interface_analyser)**

- **BioPython:** `PDBParser.get_structure('crystal', pdb_file)`; `structure.get_chains()` → list of chains; `chain.get_atoms()` for each chain. **Interface residues:** `NeighborSearch(atom_list).search_all(threshold, 'R')` → residue pairs within threshold, different chains. **Contact residues:** interface residues filtered by applying the same criteria (distance + type-specific limits) to each residue’s atoms vs the other chain (`_residue_passes_contact_criteria`). **BSA:** when `Bio.PDB.SASA.ShrakeRupley` is available: `Structure`, `Model` to build single-chain structures and a **two-chain-only** `Model` for the complex term; `ShrakeRupley(...).compute(entity, level='S')`; BSA = SASA(chain1 alone) + SASA(chain2 alone) − SASA(two-chain complex). Fallback: per-atom estimate (~15 Å² per contact atom). Per atom: `atom.coord`, `atom.element`, `atom.name`, `atom.parent`.
- **NumPy:** `np.linalg.norm()`, `np.mean()`, `np.min()`, `np.max()`, `np.std()` on distance lists.

### 2.1 Contacts

**Reader-friendly formula**  
A *contact* is a pair of atoms (one from each chain) whose distance is within the global cutoff (5.0 Å) and also within the limit for their contact type (e.g. H-bond ≤ 3.5 Å). **Interface residues** are identified by NeighborSearch: residue pairs from different chains that have at least one atom–atom pair within the distance threshold. **Contact residues** are those interface residues that pass the same criteria as the contact list: the residue must have at least one atom that, with some atom on the other chain, satisfies both the distance cutoff and the type-specific distance limit for that pair. Filtering is done by applying these criteria (not by membership in the contact list). The *contact count* is the number of atom pairs that pass the criteria.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher `contact_count` | More atom–atom proximity across the interface (often larger and/or tighter interfaces; can also reflect larger proteins and more atoms). |
| Lower `contact_count` | Fewer cross-chain atom contacts (smaller interface or weaker/looser contact). |
| Contact type counts dominated by one type | Indicates the interface chemistry is skewed (e.g. many C–C → hydrophobic; many charged → electrostatic-heavy). Counts are **atom-pair counts**, so one residue pair can contribute many contacts. |

**Contact count**

**Definition:** Number of atom pairs (chain1, chain2) with:

```
||x1 - x2|| <= d_contact
```

**and** within the type-specific distance limit for their contact type.

**Implementation (simplified):** For each atom pair (one from each chain), compute distance. If distance ≤ contact_distance (5 Å), classify the pair by atom types (N/O/S–N/O/S → H-bond; C–C → hydrophobic; N/O/S/P → electrostatic; else van der Waals). Keep the contact only if distance is within the limit for that type (e.g. H-bond ≤ 3.5 Å, electrostatic ≤ 5 Å). Contact count = number of kept pairs. Summary counts per type are reported.

- **Script:** `_find_contacts()` builds raw list; `_filter_contacts_by_limits()` classifies each contact and keeps only those with `_contact_within_limit(contact_type, distance)`.  
- **Output:** Summary counts per type (H-bonds, electrostatic, hydrophobic, van der Waals); no full contact list.

### 2.2 Buried surface area (BSA)

**Reader-friendly formula**  
*Buried surface area* is the surface area of each molecule that becomes hidden when the two chains form a binary complex. The script follows the **Biopython Interface Analysis** approach ([https://biopython.org/wiki/Interface_Analysis](https://biopython.org/wiki/Interface_Analysis)): when `Bio.PDB.SASA.ShrakeRupley` is available, BSA = SASA(chain1 alone) + SASA(chain2 alone) − SASA(two-chain complex). SASA is computed with the Shrake–Rupley rolling-ball algorithm (probe 1.4 Å). Isolated SASA uses deep-copied single-chain `Structure`/`Model` builds; **complex SASA uses a `Model` that contains only those two chains** (deep copies at their current coordinates). In a multi-chain assembly, each reported pair (e.g. A–B, A–C, A–D) therefore receives a **pairwise** BSA that does not subtract the SASA of the full lattice. If SASA is not available (or computation fails), BSA is estimated as the number of atoms in contact (within contact_distance) × 15 Å².

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher BSA | Larger interface footprint (more surface becomes buried upon complex formation). Often correlates with stronger/biologically relevant interfaces, but context matters. |
| Lower BSA | Smaller interface footprint (less buried surface). |
| BSA from fallback (15 Å²/contact atom) | A rough proxy only; absolute values are less reliable than when computed from Shrake–Rupley SASA. |

**Mathematical definition**

```
BSA = SASA(chain1) + SASA(chain2) - SASA( dimer: chain1 + chain2 only )
```

(when using ShrakeRupley); otherwise (fallback):

```
BSA ≈ |contact_atoms| * 15 Å^2
```

**Implementation (simplified):** When ShrakeRupley is available: compute SASA of chain1 alone, chain2 alone, and a **two-chain** model (deep copies of both chains only); BSA = SASA(chain1) + SASA(chain2) − SASA(that dimer). Otherwise: count atoms in contact (within contact_distance) and multiply by 15 Å² per atom.

- **Script:** If `Bio.PDB.SASA.ShrakeRupley` is available: `_compute_sasa_chain_isolated()` for each chain; build `Model(0)` with deep-copied chain1 and chain2, `ShrakeRupley().compute(..., level='S')` on that model for `sasa_complex`; BSA = sasa1 + sasa2 − sasa_complex. Else: `_estimate_buried_surface_area_fallback()` using atoms within `contact_distance` × 15 Å².  
- In `_calculate_buried_surface_area()` (and `_estimate_buried_surface_area_fallback()`).

### 2.3 Contact area

**Reader-friendly formula**  
*Contact area* is the total area attributed to the interface, computed per contact: if two atoms (with radii r₁, r₂) are closer than r₁ + r₂, they “overlap”; the script assigns an area to that overlap using a circular-cap formula: π × (smaller radius)² × (1 − d/(r₁+r₂)), then sums over all contacts.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher contact area | More/stronger geometric overlap between contacting atoms → often tighter/local packing at the interface. |
| Near-zero contact area with nonzero contacts | Contacts are within distance limits but mostly not “overlapping spheres” under this simplified model (e.g., distances near the cutoff); interpret alongside BSA and contact distances. |

**Definition:** Total area attributed to contact patches between the two chains.

**Mathematical definition:** For each contact with distance d and radii r_1, r_2: if d < r_1 + r_2 then

```
area = pi * min(r1, r2)^2 * (1 - d/(r1 + r2))     if d < (r1 + r2)
```

clamped to >= 0. Contact area = sum over contacts.

**Implementation (simplified):** For each contact, get distance d and element-based radii r₁, r₂ (default 1.7 Å). If d < r₁ + r₂, add area = π × (min(r₁,r₂))² × (1 − d/(r₁+r₂)), clamped to ≥ 0. Contact area = sum over all contacts.

- **Script:** In `_calculate_contact_area()`.

### 2.4 Average / min / max contact distance and std

**Reader-friendly formula**  
These are standard statistics over the list of contact distances (after type filtering): mean, min, max, and standard deviation.

**What it means (how to interpret)**

| Metric | Value pattern | Typical interpretation |
| --- | --- | --- |
| `average_contact_distance` | Lower | Tighter/closer contacts on average. |
| `average_contact_distance` | Higher | Looser interface contacts on average (still within type limits). |
| `contact_distance_std` | Lower | More uniform contact distances (more consistent packing). |
| `contact_distance_std` | Higher | More heterogeneous contact distances (mix of tight and loose contacts). |

**Implementation (simplified):** Take the list of contact distances (one per atom pair kept after type filtering). Average = mean of list; min = smallest; max = largest; std = standard deviation. Implemented with NumPy: `np.mean`, `np.min`, `np.max`, `np.std`.

- **average_contact_distance:** `np.mean(distances)` over all contact distances.  
- **min_contact_distance:** `np.min(distances)`.  
- **max_contact_distance:** `np.max(distances)`.  
- **contact_distance_std:** `np.std(distances)`.  
- In `_analyse_contacts()`.

### 2.5 Interface complementarity

There are two related complementarity scores, both in [0, 1].

#### 2.5.1 Shape complementarity

**Reader-friendly formula**  
*Shape complementarity* measures how “tight” the interface is in terms of contact distances. An ideal contact distance is taken as 3.5 Å. For each contact, the deviation from 3.5 Å is computed; the score is max(0, 1 − average_deviation/2), so smaller average deviation gives a higher score.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Score near 1 | Contact distances cluster tightly around 3.5 Å → interface looks “tightly packed” by this simple proxy. |
| Score near 0 | Contact distances are, on average, far from 3.5 Å → looser/less tight packing by this proxy. |

**Definition:** Score in [0, 1]; higher = tighter fit. Based on deviation from an “optimal” contact distance.

**Implementation (simplified):** Ideal distance = 3.5 Å. For each contact, 

deviation = |distance − 3.5|. 

**Mathematical definition:** For each contact distance d, deviation delta = |d - 3.5|; average deviation delta_bar. Then:

```
delta_i = abs(d_i - 3.5)
delta_bar = mean(delta_i)
complementarity_shape = max(0, 1 - delta_bar/2)
```

- In `_calculate_interface_complementarity()`; reported as `interface_complementarity`.

#### 2.5.2 Charge-tag complementarity (contact-based)

**Reader-friendly formula**  
*Charge-tag complementarity* is a contact-based diagnostic that measures how often charged contact pairs are *opposite-sign* vs *same-sign* among the filtered atom–atom contacts. The score is the fraction of charged–charged contacts that are opposite-sign; value in [0, 1]. This is distinct from **electrostatic complementarity (EC)** (McCoy et al. 1997), which correlates electrostatic potentials on facing surfaces (Section 2.5.4).

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Score near 1 | Most charged–charged contacts are opposite-sign (more favourable sign pairing among the charged contacts present). |
| Score near 0 | Mostly same-sign charged–charged contacts, or **no charged–charged contacts** at all (in which case the score is defined as 0). |
| Very small `charge_complementarity_total` | Score becomes noisy; interpret with the diagnostic counts (`opposite`, `same`, `total`). |

**Mathematical definition:** Let N_\text{opp} be the number of charged–charged contacts where the residues have opposite sign, and N_\text{same} the number where they have the same sign. Neutral residues are ignored. Then:

```
if (N_opp + N_same) == 0:
    complementarity_charge = 0
else:
    complementarity_charge = N_opp / (N_opp + N_same)
```

**Implementation (simplified):** For each contact, get the charge sign of the two residues: ASP/GLU = −1, ARG/LYS/HIS = +1, else 0. If either residue is neutral, skip the contact. If both are charged, count **opposite** when q1 × q2 < 0, **same** when q1 × q2 > 0. Score = opposite / (opposite + same), or 0 if there are no charged–charged contacts. Result is a fraction in [0, 1].

Residue charges (same set as polarity summary):  

- Negative: ASP, GLU → −1  
- Positive: ARG, LYS, HIS → +1  
- All other residues → 0 (ignored in this metric).

**Scope:** Counts are over *all* atom–atom contacts in the interface (each contact is one pair of atoms within distance; the same residue pair can contribute many contacts). So the score can take any value in [0, 1]. It will be exactly **0** when there are no charged–charged contacts or when all such contacts are same-sign (e.g. LYS–ARG). It will be exactly **1** when all charged–charged contacts are opposite-sign (e.g. only ASP–ARG, GLU–LYS). Uneven polarity counts (e.g. 21 vs 6 charged residues) do not force a fractional score: which residue pairs actually form contacts and their signs determine the fraction.

- Implemented in `_calculate_charge_complementarity()`; reported as `charge_complementarity`. Diagnostic counts: `charge_complementarity_opposite`, `charge_complementarity_same`, `charge_complementarity_total` (and in the text report when total > 0).

#### 2.5.3 Charge-tag complementarity density (per interface area)

**Reader-friendly formula**  
*Charge-tag complementarity density* is the number of opposite-sign charged–charged contacts per unit interface area (Å⁻²). It combines the count of opposite-sign charged contacts with interface size, so larger interfaces are comparable to smaller ones.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher density | More opposite-sign charged contacts per Å² of interface area. |
| Lower density | Fewer opposite-sign charged contacts per Å² (either fewer opposite contacts or a larger interface area). |
| `None` | No usable interface area denominator (both contact area and BSA are 0). |

Definition: (number of opposite-sign charged–charged contacts) / area
with area = contact_area if contact_area ≥ 0.01 Å², else BSA.
Unit: per Å² (e.g. 0.05 → 0.05 opposite-sign contacts per Å²).
Also stored: charge_complementarity_density_denominator = 'contact_area' or 'buried_surface_area' so you know which area was used.
If both contact area and BSA are 0, charge_complementarity_density is None.

It's a size-normalised measure: interfaces with the same fraction of opposite-sign contacts can be distinguished by how many such contacts there are per Å² (e.g. 2 vs 20 opposite contacts on a 100 Å² interface).

**Why contact area vs BSA:** The denominator is **contact area** when it is ≥ 0.01 Å² (same contact-based definition as the contact list and contact area). When contact area is zero or negligible, **buried surface area (BSA)** is used so the metric remains defined. Contact area is preferred because it is derived from the same atom–atom contacts that define the charged-contact counts; BSA is a different measure (solvent-accessible surface buried) but is a sensible fallback when there is no overlap-based contact area.

**Mathematical definition**

```
area = contact_area   if contact_area >= 0.01 else BSA
charge_complementarity_density = N_opp / area
```

(undefined / not reported when both contact_area and BSA are 0.)

**Implementation (simplified):** After computing charge-tag complementarity counts, set area = contact_area if contact_area ≥ 0.01, else buried_surface_area. If area > 0, charge_complementarity_density = opposite / area (units: per Å²). Also store charge_complementarity_density_denominator = `'contact_area'` or `'buried_surface_area'` to record which was used. If area is 0, report `charge_complementarity_density` = None.

- Reported as `charge_complementarity_density` (per Å²) and `charge_complementarity_density_denominator`; in the text report when density is not None.

#### 2.5.4 Electrostatic complementarity (EC; McCoy, Epa & Colman, 1997)

**Definition.** Electrostatic complementarity (EC) is reported using the formulation of McCoy, Epa & Colman (1997): EC is the Pearson correlation \(r\) between electrostatic potential values on *facing* points across the interface, evaluated in the partner field sense. In FoldKit, EC is a correlation score \(r \in [-1, 1]\); positive values indicate complementary (opposite-sign) potentials across the interface.

**Important distinction.** EC is not derived from the atom–atom contact list. It is computed from electrostatic potentials on sampled solvent-accessible surface points classified as interfacial, then paired across the interface by nearest-neighbour mapping. Do not conflate EC with the older **charge-tag** contact metrics (`charge_complementarity*`), which are based only on charged residue identities in the filtered contact set (Sections 2.5.2–2.5.3).

**Reported fields (per interface; `InterfaceAnalyserEC`):**

- **`ec_r`**: EC correlation \(r\) for the interface.
- **`ec_n_pairs`**: number of facing point pairs used in the correlation.
- **`ec_density`**: \(r/\mathrm{BSA}\) (Å\(^{-2}\)); denominator label stored as `ec_density_denominator` (currently `buried_surface_area`).

**References and full notation.** See `metrics/EC_details.md` (definitions, Fisher z aggregation for lattice summaries, density normalisation, and references).

### 2.6 Contact type counts (incl. H-bonds)

**Implementation (simplified):** When filtering contacts by type and distance, each kept contact is classified as one of: hydrogen_bond (N/O/S–N/O/S, ≤3.5 Å), electrostatic (N/O/S/P, ≤5 Å), hydrophobic (C–C, 3.5–4.5 Å), van_der_waals (else, ≤5 Å). Count how many contacts fall into each class. `hydrogen_bonds` = count for hydrogen_bond.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| High H-bond count | Many close N/O/S–N/O/S contacts within 3.5 Å; suggests polar, directional contacts (note: this is a simplified atom-type rule, not a full donor/acceptor geometry check). |
| High hydrophobic count | Many C–C contacts in 3.5–4.5 Å range; suggests a hydrophobic packing-dominated interface. |
| High electrostatic count | Many contacts involving N/O/S/P within 5 Å; suggests a polar/charged-rich interface (simplified classification). |
| High van der Waals count | Broad “other” category; can indicate many mixed/weak contacts not captured by the other definitions. |

- **Script:** `contact_type_counts` from `_filter_contacts_by_limits()`; `hydrogen_bonds` = `contact_type_counts['hydrogen_bond']`.

### 2.7 Polarity and accessibility of interface residues

**Reader-friendly formula**  
For each interface (chain pair), contact residues on each chain are classified into **charged**, **polar**, **apolar**, and **other** based on residue type. This gives both counts and fractions per category. Accessibility is quantified via per-residue SASA (Å²) computed with Shrake–Rupley: for contact residues on each chain, the script reports the **average SASA** and the **fraction of residues with SASA > 0** (accessible contact residues).

**Implementation (simplified):** For each contact residue on a chain, classify by residue name: charged (ARG, LYS, ASP, GLU, HIS), polar (SER, THR, ASN, GLN, TYR, CYS), apolar (ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, GLY), other (anything else). Count how many fall in each class per chain → polarity_counts. Fractions = count / total contact residues on that chain. For accessibility: compute per-residue SASA with ShrakeRupley; report average SASA and the fraction of contact residues with SASA > 0 per chain.

**Residue polarity classes**

- **charged:** ARG, LYS, ASP, GLU, HIS  
- **polar (uncharged):** SER, THR, ASN, GLN, TYR, CYS  
- **apolar:** ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, GLY  
- **other:** any other residue/ligand name

For the set of contact residues on a given chain:

- **polarity_countschainX:** counts in each class (charged, polar, apolar, other).  
- **polarityfractionschainX:** per-class fraction:

```
f_class = count_class / |R_c|
```

when `|R_c| > 0` (`R_c` = set of contact residues on that chain).

**Accessibility**

Per-residue SASA is computed with Shrake–Rupley (`ShrakeRupley(probe_radius=1.4, n_points=100).compute(model, level='R')`). For contact residues on each chain:

- **averagesasachainX:** mean SASA over contact residues with SASA information.  
- **accessiblefractionchainX:** fraction of those residues with SASA > 0.

Formally, for SASA values s_i over contact residues on a chain with available SASA:

```
average_sasa = (1/n) * sum_i(s_i)
accessible_fraction = count_i( s_i > 0 ) / n
```

(0 if no SASA values are available).

**What it means (how to interpret)**

| Metric | Value pattern | Typical interpretation |
| --- | --- | --- |
| `average_sasa` (contact residues) | Low | Contact residues are generally buried / less solvent exposed (often expected in a tight interface). |
| `average_sasa` (contact residues) | High | Contact residues are still quite exposed (could indicate a shallow/edge interface, or that “contact residues” include residues that barely touch). |
| `accessible_fraction` | Near 0 | Most contact residues have SASA = 0 → fully buried (given the SASA computation). |
| `accessible_fraction` | Near 1 | Most contact residues have SASA > 0 → still solvent-accessible. |
| `residues_with_sasa` | Small / 0 | SASA was unavailable or failed for many residues; treat accessibility metrics as missing/placeholder. |

**Important nuance**

- The per-residue SASA used here is **SASA in the complex/model** for those residues, not per-residue **ΔSASA** upon binding. The reported quantities describe solvent exposure in the supplied coordinates, not the change in accessible surface area associated with a particular association pathway.

**Script calculation**

- **Polarity:** `_classify_residue_polarity()` maps residue names into classes; `_summarize_polarity_and_accessibility()` accumulates `polarity_chain1`, `polarity_chain2` (counts) and `polarity_fractions_chain1`, `polarity_fractions_chain2` (fractions).  
- **Accessibility:** `_compute_residue_sasa_for_chains()` runs Shrake–Rupley with `level='R'` on a deepcopy of the model and records per-residue SASA for the two chains; `_summarize_polarity_and_accessibility()` then computes `average_sasa` and `accessible_fraction` for contact residues per chain, returned as `accessibility_chain1` and `accessibility_chain2`.

If `Bio.PDB.SASA.ShrakeRupley` is not available or SASA computation fails, accessibility fields default to 0 and only polarity counts/fractions are reported.

### 2.7.1 Lattice reference metrics (multi-copy assemblies)

**When:** Optional, from `InterfaceAnalyser.analyse_interfaces(..., reference_chain_id='A')` or `interface_analyser_lattice_charge.py --reference-chain A` / `interface_analyser_lattice_ec.py --reference-chain A` (and via `crystal_packing_analyser.py --reference-chain A`). Intended for a **pre-built multi-chain structure** (many chains, same molecule repeated), not for inferring the crystal lattice from the ASU alone.

**Purpose.** One chain ID is designated the **reference copy**. Its solvent-accessible surface is evaluated (i) in isolation and (ii) with all other supplied chains present in the same coordinate frame. A further summary gives the fraction of reference residues that have at least one atom within `contact_distance` of a **different** chain.

**Fields (in `summary`)**

- **`sasa_reference_isolated`:** Shrake–Rupley SASA (Å²) of the reference chain in isolation (same probe/settings as BSA: 1.4 Å, `level='S'` on a single-chain structure).
- **`sasa_reference_in_cluster`:** Sum of per-residue SASA (Å²) for that chain after `ShrakeRupley.compute(..., level='R')` on a **deepcopy of the full model** (all chains stay in place, so neighbouring copies occlude the probe).
- **`lattice_burial_fraction`:** When both SASA values are available and isolated SASA > 0:

```
lattice_burial_fraction = clamp( 1 - sasa_reference_in_cluster / sasa_reference_isolated, 0, 1 )
```

  Larger values indicate greater occlusion of the reference copy by neighbours in the static model (lower solvent exposure than in isolation). The quantity is an **empirical packing index**, not a thermodynamic **ΔSASA** for a specified binding or assembly reaction.

- **`lattice_contact_residue_fraction`:** Fraction of residues on the reference chain that have ≥1 atom within `contact_distance` of an atom on another chain (any chain ID ≠ reference), using `NeighborSearch` on the full model.

- **`lattice_charge_complementarity`:** **Charge-tag** complementarity between the reference chain and all other chains in the model, computed from filtered atom–atom contacts (Sections 2.4–2.6). Charged residues are classified as ASP/GLU (−) and ARG/LYS/HIS (+). The score is:

```
lattice_charge_complementarity = N_opposite / (N_opposite + N_same)
```

  where counts are over **charged–charged** atom-contact pairs (neutral residues are ignored). If no charged–charged contacts exist, the score is reported as 0.0 and totals are 0.

  Diagnostic counts are also stored:
  - `lattice_charge_complementarity_opposite`
  - `lattice_charge_complementarity_same`
  - `lattice_charge_complementarity_total`

- **`lattice_charge_complementarity_density`:** Opposite-sign charged–charged contacts per unit **reference buried area** in the lattice context:

```
lattice_charge_complementarity_density = N_opposite / (sasa_reference_isolated - sasa_reference_in_cluster)
```

  Units are opposite-sign contacts per Å². If the reference buried area cannot be computed or is ≤ 0, density is reported as `None`. The denominator label is stored as `lattice_charge_complementarity_density_denominator = 'reference_buried_area'`.

- **`lattice_charge_metrics`:** Extended lattice-wide charge metrics (a structured dict). This block is present when `reference_chain_id` is supplied. It contains four complementary views:

  - **Atom-contact based (overall)**: `lattice_charge_metrics['atom']` reproduces the scalar `lattice_charge_complementarity` definition, with the diagnostic counts.

  - **Residue-pair based (overall)**: `lattice_charge_metrics['residue_pair']` counts **unique charged residue pairs** between the reference chain and any neighbour chain. Each distinct pair contributes at most 1 to the counts, even if many atom–atom contacts exist between the same residues. This reduces “atom-contact multiplicity” bias and is often closer to “how many distinct salt-bridge-like pairings exist”.

  - **Distance-weighted residue-pair complementarity (overall)**: `lattice_charge_metrics['residue_pair_weighted']` uses the same unique residue pairs, but weights each pair by \(1/d^2\), where \(d\) is the **minimum** atom–atom distance observed between the two residues among the filtered contacts. This emphasises tighter charged approaches.

  - **Partner-resolved metrics**: `lattice_charge_metrics['by_partner_chain']` is a dict keyed by neighbour chain ID (e.g. `'B'`, `'C'`). For each partner chain, it reports residue-pair and weighted residue-pair complementarity for contacts between the reference chain and that partner only. This shows which crystal contacts dominate the lattice electrostatics.

- **McCoy EC lattice metrics (if present):** When the electrostatic complementarity (EC) analyser is used (`interface_analyser_lattice_ec.py`), lattice-wide EC summaries are reported instead of charge-tag lattice metrics:
  - `lattice_ec_r_bsa_weighted`: lattice EC \(r\) aggregated by BSA-weighted Fisher z.
  - `lattice_ec_r_npairs_weighted`: lattice EC \(r\) aggregated by n\_pairs-weighted Fisher z.
  - `lattice_ec_density_bsa_weighted`, `lattice_ec_density_npairs_weighted`: lattice EC density (Å\(^{-2}\)) normalised by the reference buried area.
  - `lattice_ec_by_partner_chain`: per-partner records with \(r_k\) and \(n_{\mathrm{pairs},k}\).

  Definitions and aggregation details are in `metrics/EC_details.md`.

**What it means (how to interpret)**

| Pattern | Typical interpretation |
| --- | --- |
| Higher `lattice_burial_fraction` | The focal chain’s surface is more occluded by other copies in this assembly (tighter static packing around that copy). |
| Lower `lattice_burial_fraction` | The focal chain remains relatively exposed despite neighbours (sparse placement or few contacts). |
| Higher `lattice_contact_residue_fraction` | More of the focal chain’s residues are near another chain (contact-based, not SASA-based). |
| Higher `lattice_charge_complementarity` | Among charged–charged lattice contacts involving the focal chain, a higher fraction are opposite-sign (more electrostatic “matching” in the contact set). |
| Very small `lattice_charge_complementarity_total` | Score becomes noisy; interpret with the diagnostic counts (`opposite`, `same`, `total`). |
| `lattice_metrics_error` | Reference chain ID absent from the parsed structure (verify chain identifiers in the input file). |

**Implementation (simplified):** `_compute_lattice_reference_metrics()` → `_compute_sasa_chain_isolated()`, `_sum_residue_sasa_for_chain_in_full_model()`, `_lattice_contact_residue_fraction()`, `_lattice_charge_complementarity()`.

For the crystallographic workflow and interpretation of these quantities, see Section 2.7.2.

### 2.7.2 Multi-copy coordinate models for crystallographic interface analysis

**Scope.** In macromolecular crystallography, contacts that define the crystal packing are often discussed in terms of **crystal contacts** between symmetry-related copies of the asymmetric unit (ASU). The present implementation does **not** apply space-group operators or rebuild the lattice from a single ASU. It operates on **Cartesian coordinates** exactly as supplied. A **multi-copy model** is therefore a single structure file (PDB/mmCIF) in which **multiple chains** (typically distinct chain IDs) place several copies of one or more molecules in the same frame—e.g. an oligomer, a user-generated symmetry expansion, or an artificial “supercell” fragment exported from a graphics or modelling program. All metrics are **geometric summaries** of that static ensemble.

**Standard pairwise output.** For every pair of distinct chains that share a model and satisfy the contact criteria, the analyser reports the same quantities as for a binary complex (Sections 2.1–2.7): filtered atom–atom contacts, **pairwise** BSA (Section 2.2: isolated SASA for each chain minus SASA of a **two-chain-only** complex), contact area, complementarity measures, polarity and accessibility of contact residues, and optional interface RMSD where applicable. In a multi-copy file, **each distinct chain pair** (e.g. A–B, A–C, A–D) is evaluated independently with that definition. **Summed** BSA over all reported pairs remains **not** equivalent to a single “total burial” of one molecule in the full lattice, because overlapping interface patches are counted in more than one pair. Treat the sum as an **aggregate over pairwise interfaces**, not a global lattice thermodynamic quantity.

**Isolated SASA block (summary).** For every chain ID that appears in any reported interface, the structure-level summary can list **isolated** Shrake–Rupley SASA (each chain alone, probe radius 1.4 Å) and their sum. That sum is a **bookkeeping** aggregate over those chains; it does not assign a single physical “monomer in crystal” state unless the modelling context supports that interpretation.

**Chain filtering (`focus_chains` / `--chains`).** Restricting analysis to a set of chain IDs limits which **pairs** are evaluated: a pair is retained if **at least one** of its chains lies in the filter. This is used to centre reporting on one molecule (e.g. all interfaces involving a chosen focal chain ID) without discarding the partner chains required for pairwise BSA and contact definitions.

**Focal-copy metrics (`reference_chain_id` / `--reference-chain`).** Section 2.7.1 defines the reported fields. Conceptually, they compare one designated chain (the **reference copy**) in two geometric settings:

1. **Isolated reference SASA** (`sasa_reference_isolated`): the solvent-accessible surface area of the reference chain alone—i.e. the chain embedded in a structure that contains **no other macromolecular chains**, with the same conformation. This matches the **isolated-chain** arm of the standard BSA construction (Section 2.2): a rolling-sphere (probe) exposure measure in Å².

2. **Embedded reference SASA** (`sasa_reference_in_cluster`): the **sum of per-residue** SASA values for the reference chain computed after Shrake–Rupley at residue level on a **deep copy of the full multi-chain model**. All other chains remain as **hard obstacles** to the probe: wherever another copy packs against the reference chain, residue-level SASA on the reference chain typically decreases.

3. **Lattice burial fraction** (`lattice_burial_fraction`): when both SASAs are defined and the isolated value is positive,

```
lattice_burial_fraction = clamp( 1 - sasa_reference_in_cluster / sasa_reference_isolated, 0, 1 )
```

   This **dimensionless occlusion index** approximates the **fraction of the reference chain’s isolation SASA that is no longer solvent-accessible** when neighbour chains occupy the same coordinate set. It is **not** a thermodynamic **ΔSASA** for a defined reaction, nor an experimental or computational “surface burial energy.” Isolated SASA uses **structure-level** output (`level='S'`) on a one-chain construct; embedded SASA is the **sum of per-residue** terms (`level='R'`) on the full assembly. The ratio serves as a **comparative proxy** across packings or models; it need not agree numerically with a single global SASA for the reference chain computed in the full model.

4. **Contact residue fraction** (`lattice_contact_residue_fraction`): the proportion of residues on the reference chain that have **at least one atom** within `contact_distance` of **any atom** belonging to a **different** chain ID. This is a **pure distance-threshold** statistic: it does not use SASA and does not distinguish biological interfaces from crystal contacts by symmetry class. Any hetero-chain partner counts, including solvent or ligand if they appear as separate chains.

**Physical meaning.** The SASA-based quantities use the same geometric construction as conventional rolling-probe (Shrake–Rupley) protein surfaces: probe radius 1.4 Å and a discrete solvent-excluded surface approximation. **Isolated** SASA is the reference chain’s accessible area with no other macromolecular chains present. **Embedded** SASA is the sum of per-residue accessible areas for that chain when every supplied chain is present, so neighbouring molecules exclude the probe from parts of the reference surface. **Lattice burial fraction** reduces that contrast to one scalar for comparing models. **Contact residue fraction** is independent of SASA: it is the proportion of reference residues with any hetero-chain atom within the distance cut-off, and may diverge from SASA-based occlusion (for example when contacts are long-ranged or geometrically shallow).

**Mathematical meaning.** Let `S_iso` = `sasa_reference_isolated` and `S_emb` = `sasa_reference_in_cluster`. The burial fraction is a **clipped linear contrast** `1 - S_emb / S_iso` on `[0, 1]`. For fixed `S_iso`, smaller embedded SASA implies a larger burial fraction. The index is **not** additive across chains and **not** an enthalpy or free energy without a calibrated energy model. The contact fraction is a **ratio of finite counts** and depends explicitly on the chosen cutoff and on which entities are assigned distinct chain identifiers in the file.

**Biological and crystallographic interpretation.** In a multi-copy model that reproduces or approximates crystal or predicted packing contacts, pairwise BSA and contact-type summaries permit **qualitative comparison** of interface extent and chemistry between coordinate sets, by the same logic as for binary complexes. Focal-copy metrics permit comparison of steric occlusion of the reference chain across models, and comparison of the fraction of reference residues within the contact cut-off of another chain. They do **not** determine space group, validate symmetry, or measure lattice stability; those require diffraction evidence, symmetry treatment appropriate to the experiment, and thermodynamic or simulation methods where relevant. In `crystal_packing_analyser.py`, `reference_chain_id` is passed to the interface stage so the focal-copy fields appear under `interface_analysis.summary` in each per-structure JSON file.

**Limitations.** (i) Interfaces are not classified as biological or crystallographic; that distinction uses external crystallographic criteria. (ii) **Static coordinates:** no conformational entropy, flexibility, or explicit solvent structure. (iii) **Input quality:** missing atoms, severe clashes, or unrealistic expansions bias SASA and contact statistics. (iv) **Chain identifiers:** duplicate or non-standard labelling changes which pairs are evaluated and which atoms are treated as hetero-chain. (v) **Summed pairwise BSA** over many chain pairs is not equivalent to a single global burial metric for one molecule.

### 2.8 Interface RMSD (contact residues, CA atoms)

**Reader-friendly formula**  
For each interface between two chains, a simple **interface RMSD** is computed by superposing the CA atoms of **matching contact residues** on the two chains and reporting the root-mean-square deviation (RMSD, Å). Matching is based on residue name and residue number (and insertion code) among residues that are in the contact-residue sets on both chains. This is most meaningful for homodimeric interfaces where the two chains are sequence-identical and aligned.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Lower RMSD | The two sides of the interface (matched contact residues) superpose well → more symmetric / similar interface geometry (often expected for homodimers). |
| Higher RMSD | Less geometric similarity between the matched contact residues (may reflect heterodimers, conformational differences, or mismatched residue mapping). |
| N/A / missing | Not computable (e.g. too few matched residues, missing CA atoms, or superposition not available). |

**Mathematical definition**

Let `M` be the set of matched contact-residue pairs between chain 1 and chain 2. For each pair `i` in `M`, take the Cα coordinates `x_i` (chain 1) and `y_i` (chain 2). After an optimal rigid-body superposition (rotation `R`, translation `t`) that minimises RMSD:

```
RMSD = sqrt( (1/|M|) * sum over i in M of ||R*y_i + t - x_i||^2 )
```

(vector norm in Å; `|M|` = number of matched pairs.)

**Implementation (simplified):** Match contact residues between the two chains by (resname, residue number, insertion code). For each matched pair, take the CA atom coordinates. If there are at least 3 matched CA pairs, use BioPython Superimposer to find the rotation and translation that minimise RMSD; RMSD = sqrt(mean of squared distances between superposed CA pairs). If fewer than 3 pairs or superposition unavailable, report N/A (or 0 in the summary).

**Script calculation**

- Collect `contact_residues_chain1` / `contact_residues_chain2` (as residue objects).  
- Build matches by residue `(resname, residue_number, insertion_code)` across the two contact-residue sets.  
- For each matched pair, extract CA atoms; skip pairs lacking CA.  
- If there are at least 3 CA pairs:
  - Use Biopython `Superimposer` (`Bio.PDB.Superimposer`) to compute the optimal superposition and RMSD.  
  - Store per-interface value as `interface_rmsd_ca`.
- In the structure-level summary, **average interface RMSD** is:

```
average_interface_rmsd_ca = (1/N_rmsd) * sum_over_valid_interfaces( RMSD_k )
```

where `N_rmsd` is the number of interfaces with RMSD > 0. This value appears in per-structure output as `average_interface_rmsd_ca`.

If there are no valid matches or the superposition fails, `interface_rmsd_ca` and the summary’s `average_interface_rmsd_ca` default to 0.

### 2.9 Summary (per structure)

- **total_buried_surface_area:** Sum of BSA over all chain-pair interfaces.  
- **total_contact_area:** Sum of contact_area over interfaces.  
- **average_buried_area_per_interface:** `total_buried_surface_area / interface_count`.  
- **average_contact_area_per_interface:** `total_contact_area / interface_count`.  
- **average_interface_rmsd_ca:** Mean CA-atom RMSD over interfaces with a valid RMSD estimate.  
- **total_interfaces:** Number of chain pairs with at least one contact.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher totals / averages | On average, the structure has larger/more contact-rich interfaces (by the definitions above). |
| Larger `total_interfaces` | More chain pairs are in contact (within the model and thresholds). |
| Nonzero `average_interface_rmsd_ca` | RMSD was computable for at least one interface (most meaningful for homomers). |

### 2.10 Post-processing: interface_molecule_report_csv.py

The interface analysers can write merged text reports (conventionally **`results.txt`** via `-o`). **`interface_molecule_report_csv.py`** parses those reports and emits CSV rows for each interface block, with the same numeric summaries as in Section 2 (contact counts, BSA, complementarity, polarity, accessibility, etc.). It does **not** re-read PDB files or recompute metrics. Optional filters: structure basename (as printed in the report, e.g. **`model_01.pdb`**), chain IDs (**`A`**, **`B`**), **`--group-by-chain`**, and merge options **`--combine-regex`** / **`--combine-glob`** (see **README.md**). For **atom-level** ASU contact rows (chain/residue/atom pairs), use **`contact_molecule_report_csv.py`** on **`contact_analyser.py`** output (Section 3.5).

**What it means (how to interpret)**

| Behaviour | Practical takeaway |
| --- | --- |
| Parser only (no recomputation) | The CSV exactly reflects the analyser’s text output; regenerate `results.txt` if you change any analysis settings. |
| One row per interface block | Good for plotting/filtering per-interface summaries without re-running structure parsing. |

**Derivation / references**

- Distance cutoffs and contact-type classifications are inspired by typical protein interface analyses, e.g.:
  - P. Chakrabarti & J. Janin, “Dissecting protein–protein recognition sites”, *Proteins* **47**, 334–343 (2002). DOI: `10.1002/prot.10085`.
  - J. Janin & B. Séraphin, “Genome-wide studies of protein–protein interaction”, *Curr. Opin. Struct. Biol.* **13**, 383–388 (2003). DOI: `10.1016/S0959-440X(03)00064-2`.
- **Interface definition and BSA** follow the Biopython Interface Analysis wiki ([https://biopython.org/wiki/Interface_Analysis](https://biopython.org/wiki/Interface_Analysis)): interface residues from NeighborSearch residue pairs within threshold; BSA = SASA(chain1) + SASA(chain2) − SASA(complex). The implementation uses `Bio.PDB.SASA.ShrakeRupley` (Shrake & Rupley, *J. Mol. Biol.* **79**, 351–371, 1973). DOI: `10.1016/0022-2836(73)90011-9` when available; otherwise a 15 Å² per contact-atom fallback.
- The simple complementarity score is a normalised measure of how close contact distances are to an “ideal” 3.5 Å, analogous to empirical interface “tightness” measures used in docking and packing quality checks.

---

## 3. contact_analyser.py

**Function selection / usage**

- Contact metrics are computed in `contact_analyser.py`, which can be:
  - Run directly on a PDB/CIF to analyse ASU contacts and crude crystal contacts, or
  - Invoked by `crystal_packing_analyser.py` for each structure.
- The main API entry is **`ContactAnalyser.analyse_contacts(pdb_file)`** (CLI uses **`_run_analysis()`** over input paths). It calls:
  - `_find_chain_contacts()` → ASU contacts (Section 3.1).
  - `_calculate_contact_density()` → Section 3.2.
  - `_estimate_crystal_contacts()` and `_calculate_cell_volume()` → Sections 3.3–3.4.
- **Contact distance:** Default 4.5 Å (candidate pairs); only contacts within the same **type-specific distance limits** as interface_analyser (H-bonds ≤3.5 Å, electrostatic ≤5 Å, hydrophobic 3.5–4.5 Å, van der Waals ≤5 Å) are kept.

**Libraries / functions used (contact_analyser)**

- **BioPython:** `PDBParser.get_structure()`, `structure.get_chains()`, `chain.get_atoms()`; `atom.coord`, `atom.element`, `atom.name`, `atom.parent` (residue). Distance: `np.linalg.norm(atom.coord - other_atom.coord)`.
- **NumPy:** `np.linalg.norm`, `np.mean`, `np.min`, `np.max` on distance lists; sum over chain atom counts for contact density.

### 3.1 ASU contacts

**Reader-friendly formula**  
*ASU contacts* are atom pairs from different chains within the ASU that are within the global distance (4.5 Å) and within the type-specific limit for their classified interaction type. Each contact stores chains, residues, atoms, distance, and type. Totals and averages are taken over this filtered list.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher ASU `total_contacts` | More inter-chain proximity within the asymmetric unit (could indicate strong crystal contacts or biological assembly depending on context). |
| Lower ASU `total_contacts` | Fewer inter-chain contacts in the ASU. |
| Type distribution skew | Interface/contact chemistry is dominated by that interaction class (with the same caveat: these are atom-pair counts). |

**Implementation (simplified):** For each atom pair from different chains, compute distance. If distance ≤ 4.5 Å, classify type (H-bond, electrostatic, hydrophobic, van der Waals) from atom types; keep only if distance is within the limit for that type. Each kept pair is one contact. total_contacts = len(contact list); average/min/max distance and counts per type from this list.

- **Contacts:** Pairs of atoms from different chains with distance ≤ `contact_distance` and within the type-specific limit for their classified type; each contact stores chain ids, residues, atoms, distance, and contact_type.  
- **Script:** In `_find_chain_contacts()`, after classifying with `_classify_contact_type()`, only append if `_contact_within_limit(contact_type, distance)`.  
- **total_contacts:** `len(asu_contact_list)` (only contacts within limits).  
- **average_distance, min_distance, max_distance:** From the filtered contact distances.  
- **contact_type_distribution:** Counts per `contact_type`.

**CLI text listing:** For each structure, the analyser can print a **Contact details** block: one line per kept contact with `chain1`, `chain2`, residue labels, atom names, distance (Å), and `contact_type`. If there are more than 200 contacts, the full list is written to a sidecar file **`{output_stem}_{structure_stem}_asu_contacts.txt`** (when `-o` is set, `{output_stem}` is derived from that path; otherwise the structure stem only). The sidecar header line starts with `# All N ASU contacts:`.

### 3.2 Contact density

**Reader-friendly formula**  
*Contact density* = (number of contacts) / (total number of atoms in all chains). It measures how many contacts there are per atom in the structure.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher contact density | Many contacts relative to atom count → generally more inter-chain connectivity (within the ASU contact definition). |
| Lower contact density | Fewer contacts per atom → less inter-chain connectivity by this definition. |

**Mathematical definition:**
```
contact_density = (number_of_contacts) / (total_number_of_atoms_in_all_chains)
```

**Implementation (simplified):** total_atoms = sum of atom counts over all chains. Contact_density = len(contacts) / total_atoms. In `_calculate_contact_density()`.

### 3.3 Crystal contact estimates (simplified)

**Reader-friendly formula**  
*Surface atoms* are those with fewer than 12 neighbours within 5.0 Å (low coordination). The number of *potential crystal contacts* is set equal to the number of surface atoms. *Estimated contact area* is surface_atoms × 4 Å². *Packing efficiency* is (molecular volume) / (cell volume), with molecular volume ≈ atom_count × 20 Å³.

**What it means (how to interpret)**

| Metric | Value pattern | Typical interpretation |
| --- | --- | --- |
| `surface_atoms` (heuristic) | Higher | More atoms are classified as “surface-like” by the neighbour-count rule (often larger proteins and/or more exposed structure). |
| `potential_crystal_contacts` | Higher | More potential sites for crystal contacts under this simplified proxy. |
| `estimated_contact_area` | Higher | Larger estimated contact area proxy (linear in `surface_atoms`). |
| `packing_efficiency` (proxy) | Higher | Larger fraction of the unit cell filled by the crude molecular volume estimate; interpret qualitatively only. |

**Implementation (simplified):** For each atom, count how many other atoms are within 5.0 Å. If count < 12, it is a surface atom. surface_atoms = number of such atoms. potential_crystal_contacts = surface_atoms. estimated_contact_area = surface_atoms × 4.0 (Å²). molecular_volume = atom_count × 20.0 (Å³); packing_efficiency = molecular_volume / cell_volume.

- **surface_atoms:** Atoms with < 12 neighbours within 5.0 Å (neighbour count over all atoms).  
- **potential_crystal_contacts:** Set equal to number of surface atoms.  
- **estimated_contact_area:** `surface_atoms * 4.0` (Å²).  
- **packing_efficiency:**  
  - Cell volume from unit cell (same formula as in packing_metrics).  
  - Molecular volume ≈ `atom_count * 20.0` (Å³).  
  - `packing_efficiency = molecular_volume / cell_volume`.

### 3.4 Unit cell volume (contact_analyser)

- Same parallelepiped formula as in Section 1.1; implemented in `_calculate_cell_volume()`.

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Same as Section 1.1 | This is the same unit-cell volume concept in the contact-analysis pipeline; interpret identically. |

### 3.5 Post-processing: contact_molecule_report_csv.py

The CLI of `contact_analyser.py` can write a merged text report (conventionally **`contact_results.txt`** via `-o`, distinct from interface reports (often **`results.txt`** from `interface_analyser_asu_charge.py` / `interface_analyser_asu_ec.py` and the lattice variants) or `packing_metrics.py`). **`contact_molecule_report_csv.py`** parses that report (and optionally a standalone **`contact_results_model_01_asu_contacts.txt`**-style sidecar) into CSV with one row per atom–atom contact: **`chain1`**, **`chain2`**, **`res1`**, **`atom1`**, **`res2`**, **`atom2`**, **`distance_A`**, **`contact_type`**, plus context columns **`set_label`** and **`structure_basename`**. It does **not** re-read PDB files or recompute distances. Filters and merge options mirror **`interface_molecule_report_csv.py`** (structure basename patterns, **`-m` / `--chains`**). Sidecar-only files lack progress lines; use **`--structure-basename model_01.pdb`** or the filename heuristic described in **README.md**. Multi-set output uses the same **`Set 'label' (patterns: …)`** headers as other analysers.

**What it means (how to interpret)**

| Behaviour | Practical takeaway |
| --- | --- |
| One row per atom–atom contact | Great for contact-level histograms; remember one residue pair can generate many rows. |
| Parser only (no recomputation) | If thresholds/definitions change, regenerate `contact_results.txt` first. |

**Derivation / references**

- The ASU contact criteria mirror the interface criteria in Section 2 with a slightly smaller candidate distance (4.5 Å) to focus on proximate contacts; this is similar in spirit to contact definitions in crystallographic tools such as PISA.
- The **crystal contact estimates** are heuristic:
  - “Surface atoms” are identified via reduced neighbour counts (cf. solvent-exposed residues defined by low coordination in packing analyses).
  - Per-atom area (~~4 Å²) and per-atom volume (~~20 Å³) reflect average SASA and volume for protein atoms in crystal structures.

---

## 4. crystal_packing_analyser.py

This script orchestrates `packing_metrics.py`, the **interface analysers** (`interface_analyser_base.py` via `interface_analyser_asu_charge.py` / `interface_analyser_asu_ec.py` and lattice variants), and `contact_analyser.py` on each structure. With `--compare`, it writes `batch_analysis_results.json` combining all per-structure outputs. It does not define new formulas; all metrics are as in Sections 1–3.

**Function selection / usage**

- When you run `crystal_packing_analyser.py` from the command line, you choose which input structures (files or directories) to analyse and optional `--compare` for a single combined JSON file.
- **Interface mode selection:** use `--interface-metrics charge` (charge-tag metrics) or `--interface-metrics ec` (McCoy electrostatic complementarity). Use `--reference-chain` when you want lattice reference metrics (Section 2.7.1) included in each structure’s `interface_analysis.summary`.
- Internally, it iterates over structures and calls the per-structure modules from Sections 1–3.

**What it means (how to interpret)**

| Behaviour | Practical takeaway |
| --- | --- |
| Orchestrator only | Use it to run the full pipeline; interpret numbers using the metric definitions in Sections 1–3. |
| Combined JSON (`--compare`) | Best starting point for plotting/summary statistics across a batch. |

---

## 5. structure_phylogeny.py

This script builds structure-based phylogenies from pairwise structural comparisons provided as an RMSD matrix (precomputed by your Coot LSQ/SSM pipeline, or optionally computed via `TM-align`).

It supports two related “distance matrix” modes for tree inference:

1. Default: use RMSD directly as the distance.
2. DALI-like “per-query pseudo Z-score” mode: convert RMSD to similarity, then to per-query z-scores, then to distances.

### 5.1 Default mode: RMSD distance matrix

**Mathematical definition**

Let `D` be the RMSD matrix with entries `D_ij >= 0`. The distance matrix for phylogeny construction is:

`d_ij = D_ij` for `i != j`, and `d_ii = 0`.

**Math (rendered)**

```
d_ij = D_ij   if i != j
d_ii = 0
```

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Smaller d_ij | Structures i and j are more similar (lower RMSD). |
| Larger d_ij | Structures are more dissimilar (higher RMSD). |

**Tree construction**

The script feeds this distance matrix to neighbour-joining (NJ) when available (otherwise it uses an internal UPGMA fallback).

---

### 5.2 DALI-like per-query pseudo Z-scores (RMSD -> distance)

This mode is activated by `structure_phylogeny.py --pseudo-z per_query`.

Let `D_ij` be the RMSD distance between structures `i` and `j`. Define a similarity matrix `S`:

#### Step 1: RMSD -> similarity (`--similarity`)

1. `--similarity neg_rmsd`:

`S_ij = -D_ij`

1. `--similarity exp` (default):

`S_ij = exp(-D_ij / tau)`

where `tau` is provided by `--tau`.

**Math (rendered)**

```
S_ij = -D_ij                  (neg_rmsd)
S_ij = exp( -D_ij / tau )      (exp)
```

**What it means (how to interpret)**

| Similarity choice | Practical interpretation |
| --- | --- |
| `neg_rmsd` | Similarity decreases linearly with RMSD (simple, unbounded). |
| `exp` | Similarity decays smoothly with RMSD; tau sets the length scale (larger tau = slower decay). |

#### Step 2: per-query z-score standardisation

For each query `i`, compute the mean and population standard deviation over all `j != i`:

`mu_i = mean_{j != i}(S_ij)`

`sigma_i = sqrt( mean_{j != i}( (S_ij - mu_i)^2 ) )`

Then define the (asymmetric) query-centric z-score:

`z_ij = (S_ij - mu_i) / sigma_i`

(If `sigma_i` is numerically ~0, the code leaves z-scores at 0 for that row.)

**Math (rendered)**

```
mu_i = mean_{j != i}( S_ij )
sigma_i = sqrt( mean_{j != i}( (S_ij - mu_i)^2 ) )
z_ij = (S_ij - mu_i) / sigma_i
```

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher z_ij | Structure j is unusually similar to query i relative to i’s typical similarities. |
| z_ij ≈ 0 | About average similarity for query i. |
| Negative z_ij | Less similar than average for query i (later clamped away in z_plus). |

#### Step 3: symmetrise the z-score matrix

The code symmetrises z-scores by averaging:

`z^sym_ij = (z_ij + z_ji) / 2`

so the downstream distance matrix is symmetric.

**Math (rendered)**

```
z_sym_ij = (z_ij + z_ji) / 2
```

#### Step 4: z-score -> distance for tree inference (`--zdist`)

First clamp negative evidence:

`z^+_ij = max(0, z^sym_ij)`

Then apply the selected transform:

1. `--zdist inv`:

`d_ij = 1 / (1 + z^+_ij)`

1. `--zdist exp`:

`d_ij = exp( - z^+_ij / z_scale )`

where `z_scale` is `--zscale`.

1. `--zdist maxminus`:

First compute:

`z_max = max_{i<j}( z^+_ij )`

then:

`d_ij = (z_max - z^+_ij) / z_max`

(with a small epsilon guard in the implementation).

This yields `d_ii = 0` and `d_ij > 0` for off-diagonal pairs.

**Math (rendered)**

```
z_plus_ij = max(0, z_sym_ij)

d_ij = 1 / (1 + z_plus_ij)                 (inv)
d_ij = exp( - z_plus_ij / z_scale )         (exp)
d_ij = (z_max - z_plus_ij) / z_max          (maxminus)
```

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Smaller d_ij | More similar structures (higher z). |
| Larger d_ij | Less similar structures (lower z, especially after clamping). |

#### Query-centric ranking (optional)

When `--ranking-csv` is provided, the script also writes a ranking table per structure `i` computed from the *asymmetric* z-scores `z_ij` (before symmetrisation):

`avg_z_i = mean_{j != i}( z_ij )`

`max_z_i = max_{j != i}( z_ij )`

and ranks by descending `(avg_z_i, max_z_i)`.

**Math (rendered)**

```
avg_z_i = mean_{j != i}( z_ij )
max_z_i = max_{j != i}( z_ij )
```

**What it means (how to interpret)**

| Ranking statistic | Practical interpretation |
| --- | --- |
| Higher avg_z_i | Structure i is broadly similar to many others (a “central” representative). |
| Higher max_z_i | Structure i has at least one very close neighbour (strong closest match). |

### 5.3 Derivation / references (why these formulas)

- The DALI query-centric Z-score framework is rooted in Holm & Sander’s distance-matrix alignment framework:
  - Holm, L.; Sander, C. *Protein structure comparison by alignment of distance matrices.* *J. Mol. Biol.* (1993). DOI: `10.1006/jmbi.1993.1489`
- The `exp(-D/tau)` similarity is a common distance-kernel/RBF-style choice:
  - Micchelli, C. *Interpolation of scattered data: Distance matrices and conditionally positive definite functions.* *Constructive Approx.* (1986). DOI: `10.1007/BF01893414`
- Neighbour-joining is the primary tree-building method:
  - Saitou, N.; Nei, M. *The neighbour-joining method: a new method for reconstructing phylogenetic trees.* *Mol. Biol. Evol.* (1987). DOI: `10.1093/oxfordjournals.molbev.a040454`
- The per-query z-score standardisation uses conventional z-score definitions. The final `z^+ -> distance` mapping is a practical heuristic to produce non-negative distances usable by NJ/UPGMA.

---

## 6. dali_score.py

**Function selection / usage**

- Compute a Dali-like structural similarity score from two structures using residue equivalences.
- **Pairwise:** `python dali_score.py pdb_a pdb_b [-a alignment_file] [--chain-a ID] [--chain-b ID] [--dalilite-path DIR] [-o output.csv]`
- **All-vs-all:** `python dali_score.py --all-vs-all dir_or_file [dir_or_file ...] [--filter PATTERN] [--root-only-dirs] [--output-tree FILE] [--output-plot FILE] [--output-matrix FILE] [--output-ranking FILE]`
- **Directory scan (`_collect_pdb_files`):** default is **recursive** (`*.pdb`, `*.cif`, `*.ent` under each directory). With **`--root-only-dirs`**, only the top level of each directory is used (no subfolders).
- Alignment sources (in order): (1) **DaliLite** (if `--dalilite-path` or `DALILITE_HOME` is set by user), (2) `--alignment` file, (3) biotite structural alignment (if installed), (4) sequence-order matching (same residue numbering). When DaliLite is used, its canonical Z-score and RMSD are reported.

**Libraries / functions used (dali_score)**

- **BioPython:** `PDBParser.get_structure()`, structure iteration for CA atoms; `atom.coord`.
- **NumPy:** `np.linalg.norm`, `np.exp`, `np.atleast_1d`; distance matrices for pairwise Cα–Cα distances.
- **biotite (optional):** `PDBFile.read()`, `atom_array[atom_array.atom_name == 'CA']`, `superimpose_structural_homologs()` returning `fixed_indices`, `mobile_indices` for residue equivalences.

### 6.1 Raw Dali score

**Reader-friendly formula**  
The *raw Dali score* measures structural similarity by comparing intramolecular Cα–Cα distances between two structures. For each pair of equivalent residues (i, j) in the aligned core, compute the distance between residues i and j within structure A and within structure B. The score φ(i,j) rewards similar distances and penalizes deviations; an envelope downweights long-range pairs. The total score is the sum over all residue pairs in the core.

**Mathematical definition**

```
S(A,B) = sum_{i in core} sum_{j in core} phi(i,j)
```

where

```
phi(i,j) = (theta - diff(i,j)) * exp( - (d_star_ij / R0)^2 )
```

```
theta = 0.2
R0 = 20 Å
d_star_ij = (d_ij_A + d_ij_B) / 2
diff(i,j) = abs(d_ij_A - d_ij_B) / d_star_ij
```

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher raw score S(A,B) | More internal-distance agreement across the aligned core → more structurally similar by this scoring model. |
| Lower raw score | Less agreement and/or smaller/poorer aligned core. |

Only terms with φ > 0 contribute. d_{ij}^A and d_{ij}^B are intramolecular Cα–Cα distances in structures A and B.

**Implementation (simplified):** For each pair (i,j) in the core, compute d_{ij}^A and d_{ij}^B from CA coordinates. d* = (d_A + d_B)/2; diff = |d_A − d_B|/d*; φ = max(0, (0.2 − diff) × exp(−(d*/20)²)). Sum all φ. Script: `compute_dali_score()`, `_residue_pair_score()`.

### 6.2 Z-score

**Reader-friendly formula**  
The *Z-score* normalises the raw score for protein length so that scores are comparable across pairs of different sizes. It uses empirically fitted mean m(L) and standard deviation σ(L) for random structure pairs of effective length L = √(L_A × L_B).

**Mathematical definition**

```
Z(A,B) = ( S(A,B) - m(L) ) / sigma(L)
```

```
L = sqrt( L_A * L_B )
```

m(L) ≈ 7.95 + 0.71L − 2.59×10⁻⁴ L² − 1.92×10⁻⁶ L³   for L ≤ 400

m(L) = m(400) + (L − 400)   for L > 400

```
sigma(L) = 0.5 * m(L)
```

**What it means (how to interpret)**

| Value pattern | Typical interpretation |
| --- | --- |
| Higher Z-score | Stronger-than-expected similarity for effective length L; generally more confident structural relationship. |
| Z-score near 0 | Similarity near fitted random baseline for that length. |
| Negative Z-score | Below-baseline similarity (often treated as “no meaningful similarity”). |

**Implementation (simplified):** L = sqrt(L_A * L_B); m = polynomial or linear extrapolation; σ = 0.5*m; Z = (S − m) / σ. Script: `compute_z_score()`, `_mean_score_fit()`, `_std_score_fit()`.

### 6.3 DaliLite integration

When **DaliLite** is available (user specifies its installation directory), it is used as the primary alignment source. In this repository, `**dali_score.py` does not rely on bare `dali.pl --pdbfile1/--pdbfile2`** (which often yields empty internal data). It uses the same pattern as DaliLite v5 batch workflows:

1. **Short working directory:** comparisons run under a short temporary directory (typically a randomly named subdirectory of the system temp folder; DaliLite Fortran binaries enforce a maximum path length of about 80 characters).
2. **Staging:** PDB/mmCIF inputs may be copied or rewritten; PDB files without a `HEADER` record can get a dummy `HEADER` so `mkdssp` treats them as PDB, not mmCIF (`_stage_pdb_for_dalilite_import`).
3. `**import.pl`** builds `DAT/*.dat` per structure: `import.pl --pdbfile … --pdbid xxxx --dat DAT/`.
4. `**dali.pl**` pairwise comparison: `dali.pl --cd1 xxxxX --cd2 yyyyY --dat1 DAT` (with `--outfmt` for summary, equivalences, transrot as needed).

The script parses `**dali.pl` text output** for:

- **Z-score** (canonical Dali Z-score, used in preference to the empirical fit)
- **RMSD** of structurally equivalent Cα atoms where reported
- **Residue equivalences** from the structural equivalences block

**Configuration:** set `--dalilite-path` or `**DALILITE_HOME`** to the DaliLite installation root (directory containing `bin/dali.pl`).

**mkdssp (required for `import.pl`):** DaliLite calls `**mkdssp`** using the path in `**bin/mpidali.pm**` (`$DSSP_EXE`; upstream defaults often use `/usr/local/bin/mkdssp`). That file must exist and be executable—otherwise `DAT/*.dat` can be empty and comparisons fail silently. Fix by installing a compatible `**mkdssp**`, symlinking it to the path in `mpidali.pm`, or editing `$DSSP_EXE`. Perl modules under `bin/` (e.g. `FSSP.pm`) must be readable by the user. Optionally set `**MKDSSP**` to the full path of `mkdssp` so this script’s diagnostics match your install.

**Manual all-vs-all** (outside `dali_score.py`): batch-`import.pl` for each model into a shared `DAT/`, build `query.list` of five-character ids (`xxxxX` per `DAT/*.dat`), then `dali.pl --matrix --query query.list --dat1 DAT` from a short working directory path.

### 6.3.1 dalilite_superpose_scores.py

`**dalilite_superpose_scores.py`** runs the DaliLite pairwise path via `dali_score._dalilite_pair_via_dat`, then applies **translation–rotation** from `--outfmt transrot` to write a superposed target PDB (BioPython). It supports **all-vs-all**, pairwise CSV, Z-matrix CSV, Newick tree, optional dendrogram plot, and an optional **matplotlib Z-score heatmap** (`--heatmap` and related options; colour bar label “Dali Z”; implementation delegates to `rmsd_to_csv.plot_heatmap` with `autoscale_positive_offdiag_only=False`).

Requirements and environment variables match `**dali_score.py`** (DaliLite path, mkdssp as above). Use `**--fallback-biotite**` when DaliLite reports no significant hit (often below Dali’s usual reporting threshold, Z ~ 2) but a structural alignment is still desired.

**All-vs-all structure list:** by default, `_collect_pdb_files(..., recursive=False)` is used: only `*.pdb` / `*.cif` / `*.ent` in the **root** of each directory argument. Pass **`--recursive`** to include subdirectories (same recursive behaviour as `dali_score.py` all-vs-all by default).

### 6.3.2 run_all_superpositions.py (Coot LSQ batch)

`**run_all_superpositions.py**` is a small driver that loops over **condition** subdirectories (e.g. `condition_1`) and **tags** (same idea as `--filter=set_a` in the Coot LSQ/SSM scripts), calling `**superimpose_coot_LSQ.py`** with `**--pattern**`, `**--divider=LSQ_**`, per-tag reference globs, and `**CONDITION_PREFIX**`-derived filename tokens. Defaults such as `**REFERENCE_SUBDIR**` and `**REFERENCE_STEM**` are documented in the script and in the main `**README.md**` (superposition section).

### 6.4 All-vs-all and tree output (`dali_score.py`)

When `dali_score.py` is run with `--all-vs-all` and one or more directories or PDB/CIF files, it compares every pair of structures (after `_collect_pdb_files`, recursive unless `--root-only-dirs`), collects Z-scores, and optionally generates:

- **Pairwise table** (`-o`): CSV with pdb_a, pdb_b, raw_score, z_score, n_core, alignment_source, dalilite_rmsd, lali, nres, pct_id, dalilite_hit_id, description (DaliLite summary fields empty when alignment is not from DaliLite)
- **Ranking** (`--output-ranking`, default `dali_ranking.csv`): structures ranked by average and max Z-score
- **Z-score matrix** (`--output-matrix`): symmetric matrix CSV for downstream use
- **Newick tree** (`--output-tree`, default `dali_tree.nwk`): phylogenetic tree from Z-score–derived distances (neighbour-joining or UPGMA)
- **Dendrogram plot** (`--output-plot`): PNG/SVG tree visualisation (requires ete3 or biopython+matplotlib)

**Z-score → distance transform** (`--transform`): Higher Z = more similar. For trees, Z is converted to distance: `inv` (d = 1/(1+Z)), `maxminus` (normalised), or `exp` (d = exp(-Z/scale)). Tree building uses scikit-bio or Biopython for neighbour-joining when available, else a pure-Python UPGMA fallback. See Saitou & Nei (1987) for NJ.

**Filtering:** `--filter` limits which files are included: plain text matches as a substring on the basename; patterns with `*`, `?`, or `[` use `fnmatch` on basename or stem. Implementation: `_collect_pdb_files(paths, filter_pattern, recursive=...)` in `dali_score.py` (`recursive` is true unless `--root-only-dirs`; `dalilite_superpose_scores.py` passes `recursive` according to its `--recursive` flag).

### 6.5 Residue equivalences

**Alignment file format (optional):** TSV or CSV with either:

- Two columns: `resnum_A`, `resnum_B` (for single chain per structure)
- Four columns: `chain_A`, `resnum_A`, `chain_B`, `resnum_B`

Lines starting with `#` are comments. Header row is auto-skipped.

**Automatic alignment:** When no alignment file is provided:

1. **biotite** (if installed): `superimpose_structural_homologs()` returns fixed/mobile indices of equivalent CA atoms; mapped to (chain, resseq, icode).
2. **sequence-order:** For structures sharing residue numbering (e.g. same chain IDs and residue numbers), use 1:1 correspondence by (chain, resseq).

**Derivation / references**

- Holm, L.; Sander, C. *Protein structure comparison by alignment of distance matrices.* *J. Mol. Biol.* **233**, 123–138 (1993). DOI: [10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489) — Original Dali method: raw score φ(i,j), Z-score normalisation, empirical m(L) and σ(L).
- Holm, L. *Dali and the persistence of protein shape.* *Protein Sci.* **29**, 128–140 (2020; first published online 2019). DOI: [10.1002/pro.3749](https://doi.org/10.1002/pro.3749) — DaliLite standalone implementation.
- The Dali score uses intramolecular Cα–Cα distances only; superposition is not required. Equivalences can come from DaliLite, Coot/SSM (export alignment if available), biotite, or any structural alignment tool.

---

## Reference: Default parameters


| Script / class                       | Parameter                                             | Default                         |
| ------------------------------------ | ----------------------------------------------------- | ------------------------------- |
| interface_analyser                   | contact_distance                                      | 5.0 Å                           |
| contact_analyser                     | contact_distance                                      | 4.5 Å                           |
| packing_metrics                      | protein density (for solvent)                         | 0.81 Da/Å³                      |
| interface_analyser                   | BSA (fallback when SASA unavailable)                  | 15 Å² per contact atom          |
| contact_analyser                     | surface neighbour cutoff                               | 5.0 Å; < 12 neighbours → surface |
| contact_analyser                     | estimated area per contact                            | 4 Å²                            |
| contact_analyser                     | avg atomic volume (packing)                           | 20 Å³                           |
| interface_analyser, contact_analyser | H-bond max distance                                   | 3.5 Å                           |
| interface_analyser, contact_analyser | Electrostatic max distance                            | 5.0 Å                           |
| interface_analyser, contact_analyser | Hydrophobic range                                     | 3.5–4.5 Å                       |
| interface_analyser, contact_analyser | Van der Waals max distance                            | 5.0 Å                           |


---

## References (selected)

*Bibliographic metadata and DOIs below were checked against the [Crossref](https://www.crossref.org/) REST API (works endpoint). PubMed-indexed articles align with the same DOIs where applicable.*

- **Holm, L.; Sander, C.** (1993). Protein structure comparison by alignment of distance matrices. *J. Mol. Biol.* **233**, 123–138. DOI: [10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489)
- **Holm, L.** (2020). Dali and the persistence of protein shape. *Protein Sci.* **29**, 128–140 (print; first published online 2019). DOI: [10.1002/pro.3749](https://doi.org/10.1002/pro.3749)
- **Saitou, N.; Nei, M.** (1987). The neighbour-joining method: a new method for reconstructing phylogenetic trees. *Mol. Biol. Evol.* **4**, 406–425. DOI: [10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)
- **Matthews, B.W.** (1968). Solvent content of protein crystals. *J. Mol. Biol.* **33**, 491–497. DOI: [10.1016/0022-2836(68)90205-2](https://doi.org/10.1016/0022-2836(68)90205-2)
- **Shrake, A.; Rupley, J.A.** (1973). Environment and exposure to solvent of protein atoms. *J. Mol. Biol.* **79**, 351–371. DOI: [10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)

