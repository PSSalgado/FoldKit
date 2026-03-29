# Metrics Details: Formulas and Script Calculations

This document lists all metrics computed by the FoldKit analyzer and metrics scripts, with:

**Related (not covered below):** **`trim_models.py`** only harmonizes residue ranges (standalone script; shared with **`trim_superimposeLSQ.py`** trimming); it does not compute packing or similarity metrics. See the main **README.md** (File management).

- **Mathematical definition** and **script calculation** for each metric.
- **Reader-friendly formula**: a short, plain-language description of what is being computed and how.
- **Libraries and functions**: BioPython, NumPy, SciPy, NetworkX, pandas, scikit-learn (and fallbacks) used in the implementation.
- **Derivation / references** where applicable.

---

## Libraries used (summary)


| Script                  | Libraries                                                                                                                        | Main functions used                                                                                                                                                                                                                                                                                                                                                                                    |
| ----------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| packing_metrics.py      | **BioPython** (Bio.PDB.PDBParser, PDBIO, PDBExceptions), **NumPy**, **SciPy** (optional)                                         | `PDBParser.get_structure()`, manual CRYST1 read; `structure` → model → chain → residue → atom; `atom.coord`, `atom.element`, `atom.name`; `np.radians`, `np.sqrt`, `np.cos`; `scipy.spatial.distance.pdist()` or fallback loop with `np.linalg.norm`                                                                                                                                                   |
| interface_analyzer.py   | **BioPython** (PDBParser, NeighborSearch, Structure, Model, PDBExceptions), **Bio.PDB.SASA** (ShrakeRupley, optional), **NumPy** | `parser.get_structure()`, `structure.get_chains()`, `chain.get_atoms()`; `NeighborSearch(atom_list).search_all(threshold, 'R')` → interface residues; contact residues = interface residues filtered by same criteria (distance + type) via `_residue_passes_contact_criteria()`; `ShrakeRupley().compute()` for BSA; `np.linalg.norm`, `np.mean`, etc.                                                |
| contact_analyzer.py     | **BioPython** (PDBParser, PDBExceptions), **NumPy**                                                                              | Same structure traversal; `chain.get_atoms()`, `atom.coord`; `np.linalg.norm`; manual contact loop and type classification                                                                                                                                                                                                                                                                             |
| channel_analyzer.py     | **BioPython** (PDBParser), **NumPy**                                                                                             | `parser.get_structure()`, `structure.get_atoms()`, `atom.coord`; `np.array`, `np.min`/`np.max` (axis=0), `np.random.uniform`, `np.linalg.norm`, `np.prod`, `np.std`, `np.mean`, `np.sort`, `np.random.choice`                                                                                                                                                                                          |
| graph_analyzer.py       | **BioPython** (PDBParser), **NumPy**, **NetworkX**                                                                               | `parser.get_structure()`, `structure.get_chains()`/`get_residues()`, `residue.get_atoms()`; `nx.Graph`, `nx.density`, `nx.average_shortest_path_length`, `nx.diameter`, `nx.average_clustering`, `nx.betweenness_centrality`, `nx.closeness_centrality`, `nx.degree_centrality`, `nx.algorithms.community.greedy_modularity_communities`, `nx.algorithms.community.modularity`, `nx.erdos_renyi_graph` |
| comparative_analyzer.py | **NumPy**, **pandas**, **scikit-learn**                                                                                          | `pd.DataFrame()`, `df.select_dtypes()`, `np.mean`/`np.std`/`np.min`/`np.max`/`np.median`, `np.corrcoef`, `StandardScaler().fit_transform()`, `KMeans()`, `PCA().fit_transform()`, `np.cumsum`, `np.abs`, `np.where`                                                                                                                                                                                    |
| dali_score.py          | **BioPython** (PDBParser), **NumPy**, **biotite** (optional)                                                                      | `PDBParser.get_structure()`, structure traversal for CA atoms; `np.linalg.norm`, `np.exp` for distance matrices and Dali formula; `biotite.structure.superimpose_structural_homologs()` for automatic residue equivalences when available                                                                                                                          |


---

## 1. packing_metrics.py

**Function selection / usage**

- These metrics are computed in `packing_metrics.py`, typically called by `crystal_packing_analyzer.py` for each structure.
- At the CLI level, you select whether to run packing metrics by:
  - Running `packing_metrics.py` directly on one or more PDB/CIF files, or
  - Enabling the “packing” stage in `crystal_packing_analyzer.py` (it runs packing by default unless explicitly skipped in that orchestrator).
- Within `packing_metrics.py`, the main driver function (e.g. `analyze_packing_for_structure`) calls:
  - `_calculate_unit_cell_volume()` → §1.1
  - `_analyze_molecular_content()` → §1.2
  - `_calculate_matthews_metrics()` → §1.3
  - `_calculate_packing_density()` and `_analyze_void_spaces()` → §1.4–1.5

### 1.1 Unit cell volume

**Reader-friendly formula**

*Unit cell volume* is the volume of the crystallographic box. For a general box with side lengths *a*, *b*, *c* and angles α, β, γ between them, the volume is: *(a × b × c)* times a geometric factor that depends only on the three angles. The script reads *a*, *b*, *c*, α, β, γ from the PDB `CRYST1` line, converts angles to radians, then applies this formula.

**Mathematical definition**

For a general parallelepiped with edge lengths a, b, c (Å) and angles \alpha, \beta, \gamma (degrees):

V = abc\sqrt{1 + 2\cos\alpha\cos\beta\cos\gamma - \cos^2\alpha - \cos^2\beta - \cos^2\gamma}

**Implementation (simplified):** Read a, b, c (Å) and angles α, β, γ (degrees) from the PDB CRYST1 line. Convert angles to radians. Volume = a × b × c × sqrt(1 + 2×cos(α)×cos(β)×cos(γ) − cos²(α) − cos²(β) − cos²(γ)).

**Script calculation**

- Unit cell parameters are read from the PDB `CRYST1` line (defaults: a=b=c=1, α=β=γ=90°, space_group P1 if missing).
- Angles are converted to radians, then:
  - `volume = a * b * c * np.sqrt(1 + 2*cos(α)*cos(β)*cos(γ) - cos(α)**2 - cos(β)**2 - cos(γ)**2)`
- Implemented in `_calculate_unit_cell_volume()` (and similarly in `contact_analyzer._calculate_cell_volume()`).

**Libraries / functions used**

- **Input:** Unit cell is read manually by parsing the file for a line starting with `CRYST1` and slicing fields (no BioPython for CRYST1). Structure is loaded with **Bio.PDB.PDBParser** → `parser.get_structure('crystal', pdb_file)` for use in later steps.
- **NumPy:** `np.radians()` to convert angles, `np.cos()`, `np.sqrt()` for the volume expression.

---

### 1.2 Molecular content

**Reader-friendly formula**  
*Molecular content* is a set of counts and one sum: how many atoms, residues, and chains; molecular weight (sum of atomic masses over all atoms); and how many atoms per element. The script walks the structure model → chain → residue → atom and uses each atom’s element (or name) to look up mass and to update counts.

**Implementation (simplified):** Loop over all models → chains → residues → atoms. For each atom: increment total_atoms; add atomic_weights[element] to molecular_weight (use 12.011 for unknown); increment element_composition[element]. Count chains and residues. Result: total_atoms, molecular_weight, element_composition, chain_count, residue_count.


| Metric                  | Mathematical definition                         | Script calculation                                                      |
| ----------------------- | ----------------------------------------------- | ----------------------------------------------------------------------- |
| **total_atoms**         | Count of all atoms in the structure             | Loop over all models → chains → residues → atoms; increment counter.    |
| **molecular_weight**    | M = \sum_{\text{atoms}} m_{\text{element}} (Da) | Sum of `atomic_weights[element]` per atom (default 12.011 for unknown). |
| **element_composition** | Count per element                               | Dict of element → count over all atoms.                                 |
| **chain_count**         | Number of chains                                | Count of chains across all models.                                      |
| **residue_count**       | Number of residues                              | Count of residues across all models.                                    |


- Implemented in `_analyze_molecular_content()`.
- Atomic weights (Da): C 12.011, N 14.007, O 15.999, S 32.06, P 30.974, H 1.008, etc.

**Libraries / functions used**

- **BioPython:** `parser.get_structure()` returns a `Structure`; iteration `for model in structure`, `for chain in model`, `for residue in chain`, `for atom in residue`. Per atom: `atom.element` (string), `atom.name` (fallback if element missing). No PDBIO/Select used for this metric.
- **Script:** Built-in dict `atomic_weights` keyed by element symbol; default carbon mass for unknown.

---

### 1.3 Matthews coefficient and solvent content

**Reader-friendly formula**  
*Matthews coefficient* is unit-cell volume per dalton of protein (Å³/Da). *Protein volume* is estimated as mass divided by a fixed density (≈0.81 Da/Å³). *Solvent volume* is unit-cell volume minus protein volume. *Solvent content* is that solvent volume as a percentage of unit-cell volume.

**Mathematical definition**

**Matthews coefficient** (Å³/Da): V_{\mathrm{m}} = \frac{V_{\mathrm{cell}}}{M}

**Estimated protein volume** (Å³): V_{\mathrm{protein}} = \frac{M}{\rho},\quad \rho \approx 0.81\mathrm{Da/\AA^3}

**Solvent volume** (Å³): V_{\mathrm{solvent}} = V_{\mathrm{cell}} - V_{\mathrm{protein}}

**Solvent content** (%): \text{solvent} = \frac{V_{\mathrm{solvent}}}{V_{\mathrm{cell}}} \times 100

**Implementation (simplified):** Matthews = unit_cell_volume / molecular_weight (0 if molecular_weight ≤ 0). Protein_volume = molecular_weight / 0.81 (Å³). Solvent_volume = unit_cell_volume − protein_volume. Solvent_content = (solvent_volume / unit_cell_volume) × 100, clamped to [0, 100].

- **Script:** In `_calculate_matthews_metrics()`.

**Libraries / functions used**

- **Script only:** Arithmetic using `unit_cell_volume` and `molecular_weight` from previous steps; no extra BioPython or NumPy calls in this function.

**Derivation / references**

- The **Matthews coefficient** is defined as the ratio of unit-cell volume to protein mass:
  - B. W. Matthews, “Solvent content of protein crystals”, *J. Mol. Biol.* **33**, 491–497 (1968). DOI: [10.1016/0022-2836(68)90205-2](https://doi.org/10.1016/0022-2836(68)90205-2)
- The typical protein partial specific volume (≈0.73–0.75 mL/g) and average crystal solvent content (40–60%) imply:
  - \rho \approx 1.35\ \mathrm{g/cm^3} \approx 0.81\ \mathrm{Da/\AA^3}, giving the formula for V_\mathrm{protein} used here.
- The solvent fraction follows directly from a volume balance:
  - V_\mathrm{cell} = V_\mathrm{protein} + V_\mathrm{solvent}, so \text{solvent} = V_\mathrm{solvent}/V_\mathrm{cell}.

---

### 1.4 Packing density and void fraction

**Reader-friendly formula**  
*Total atomic volume* is the sum of a fixed volume per atom type (from van der Waals–style tables). *Packing density* is that sum divided by unit-cell volume (what fraction of the cell is “filled” by atoms). *Packing efficiency* is the same as packing density expressed as a percentage. *Void fraction* is one minus packing density (the unfilled fraction).

**Mathematical definition**

**Total atomic volume** (Å³): V_{\mathrm{atoms}} = \sum_{\text{atoms}} V_{\mathrm{element}}

**Packing density** (dimensionless): \text{packingdensity} = \frac{V_{\mathrm{atoms}}}{V_{\mathrm{cell}}}

**Packing efficiency** (%): \text{packingefficiencypercent} = \text{packingdensity} \times 100

**Void fraction:** \text{voidfraction} = 1 - \text{packingdensity}

**Implementation (simplified):** For each atom, look up atomic_volumes[element] (default 16.44 Å³ for C). Total_atomic_volume = sum over all atoms. Packing_density = total_atomic_volume / unit_cell_volume (0 if volume ≤ 0). Packing_efficiency_percent = packing_density × 100. Void_fraction = 1 − packing_density.

**Libraries / functions used**

- **BioPython:** Same structure traversal (model → chain → residue → atom); `atom.element`, `atom.name` for element lookup.
- **Script:** Dict `atomic_volumes` by element; default carbon volume.
- **NumPy:** Used only for results (no np calls in this block beyond possible type handling).

**Derivation / references**

- **Atomic volumes** are approximated from empirical van der Waals volumes, e.g.:
  - A. Bondi, “van der Waals Volumes and Radii”, *J. Phys. Chem.* **68**, 441–451 (1964). DOI: [10.1021/j100785a001](https://doi.org/10.1021/j100785a001); and later tabulations used in protein packing analyses.
- **Packing density** is the fractional occupancy of the unit cell by atom-based spheres:
  - V_\mathrm{atoms} = \sum V_\mathrm{element}, so \phi = V_\mathrm{atoms}/V_\mathrm{cell}.
  - \phi close to 0.7–0.8 is typical for tightly packed proteins; higher `void_fraction` indicates looser packing.

---

### 1.5 Void / interatomic metrics (_analyze_void_spaces)

**Reader-friendly formula**  
Each atom is assigned a *radius* from its element’s volume assuming a sphere (V = ⁴⁄₃πr³ → r = (3V/(4π))^(1/3)). *Average interatomic distance* is the mean of all pairwise distances between atoms; *minimum* is the smallest such distance. *Atom density* is number of atoms divided by unit-cell volume. The script collects all atom coordinates, then either uses SciPy’s `pdist` for pairwise distances or a double loop with `np.linalg.norm`.

**Mathematical definition**

**Atomic radius from volume** (sphere model V = \frac{4}{3}\pi r^3): r = \left( \frac{3V}{4\pi} \right)^{1/3}

**Average interatomic distance** (Å): \bar{d} = \frac{1}{|\mathcal{P}|} \sum_{(i,j) \in \mathcal{P}} \mathbf{x}_i - \mathbf{x}_j, where \mathcal{P} is the set of all unordered atom pairs.

**Minimum interatomic distance** (Å): d_{\min} = \min_{(i,j)} \mathbf{x}_i - \mathbf{x}_j

**Atom density** (atoms/Å³): \rho_{\mathrm{atom}} = \frac{N_{\mathrm{atoms}}}{V_{\mathrm{cell}}}

**Implementation (simplified):** Per atom, radius = (3 × atomic_volumes[element] / (4π))^(1/3). Collect all atom coordinates; compute all pairwise distances (SciPy pdist or double loop). Average interatomic distance = mean(distances); minimum = min(distances). Atom density = number of atoms / unit_cell_volume (0 if no atoms or volume ≤ 0).

**Libraries / functions used**

- **BioPython:** Structure iteration; `atom.coord` (NumPy-like array), `atom.element`, `atom.name`.
- **NumPy:** `np.array(coords)`, `np.pi`; `np.mean(distances)`, `np.min(distances)`; fallback: `np.linalg.norm(coords[i] - coords[j])` in a double loop.
- **SciPy (optional):** `scipy.spatial.distance.pdist(coords)` to compute all pairwise distances in one call; script falls back to the loop if SciPy is not available.

**Derivation / references**

- The atomic radius from volume assumes a **spherical atom model**: V = \frac{4}{3}\pi r^3 → r = (3V/4\pi)^{1/3}.
- Average and minimum interatomic distances are standard descriptors of packing tightness and are analogous to nearest-neighbor statistics in condensed-matter physics.
- Atom-density measures mirror typical crystallographic discussions of packing fraction and nearest-neighbor spacing in protein crystals.

---

## 2. interface_analyzer.py

**Function selection / usage**

- Interface metrics are computed in `interface_analyzer.py`, usually via a top-level function like `analyze_interfaces(structure)` that is called by `crystal_packing_analyzer.py`.
- At the CLI level, you select these metrics by:
  - Running `interface_analyzer.py` directly on a PDB/CIF (specifying chains or letting the script enumerate interfaces), or
  - Allowing `crystal_packing_analyzer.py` to run the interface stage on each structure.
- Within `interface_analyzer.py`, the main pipeline:
  - **Interface residues:** `NeighborSearch(atom_list).search_all(contact_distance, 'R')` → residue pairs within threshold, different chains.
  - **Contact residues:** Those interface residues filtered by applying the same criteria (distance + type-specific limits) as the contact list: keep residue iff it has ≥1 atom that forms a pair with the other chain satisfying those criteria.
  - Calls `_find_contacts()` and `_filter_contacts_by_limits()` → §§2.1–2.4, 2.6.
  - **BSA:** `_calculate_buried_surface_area(chain1, chain2, model)` using `Bio.PDB.SASA.ShrakeRupley` when available (§2.2); else fallback estimate.
  - Calls `_calculate_contact_area()` and `_calculate_interface_complementarity()` → §2.3–2.5, 2.7.
- **Contact distance:** Default 5.0 Å; atom pairs with distance ≤ this are candidates; only those passing type-specific limits (below) are counted.

**Type-specific distance limits (contacts outside these are ignored):**


| Type          | Limit (Å) |
| ------------- | --------- |
| H-bonds       | ≤ 3.5     |
| Electrostatic | ≤ 5.0     |
| Hydrophobic   | 3.5–4.5   |
| Van der Waals | ≤ 5.0     |


- Contacts are classified by atom types (N/O/S–N/O/S → H-bond; C–C → hydrophobic; N/O/S/P → electrostatic; else van der Waals), then only kept if distance is within the limit for that type. All downstream metrics (BSA, contact area, complementarity) use only these filtered contacts.

**Libraries / functions used (interface_analyzer)**

- **BioPython:** `PDBParser.get_structure('crystal', pdb_file)`; `structure.get_chains()` → list of chains; `chain.get_atoms()` for each chain. **Interface residues:** `NeighborSearch(atom_list).search_all(threshold, 'R')` → residue pairs within threshold, different chains. **Contact residues:** interface residues filtered by applying the same criteria (distance + type-specific limits) to each residue’s atoms vs the other chain (`_residue_passes_contact_criteria`). **BSA:** when `Bio.PDB.SASA.ShrakeRupley` is available: `Structure`, `Model` to build single-chain structures; `ShrakeRupley(...).compute(entity, level='S')`; BSA = SASA(chain1 alone) + SASA(chain2 alone) − SASA(model). Fallback: per-atom estimate (~15 Å² per contact atom). Per atom: `atom.coord`, `atom.element`, `atom.name`, `atom.parent`.
- **NumPy:** `np.linalg.norm()`, `np.mean()`, `np.min()`, `np.max()`, `np.std()` on distance lists.

### 2.1 Contacts

**Reader-friendly formula**  
A *contact* is a pair of atoms (one from each chain) whose distance is within the global cutoff (5.0 Å) and also within the limit for their contact type (e.g. H-bond ≤ 3.5 Å). **Interface residues** are identified by NeighborSearch: residue pairs from different chains that have at least one atom–atom pair within the distance threshold. **Contact residues** are those interface residues that pass the same criteria as the contact list: the residue must have at least one atom that, with some atom on the other chain, satisfies both the distance cutoff and the type-specific distance limit for that pair. Filtering is done by applying these criteria (not by membership in the contact list). The *contact count* is the number of atom pairs that pass the criteria.

**Contact count**

**Definition:** Number of atom pairs (chain1, chain2) with \left\mathbf{x}_1 - \mathbf{x}*2\right \le d*{\mathrm{contact}} **and** within the type-specific distance limit for their contact type.

**Implementation (simplified):** For each atom pair (one from each chain), compute distance. If distance ≤ contact_distance (5 Å), classify the pair by atom types (N/O/S–N/O/S → H-bond; C–C → hydrophobic; N/O/S/P → electrostatic; else van der Waals). Keep the contact only if distance is within the limit for that type (e.g. H-bond ≤ 3.5 Å, electrostatic ≤ 5 Å). Contact count = number of kept pairs. Summary counts per type are reported.

- **Script:** `_find_contacts()` builds raw list; `_filter_contacts_by_limits()` classifies each contact and keeps only those with `_contact_within_limit(contact_type, distance)`.  
- **Output:** Summary counts per type (H-bonds, electrostatic, hydrophobic, van der Waals); no full contact list.

### 2.2 Buried surface area (BSA)

**Reader-friendly formula**  
*Buried surface area* is the surface area of each molecule that becomes hidden when the two chains form the interface. The script follows the **Biopython Interface Analysis** approach ([https://biopython.org/wiki/Interface_Analysis](https://biopython.org/wiki/Interface_Analysis)): when `Bio.PDB.SASA.ShrakeRupley` is available, BSA = SASA(chain1 alone) + SASA(chain2 alone) − SASA(complex). SASA is computed with the Shrake–Rupley rolling-ball algorithm (probe 1.4 Å). Single-chain structures are built with `Structure`/`Model` and deep-copied chains so that SASA in isolation is well-defined. If SASA is not available (or computation fails), BSA is estimated as the number of atoms in contact (within contact_distance) × 15 Å².

**Mathematical definition**

\text{BSA} = \text{SASA}(\text{chain}_1) + \text{SASA}(\text{chain}_2) - \text{SASA}(\text{complex})

(when using ShrakeRupley); otherwise \text{BSA} \approx |\text{contact atoms}| \times 15\mathrm{\AA^2} (fallback).

**Implementation (simplified):** When ShrakeRupley is available: compute SASA of chain1 alone, chain2 alone, and the complex; BSA = SASA(chain1) + SASA(chain2) − SASA(complex). Otherwise: count atoms in contact (within contact_distance) and multiply by 15 Å² per atom.

- **Script:** If `Bio.PDB.SASA.ShrakeRupley` is available: build single-chain structures (deepcopy chain into new `Structure`/`Model`), run `ShrakeRupley().compute(..., level='S')` on each and on the full model; BSA = sasa1 + sasa2 − sasa_complex. Else: `_estimate_buried_surface_area_fallback()` using atoms within `contact_distance` × 15 Å².  
- In `_calculate_buried_surface_area()` (and `_estimate_buried_surface_area_fallback()`).

### 2.3 Contact area

**Reader-friendly formula**  
*Contact area* is the total area attributed to the interface, computed per contact: if two atoms (with radii r₁, r₂) are closer than r₁ + r₂, they “overlap”; the script assigns an area to that overlap using a circular-cap formula: π × (smaller radius)² × (1 − d/(r₁+r₂)), then sums over all contacts.

**Definition:** Total area attributed to contact patches between the two chains.

**Mathematical definition:** For each contact with distance d and radii r_1, r_2: if d < r_1 + r_2 then
\text{area} = \pi \min(r_1,r_2)^2 \left(1 - \frac{d}{r_1+r_2}\right), clamped to \ge 0. Contact area = sum over contacts.

**Implementation (simplified):** For each contact, get distance d and element-based radii r₁, r₂ (default 1.7 Å). If d < r₁ + r₂, add area = π × (min(r₁,r₂))² × (1 − d/(r₁+r₂)), clamped to ≥ 0. Contact area = sum over all contacts.

- **Script:** In `_calculate_contact_area()`.

### 2.4 Average / min / max contact distance and std

**Reader-friendly formula**  
These are standard statistics over the list of contact distances (after type filtering): mean, min, max, and standard deviation.

**Implementation (simplified):** Take the list of contact distances (one per atom pair kept after type filtering). Average = mean of list; min = smallest; max = largest; std = standard deviation. Implemented with NumPy: `np.mean`, `np.min`, `np.max`, `np.std`.

- **average_contact_distance:** `np.mean(distances)` over all contact distances.  
- **min_contact_distance:** `np.min(distances)`.  
- **max_contact_distance:** `np.max(distances)`.  
- **contact_distance_std:** `np.std(distances)`.  
- In `_analyze_contacts()`.

### 2.5 Interface complementarity

There are two related complementarity scores, both in [0, 1].

#### 2.5.1 Shape complementarity

**Reader-friendly formula**  
*Shape complementarity* measures how “tight” the interface is in terms of contact distances. An ideal contact distance is taken as 3.5 Å. For each contact, the deviation from 3.5 Å is computed; the score is max(0, 1 − average_deviation/2), so smaller average deviation gives a higher score.

**Definition:** Score in [0, 1]; higher = tighter fit. Based on deviation from an “optimal” contact distance.

**Implementation (simplified):** Ideal distance = 3.5 Å. For each contact, 

deviation = |distance − 3.5|. 

**Mathematical definition:** For each contact distance d, deviation \delta = |d - 3.5|; average deviation \bar{\delta}. Then \text{complementarity}_\text{shape} = \max(0, 1 - \bar{\delta}/2). So smaller average deviation gives a higher score.

- In `_calculate_interface_complementarity()`; reported as `interface_complementarity`.

#### 2.5.2 Charge complementarity

**Reader-friendly formula**  
*Charge complementarity* measures how often charged contact pairs are *attractive* (opposite-signed) vs *repulsive* (same-signed). The score is the fraction of charged–charged contacts that are opposite-signed; value in [0, 1], higher = more electrostatically complementary.

**Mathematical definition:** Let N_\text{opp} be the number of charged–charged contacts where the residues have opposite sign, and N_\text{same} the number where they have the same sign. Neutral residues are ignored. Then:
\text{complementarity}*{\text{charge}} =
\begin{cases}
0 & \text{if } N*{\text{opp}} + N_{\text{same}} = 0
\dfrac{N_{\text{opp}}}{N_{\text{opp}} + N_{\text{same}}} & \text{otherwise}
\end{cases}

**Implementation (simplified):** For each contact, get the charge sign of the two residues: ASP/GLU = −1, ARG/LYS/HIS = +1, else 0. If either residue is neutral, skip the contact. If both are charged, count **opposite** when q1 × q2 < 0, **same** when q1 × q2 > 0. Score = opposite / (opposite + same), or 0 if there are no charged–charged contacts. Result is a fraction in [0, 1].

Residue charges (same set as polarity summary):  

- Negative: ASP, GLU → −1  
- Positive: ARG, LYS, HIS → +1  
- All other residues → 0 (ignored in this metric).

**Scope:** Counts are over *all* atom–atom contacts in the interface (each contact is one pair of atoms within distance; the same residue pair can contribute many contacts). So the score can take any value in [0, 1]. It will be exactly **0** when there are no charged–charged contacts or when all such contacts are same-sign (e.g. LYS–ARG). It will be exactly **1** when all charged–charged contacts are opposite-sign (e.g. only ASP–ARG, GLU–LYS). Uneven polarity counts (e.g. 21 vs 6 charged residues) do not force a fractional score: which residue pairs actually form contacts and their signs determine the fraction.

- Implemented in `_calculate_charge_complementarity()`; reported as `charge_complementarity`. Diagnostic counts: `charge_complementarity_opposite`, `charge_complementarity_same`, `charge_complementarity_total` (and in the text report when total > 0).

#### 2.5.3 Charge complementarity density (per interface area)

**Reader-friendly formula**  
*Charge complementarity density* is the number of opposite-sign charged–charged contacts per unit interface area (Å⁻²). It combines “how complementary” the interface is (in terms of count of attractive charged contacts) with interface size, so larger interfaces are comparable to smaller ones. Higher values mean more complementary charged contacts per Å².

Definition: (number of opposite-sign charged–charged contacts) / area
with area = contact_area if contact_area ≥ 0.01 Å², else BSA.
Unit: per Å² (e.g. 0.05 → 0.05 opposite-sign contacts per Å²).
Also stored: charge_complementarity_density_denominator = 'contact_area' or 'buried_surface_area' so you know which area was used.
If both contact area and BSA are 0, charge_complementarity_density is None.

It's a size-normalized measure: interfaces with the same fraction of opposite-sign contacts can be distinguished by how many such contacts there are per Å² (e.g. 2 vs 20 opposite contacts on a 100 Å² interface).

**Why contact area vs BSA:** The denominator is **contact area** when it is ≥ 0.01 Å² (same contact-based definition as the contact list and contact area). When contact area is zero or negligible, **buried surface area (BSA)** is used so the metric remains defined. Contact area is preferred because it is derived from the same atom–atom contacts that define the charged-contact counts; BSA is a different measure (solvent-accessible surface buried) but is a sensible fallback when there is no overlap-based contact area.

**Mathematical definition**

\text{chargecomplementaritydensity} = \frac{N_{\text{opp}}}{\text{area}},\quad \text{area} = \begin{cases} \text{contactarea} & \text{if contactarea} \ge 0.01 \text{BSA} & \text{otherwise} \end{cases}

(undefined / not reported when both contact_area and BSA are 0.)

**Implementation (simplified):** After computing charge complementarity, set area = contact_area if contact_area ≥ 0.01, else buried_surface_area. If area > 0, charge_complementarity_density = opposite / area (units: per Å²). Also store charge_complementarity_density_denominator = `'contact_area'` or `'buried_surface_area'` to record which was used. If area is 0, report `charge_complementarity_density` = None.

- Reported as `charge_complementarity_density` (per Å²) and `charge_complementarity_density_denominator`; in the text report when density is not None.

### 2.6 Contact type counts (incl. H-bonds)

**Implementation (simplified):** When filtering contacts by type and distance, each kept contact is classified as one of: hydrogen_bond (N/O/S–N/O/S, ≤3.5 Å), electrostatic (N/O/S/P, ≤5 Å), hydrophobic (C–C, 3.5–4.5 Å), van_der_waals (else, ≤5 Å). Count how many contacts fall into each class. `hydrogen_bonds` = count for hydrogen_bond.

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

f_{\mathrm{class}} = \frac{\text{count}_{\mathrm{class}}}{|\mathcal{R}_c|}

when |\mathcal{R}_c| > 0.

**Accessibility**

Per-residue SASA is computed with Shrake–Rupley (`ShrakeRupley(probe_radius=1.4, n_points=100).compute(model, level='R')`). For contact residues on each chain:

- **averagesasachainX:** mean SASA over contact residues with SASA information.  
- **accessiblefractionchainX:** fraction of those residues with SASA > 0.

Formally, for SASA values s_i over contact residues on a chain with available SASA:

\text{averagesasa} = \frac{1}{n}\sum_{i=1}^{n} s_i,\qquad
\text{accessiblefraction} = \frac{|i : s_i > 0|}{n}

(0 if no SASA values are available).

**Script calculation**

- **Polarity:** `_classify_residue_polarity()` maps residue names into classes; `_summarize_polarity_and_accessibility()` accumulates `polarity_chain1`, `polarity_chain2` (counts) and `polarity_fractions_chain1`, `polarity_fractions_chain2` (fractions).  
- **Accessibility:** `_compute_residue_sasa_for_chains()` runs Shrake–Rupley with `level='R'` on a deepcopy of the model and records per-residue SASA for the two chains; `_summarize_polarity_and_accessibility()` then computes `average_sasa` and `accessible_fraction` for contact residues per chain, returned as `accessibility_chain1` and `accessibility_chain2`.

If `Bio.PDB.SASA.ShrakeRupley` is not available or SASA computation fails, accessibility fields default to 0 and only polarity counts/fractions are reported.

### 2.8 Interface RMSD (contact residues, CA atoms)

**Reader-friendly formula**  
For each interface between two chains, a simple **interface RMSD** is computed by superposing the CA atoms of **matching contact residues** on the two chains and reporting the root-mean-square deviation (RMSD, Å). Matching is based on residue name and residue number (and insertion code) among residues that are in the contact-residue sets on both chains. This is most meaningful for homodimeric interfaces where the two chains are sequence-identical and aligned.

**Mathematical definition**

Let \mathcal{M} be the set of matched contact residues between chain 1 and chain 2. For each matched pair i \in \mathcal{M}, take the CA coordinates \mathbf{x}_i (chain 1) and \mathbf{y}*i (chain 2). After optimal rigid-body superposition (rotation R, translation \mathbf{t}) that minimizes RMSD:
\mathrm{RMSD} = \sqrt{\frac{1}{|\mathcal{M}|} \sum*{i \in \mathcal{M}} \left R\mathbf{y}_i + \mathbf{t} - \mathbf{x}_i \right^2}

**Implementation (simplified):** Match contact residues between the two chains by (resname, residue number, insertion code). For each matched pair, take the CA atom coordinates. If there are at least 3 matched CA pairs, use BioPython Superimposer to find the rotation and translation that minimize RMSD; RMSD = sqrt(mean of squared distances between superposed CA pairs). If fewer than 3 pairs or superposition unavailable, report N/A (or 0 in the summary).

**Script calculation**

- Collect `contact_residues_chain1` / `contact_residues_chain2` (as residue objects).  
- Build matches by residue `(resname, residue_number, insertion_code)` across the two contact-residue sets.  
- For each matched pair, extract CA atoms; skip pairs lacking CA.  
- If there are at least 3 CA pairs:
  - Use Biopython `Superimposer` (`Bio.PDB.Superimposer`) to compute the optimal superposition and RMSD.  
  - Store per-interface value as `interface_rmsd_ca`.
- In the structure-level summary, **average interface RMSD** is:

 \text{averageinterfacermsdca} =
\frac{1}{N_\mathrm{rmsd}} \sum_{\text{interfaces with valid RMSD}} \mathrm{RMSD}_k

where N_\mathrm{rmsd} is the number of interfaces with RMSD > 0. This value is exported to the comparative pipeline as `average_interface_rmsd_ca`.

If there are no valid matches or the superposition fails, `interface_rmsd_ca` and the summary’s `average_interface_rmsd_ca` default to 0.

### 2.9 Summary (per structure)

- **total_buried_surface_area:** Sum of BSA over all chain-pair interfaces.  
- **total_contact_area:** Sum of contact_area over interfaces.  
- **average_buried_area_per_interface:** `total_buried_surface_area / interface_count`.  
- **average_contact_area_per_interface:** `total_contact_area / interface_count`.  
- **average_interface_rmsd_ca:** Mean CA-atom RMSD over interfaces with a valid RMSD estimate.  
- **total_interfaces:** Number of chain pairs with at least one contact.

**Derivation / references**

- Distance cutoffs and contact-type classifications are inspired by typical protein interface analyses, e.g.:
  - P. Chakrabarti & J. Janin, “Dissecting protein–protein recognition sites”, *Proteins* **47**, 334–343 (2002). DOI: `10.1002/prot.10085`.
  - J. Janin & B. Séraphin, “Genome-wide studies of protein–protein interaction”, *Curr. Opin. Struct. Biol.* **13**, 383–388 (2003). DOI: `10.1016/S0959-440X(03)00064-2`.
- **Interface definition and BSA** follow the Biopython Interface Analysis wiki ([https://biopython.org/wiki/Interface_Analysis](https://biopython.org/wiki/Interface_Analysis)): interface residues from NeighborSearch residue pairs within threshold; BSA = SASA(chain1) + SASA(chain2) − SASA(complex). The implementation uses `Bio.PDB.SASA.ShrakeRupley` (Shrake & Rupley, *J. Mol. Biol.* **79**, 351–371, 1973). DOI: `10.1016/0022-2836(73)90011-9` when available; otherwise a 15 Å² per contact-atom fallback.
- The simple complementarity score is a normalized measure of how close contact distances are to an “ideal” 3.5 Å, analogous to empirical interface “tightness” measures used in docking and packing quality checks.

---

## 3. contact_analyzer.py

**Function selection / usage**

- Contact metrics are computed in `contact_analyzer.py`, which can be:
  - Run directly on a PDB/CIF to analyze ASU contacts and crude crystal contacts, or
  - Invoked by `crystal_packing_analyzer.py` for each structure.
- The main driver function (e.g. `analyze_contacts_for_structure`) calls:
  - `_find_chain_contacts()` → ASU contacts (§3.1).
  - `_calculate_contact_density()` → §3.2.
  - `_estimate_crystal_contacts()` and `_calculate_cell_volume()` → §§3.3–3.4.
- **Contact distance:** Default 4.5 Å (candidate pairs); only contacts within the same **type-specific distance limits** as interface_analyzer (H-bonds ≤3.5 Å, electrostatic ≤5 Å, hydrophobic 3.5–4.5 Å, van der Waals ≤5 Å) are kept.

**Libraries / functions used (contact_analyzer)**

- **BioPython:** `PDBParser.get_structure()`, `structure.get_chains()`, `chain.get_atoms()`; `atom.coord`, `atom.element`, `atom.name`, `atom.parent` (residue). Distance: `np.linalg.norm(atom.coord - other_atom.coord)`.
- **NumPy:** `np.linalg.norm`, `np.mean`, `np.min`, `np.max` on distance lists; sum over chain atom counts for contact density.

### 3.1 ASU contacts

**Reader-friendly formula**  
*ASU contacts* are atom pairs from different chains within the ASU that are within the global distance (4.5 Å) and within the type-specific limit for their classified interaction type. Each contact stores chains, residues, atoms, distance, and type. Totals and averages are taken over this filtered list.

**Implementation (simplified):** For each atom pair from different chains, compute distance. If distance ≤ 4.5 Å, classify type (H-bond, electrostatic, hydrophobic, van der Waals) from atom types; keep only if distance is within the limit for that type. Each kept pair is one contact. total_contacts = len(contact list); average/min/max distance and counts per type from this list.

- **Contacts:** Pairs of atoms from different chains with distance ≤ `contact_distance` and within the type-specific limit for their classified type; each contact stores chain ids, residues, atoms, distance, and contact_type.  
- **Script:** In `_find_chain_contacts()`, after classifying with `_classify_contact_type()`, only append if `_contact_within_limit(contact_type, distance)`.  
- **total_contacts:** `len(asu_contact_list)` (only contacts within limits).  
- **average_distance, min_distance, max_distance:** From the filtered contact distances.  
- **contact_type_distribution:** Counts per `contact_type`.

### 3.2 Contact density

**Reader-friendly formula**  
*Contact density* = (number of contacts) / (total number of atoms in all chains). It measures how many contacts there are per atom in the structure.

**Mathematical definition:**
\text{contactdensity} = \frac{\text{number of contacts}}{\text{total number of atoms in all chains}}

**Implementation (simplified):** total_atoms = sum of atom counts over all chains. Contact_density = len(contacts) / total_atoms. In `_calculate_contact_density()`.

### 3.3 Crystal contact estimates (simplified)

**Reader-friendly formula**  
*Surface atoms* are those with fewer than 12 neighbors within 5.0 Å (low coordination). The number of *potential crystal contacts* is set equal to the number of surface atoms. *Estimated contact area* is surface_atoms × 4 Å². *Packing efficiency* is (molecular volume) / (cell volume), with molecular volume ≈ atom_count × 20 Å³.

**Implementation (simplified):** For each atom, count how many other atoms are within 5.0 Å. If count < 12, it is a surface atom. surface_atoms = number of such atoms. potential_crystal_contacts = surface_atoms. estimated_contact_area = surface_atoms × 4.0 (Å²). molecular_volume = atom_count × 20.0 (Å³); packing_efficiency = molecular_volume / cell_volume.

- **surface_atoms:** Atoms with < 12 neighbors within 5.0 Å (neighbor count over all atoms).  
- **potential_crystal_contacts:** Set equal to number of surface atoms.  
- **estimated_contact_area:** `surface_atoms * 4.0` (Å²).  
- **packing_efficiency:**  
  - Cell volume from unit cell (same formula as in packing_metrics).  
  - Molecular volume ≈ `atom_count * 20.0` (Å³).  
  - `packing_efficiency = molecular_volume / cell_volume`.

### 3.4 Unit cell volume (contact_analyzer)

- Same parallelepiped formula as in §1.1; implemented in `_calculate_cell_volume()`.

**Derivation / references**

- The ASU contact criteria mirror the interface criteria in §2 with a slightly smaller candidate distance (4.5 Å) to focus on proximate contacts; this is similar in spirit to contact definitions in crystallographic tools such as PISA.
- The **crystal contact estimates** are heuristic:
  - “Surface atoms” are identified via reduced neighbor counts (cf. solvent-exposed residues defined by low coordination in packing analyses).
  - Per-atom area (~~4 Å²) and per-atom volume (~~20 Å³) reflect average SASA and volume for protein atoms in crystal structures.

---

## 4. channel_analyzer.py

**Function selection / usage**

- Channel and void metrics are computed in `channel_analyzer.py`.
- You select them by:
  - Running `channel_analyzer.py` directly on a PDB/CIF (with optional arguments for probe radius, sampling counts, etc.), or
  - Letting `crystal_packing_analyzer.py` call the corresponding `analyze_channels_for_structure` function.
- Inside `channel_analyzer.py`, the main pipeline:
  - Builds coordinate arrays and bounding boxes.
  - Calls `_estimate_void_fraction_and_volume()` → §4.1.
  - Calls `_analyze_density_variations()` → §4.2.
  - Calls `_detect_bottlenecks()` → §4.3.
- **Probe radius:** Default 1.4 Å (used for “accessible” space).

**Libraries / functions used (channel_analyzer)**

- **BioPython:** `PDBParser.get_structure()`, `structure.get_atoms()`; `atom.coord` collected into a NumPy array.
- **NumPy:** `np.array([atom.coord for atom in atoms])`, `np.min(coords, axis=0)`, `np.max(coords, axis=0)`; `np.random.uniform` for sampling points; `np.linalg.norm(atom_coords - point, axis=1)`, `np.min(distances)` in `_is_point_accessible`; `np.prod(max_coords - min_coords)` for box volume; `np.random.choice`, `np.sum(distances <= 8.0)` for density variations; `np.std`, `np.mean` for CV and connectivity; `np.sort(distances)` for bottlenecks; `np.mean`, `np.min` on bottleneck radii.

### 4.1 Void fraction and accessible volume

**Reader-friendly formula**  
*Void fraction* is estimated by Monte Carlo sampling: draw many random points in a box around the structure (atom extent ± 5 Å). A point is *accessible* if it is farther than (probe_radius + atom_radius) from every atom. Void fraction = (number of accessible points) / (number of points). *Accessible volume* = (box volume) × void fraction. Probe radius 1.4 Å approximates a water molecule.

**Mathematical definition**

A point is **accessible** if for every atom, \left\mathbf{p}-\mathbf{x}*a\right > r*{\mathrm{probe}} + r_{\mathrm{atom}} (r_{\mathrm{atom}} = 1.7\text{Å}). Then:
\text{voidfraction} = \frac{\text{number of accessible sampled points}}{\text{number of sampled points}}

V_{\mathrm{accessible}} = V_{\mathrm{box}} \times \text{voidfraction}

- Bounding box: min/max of atom coordinates ± 5 Å. Up to 10,000 random points sampled uniformly in this box.

**Implementation (simplified):** Build a box from atom coordinates (min/max ± 5 Å). Sample many random points (e.g. up to 10,000) uniformly in the box. A point is accessible if, for every atom, distance(point, atom) > probe_radius + atom_radius (atom_radius = 1.7 Å). void_fraction = (number of accessible points) / (number of points). accessible_volume = box_volume × void_fraction.

### 4.2 Density variations and connectivity

**Reader-friendly formula**  
For a sample of atoms, count how many neighbors each has within 8 Å. *Density variation score* is the standard deviation of these counts (spread in local density). *Connectivity score* is the coefficient of variation (std/mean) of the counts, capped at 1 by dividing by 2 and taking min(1, cv/2). *Channel count* is a heuristic: max(1, int(connectivity_score × 5)).

**Mathematical definition**

\mathrm{cv} = \frac{\sigma}{\mu},\quad \text{connectivityscore} = \min\left(1, \frac{\mathrm{cv}}{2}\right)

**Implementation (simplified):** For each of up to 100 random atoms, count how many other atoms are within 8 Å → list of counts. density_variation_score = std(counts). cv = std(counts) / mean(counts); connectivity_score = min(1, cv/2). channel_count = max(1, int(connectivity_score × 5)) (heuristic).

### 4.3 Bottlenecks

**Reader-friendly formula**  
For a sample of atoms, the “bottleneck” at each is the *nearest-neighbor distance* to another atom (second smallest distance, excluding self). If that distance is < 4 Å, a bottleneck is recorded with radius = half that distance. Only bottlenecks with radius ≥ 5 Å are kept. *Bottleneck count*, *average* and *minimum bottleneck radius* are computed from this filtered set.

**Implementation (simplified):** For each of up to 50 atoms, compute distances to all other atoms; take the second smallest (nearest neighbor). If that distance < 4.0 Å, record a bottleneck with radius = distance/2. Keep only bottlenecks with radius ≥ 5 Å. bottleneck_count = number kept; average/minimum_bottleneck_radius = mean/min of those radii.

- Sample up to 50 atoms; for each, take second smallest distance to other atoms (nearest neighbor, excluding self).  
- If that distance < 4.0 Å, record a “bottleneck” with position = atom center, radius = nearest_distance/2.  
- **Minimum channel size:** Only bottlenecks with radius ≥ 5 Å are kept; smaller ones are discarded.  
- **bottleneck_count:** Number of such points (after filtering).  
- **average_bottleneck_radius, minimum_bottleneck_radius:** From the filtered list of radii.

**Derivation / references**

- The sampling-based void fraction is a Monte Carlo approximation to accessible volume, conceptually similar to probe-based cavity finders (e.g. void and channel analyses in tools like CAVER or MOLE).
- The probe radius of 1.4 Å matches the standard water-molecule probe used in SASA calculations (e.g. Lee & Richards, *J. Mol. Biol.* **55**, 379–400, 1971). DOI: [10.1016/0022-2836(71)90324-x](https://doi.org/10.1016/0022-2836(71)90324-x)
- Bottleneck detection via nearest-neighbor distances is a simplified proxy for minimum pore radii along putative channels; the 5 Å minimum radius is chosen so that only biologically relevant or solvent-accessible channels are retained.

---

## 5. graph_analyzer.py

**Function selection / usage**

- Graph-based metrics are computed in `graph_analyzer.py`.
- You select them by:
  - Running `graph_analyzer.py` on one or more structures, or
  - Using `crystal_packing_analyzer.py`, which calls a function such as `analyze_graphs_for_structure` for each input.
- Inside `graph_analyzer.py`, the main workflow:
  - Builds residue-level and chain-level graphs from contact information (using `contact_threshold`).
  - Calls helper functions that wrap NetworkX routines:
    - Graph density, degrees, paths, clustering → §§5.1–5.3.
    - Centralities, communities, modularity → §§5.4–5.5.
    - Random-graph comparison and small-world coefficient → §5.6.
- Residues are nodes; two residues are connected if minimum atom–atom distance ≤ **contact_threshold** (default 5.0 Å).  
- Chain graph: chains as nodes; edge if they have at least one residue–residue contact.

**Libraries / functions used (graph_analyzer)**

- **BioPython:** `PDBParser.get_structure()`, `structure.get_chains()`, `structure.get_residues()`; `chain.get_residues()`, `residue.get_atoms()`; `atom.coord`. Residue–residue distance = minimum over atom pairs of `np.linalg.norm(a1.coord - a2.coord)`.
- **NumPy:** `np.linalg.norm` for distances.
- **NetworkX:** `nx.Graph()`, `graph.add_node()`, `graph.add_edge()`; `nx.density(graph)`; `nx.is_connected(graph)`; `nx.average_shortest_path_length(graph)`, `nx.diameter(graph)`; `nx.average_clustering(graph)`; `nx.betweenness_centrality(graph)`, `nx.closeness_centrality(graph)`, `nx.degree_centrality(graph)`; `nx.connected_components(graph)`; `nx.algorithms.community.greedy_modularity_communities(residue_graph)`, `nx.algorithms.community.modularity(residue_graph, communities)`; `nx.erdos_renyi_graph(n, density)` for random graph.

### 5.1 Residue and chain graphs

**Reader-friendly formula**  
The *residue graph* has one node per residue; an edge exists if the minimum distance between any atom of one residue and any atom of the other is ≤ contact_threshold (5 Å). The *chain graph* has one node per chain; an edge exists if at least one residue–residue contact exists between the two chains. *Graph density* = 2×|edges| / (n(n−1)). *Average degree* = mean of each node’s number of edges.

**Mathematical definition**

\text{density} = \frac{2|E|}{n(n-1)},\quad \text{averagedegree} = \frac{1}{n}\sum_i \deg(i)

**Implementation (simplified):** Residue graph: one node per residue; edge between two residues if the minimum atom–atom distance between them is ≤ contact_threshold (5 Å). Chain graph: one node per chain; edge if at least one residue–residue contact between those chains. node_count = number of nodes; edge_count = number of edges. density = 2×edge_count / (n×(n−1)). average_degree = sum of degrees / n. Script: `nx.density(graph)`, etc.

### 5.2 Shortest paths and diameter

**Reader-friendly formula**  
*Average path length* is the mean number of steps along the shortest path over all pairs of nodes (on the largest connected component if the graph is disconnected). *Diameter* is the maximum shortest-path length over those pairs.

**Implementation (simplified):** For each pair of nodes, shortest path length = number of edges along the shortest path. average_path_length = mean of these (over pairs in the largest component if disconnected). diameter = maximum of these. Script: `nx.average_shortest_path_length()`, `nx.diameter()`.

- **average_path_length:** Mean shortest-path length over all node pairs (or on largest connected component if graph disconnected).  
- **diameter:** Maximum shortest-path length (or on largest component).  
- **Script:** `nx.average_shortest_path_length()`, `nx.diameter()`.

### 5.3 Clustering coefficient

**Reader-friendly formula**  
For each node, the *local clustering coefficient* is (number of triangles that include that node) / (number of possible triangles among its neighbors). The *average clustering coefficient* is the mean over all nodes.

**Mathematical definition**

C = \frac{1}{n}\sum_i C_i,\quad C_i = \frac{\text{number of triangles involving node } i}{\binom{\deg(i)}{2}}

**Implementation (simplified):** For each node i, C_i = (number of triangles involving i) / (number of pairs of i's neighbors). Average clustering = mean of C_i over all nodes. Script: `nx.average_clustering(graph)`.

### 5.4 Centrality (residue/chain graph)

**Reader-friendly formula**  
*Betweenness:* fraction of shortest paths between pairs of nodes that pass through this node. *Closeness:* inverse of the mean distance from this node to all others. *Degree centrality:* node degree / (n−1).

**Implementation (simplified):** For each node: betweenness = (number of shortest paths that pass through this node) / (total number of shortest paths); closeness = 1 / (mean distance from this node to others); degree_centrality = degree / (n−1). Script reports max and mean over nodes: `nx.betweenness_centrality()`, `nx.closeness_centrality()`, `nx.degree_centrality()`.

- **betweenness_centrality:** Fraction of shortest paths passing through each node.  
- **closeness_centrality:** Inverse of mean distance from node to others.  
- **degree_centrality:** degree(node) / (n−1).  
- **Script:** `nx.betweenness_centrality()`, `nx.closeness_centrality()`, `nx.degree_centrality()`; reported as max and mean over nodes.

### 5.5 Modularity and communities

**Reader-friendly formula**  
*Communities* are found by greedy modularity maximization (NetworkX). *Modularity* measures how much the edge density within communities exceeds that expected in a random graph with the same degrees.

**Implementation (simplified):** Run greedy modularity maximization to partition nodes into communities. modularity = measure of how clustered the edges are within communities vs random. community_count = number of communities. Script: `nx.algorithms.community.greedy_modularity_communities(residue_graph)`, `nx.algorithms.community.modularity(residue_graph, communities)`.

- Communities: `nx.algorithms.community.greedy_modularity_communities(residue_graph)`.  
- **modularity:** `nx.algorithms.community.modularity(residue_graph, communities)`.  
- **community_count:** Number of communities.

### 5.6 Small-world coefficient

**Reader-friendly formula**  
A *small-world* graph has high clustering (like a regular lattice) but short average path length (like a random graph). The script builds an Erdős–Rényi random graph with the same number of nodes and same edge density, then computes: *clustering_ratio* = (clustering of residue graph) / (clustering of random graph); *path_length_ratio* = (avg path of residue graph) / (avg path of random graph). *Small-world coefficient* = clustering_ratio / path_length_ratio (0 if path_ratio = 0).

**Mathematical definition**

\text{smallworldcoefficient} = \frac{\text{clusteringratio}}{\text{pathlengthratio}} \quad (0 \text{ if } \text{pathlengthratio}=0).

**Implementation (simplified):** Build a random graph (Erdős–Rényi) with the same number of nodes and same edge density. clustering_ratio = (clustering of residue graph) / (clustering of random graph). path_length_ratio = (average path length of residue graph) / (average path length of random graph). small_world_coefficient = clustering_ratio / path_length_ratio (0 if path_length_ratio = 0).

- Build Erdős–Rényi random graph with same n and same edge density.  
- **clustering_ratio:** `average_clustering(residue_graph) / average_clustering(random_graph)`.  
- **path_length_ratio:** `average_shortest_path_length(residue_graph) / average_shortest_path_length(random_graph)`.

**Derivation / references**

- Graph-theoretic measures follow standard NetworkX and network-science definitions:
  - M. E. J. Newman, *Networks: An Introduction* (Oxford University Press, 2010).
- Small-world analysis mirrors the Watts–Strogatz framework, comparing real-graph clustering/path lengths to an Erdős–Rényi random graph:
  - D. J. Watts & S. H. Strogatz, “Collective dynamics of ‘small-world’ networks”, *Nature* **393**, 440–442 (1998). DOI: `10.1038/30918`.
- Community detection and modularity use greedy modularity maximization implemented in NetworkX; modularity is as defined by Newman–Girvan.

---

## 6. comparative_analyzer.py

**Function selection / usage**

- Comparative metrics are computed in `comparative_analyzer.py` and operate on a combined table (DataFrame) of per-structure metrics.
- You select them by:
  - Running `comparative_analyzer.py` directly and pointing it to a metrics CSV/JSON, or
  - Running the “comparative” step in `crystal_packing_analyzer.py`, which:
    - Aggregates outputs from `packing_metrics`, `interface_analyzer`, `contact_analyzer`, `channel_analyzer`, and `graph_analyzer`.
    - Passes the resulting DataFrame into summary, correlation, clustering, PCA, and outlier routines.

Operates on a **DataFrame** of metrics (one row per structure) extracted from packing_metrics, interface_analysis, contact_analysis, and graph_analysis.

**Libraries / functions used (comparative_analyzer)**

- **pandas:** `pd.DataFrame(metrics_list)`, `df.select_dtypes(include=[np.number])` to get numeric columns.
- **NumPy:** `np.mean`, `np.std`, `np.min`, `np.max`, `np.median` on column arrays; `np.corrcoef(col1, col2)[0, 1]` for Pearson correlation; `np.isnan`; `np.cumsum` for cumulative variance; `np.abs`, `np.where(z_scores > 2.5)[0]` for outliers.
- **scikit-learn:** `StandardScaler().fit_transform(numeric_df)`; `KMeans(n_clusters=k, random_state=42, n_init=10).fit()`; `PCA(n_components=...).fit_transform()`; `pca.components`*, `pca.explained_variance_ratio`*.

### 6.1 Summary statistics (per metric column)

**Reader-friendly formula**  
For each numeric metric column, the script computes the **mean**, **standard deviation**, **min**, **max**, and **median** over the structures (rows). These are the usual descriptive statistics.

**Mathematical definition**

For each numeric column x: \text{mean} = \bar{x},\quad \text{std} = \sigma_x,\quad \text{min},\ \text{max},\quad \text{median}

**Implementation (simplified):** For each numeric column (one metric), take the values across all rows (structures). mean = average of values; std = standard deviation; min = smallest; max = largest; median = middle value. Script: `np.mean`, `np.std`, `np.min`, `np.max`, `np.median` on column values.

### 6.2 Correlations

**Reader-friendly formula**  
*Pairwise Pearson correlation* between every pair of numeric metric columns: for each pair (i, j), the script computes the linear correlation coefficient between the two columns and stores it (e.g. as "metric1_vs_metric2": value).

**Implementation (simplified):** For each pair of numeric columns (i < j), compute the Pearson correlation between the two lists of values. Store as e.g. "metric1_vs_metric2": correlation. Script: `np.corrcoef(col1, col2)[0, 1]`.

- **Definition:** Pairwise Pearson correlation between metric columns.  
- **Script:** For each pair of columns (i < j), `np.corrcoef(col1, col2)[0,1]`; stored as e.g. `"metric1_vs_metric2": corr`.

### 6.3 Clustering (K-means)

**Reader-friendly formula**  
Data are standardized (zero mean, unit variance per column) with scikit-learn’s `StandardScaler`. The script tries K = 2 up to min(5, n_structures−1) and picks the K with the smallest K-means inertia (within-cluster sum of squared distances). It then reports optimal K, cluster label per structure, and inertia. Implemented with `KMeans(n_clusters=k, random_state=42, n_init=10)`.

- Data standardized with `sklearn.preprocessing.StandardScaler`.  
- K from 2 to min(5, n_structures−1); choose K with smallest K-means inertia.  
- **Script:** `KMeans(n_clusters=k, random_state=42, n_init=10)`; report `optimal_clusters`, `cluster_assignments`, `inertia`.

### 6.4 PCA

**Reader-friendly formula**  
The same standardized data are projected onto principal components. The script uses as many components as min(n_features, n_structures). It reports *explained variance ratio* per component, *cumulative variance* (cumsum of those), and the first two rows of `pca.components`_ as PC1 and PC2 loadings.

**Implementation (simplified):** Use the same standardized numeric table. Fit PCA with n_components = min(n_features, n_structures). explained_variance_ratio = fraction of total variance per component; cumulative_variance = cumsum of those. PC1 and PC2 loadings = first two rows of pca.components_. Script: `sklearn.decomposition.PCA`, `np.cumsum`.

- Same standardized data; PCA with `n_components = min(n_features, n_structures)`.  
- **explained_variance_ratio:** Per-component.  
- **cumulative_variance:** `np.cumsum(explained_variance_ratio)`.  
- **pc1_loadings, pc2_loadings:** First two rows of `pca.components`_.

### 6.5 Outliers

**Reader-friendly formula**  
For each metric column, a *z-score* is computed for each value: z = (x − mean) / std. Values with |z| > 2.5 are flagged as outliers. The script reports structure id, metric name, value, and z-score for each such case.

**Mathematical definition:** For each metric, z = (x - \bar{x})/\sigma_x. Values with |z| > 2.5 are flagged as outliers.

**Implementation (simplified):** For each numeric column, compute mean and std over the values. For each value x, z = (x − mean) / std. Flag as outlier if |z| > 2.5. Report structure id, metric name, value, z-score for each flagged value. Script: `z_scores = np.abs((values - np.mean(values)) / np.std(values))`; `np.where(z_scores > 2.5)[0]`.

- **Definition:** For each metric, values with |z-score| > 2.5, where z = (x - mean)/std.  
- **Script:** `z_scores = np.abs((values - np.mean(values)) / np.std(values))`; `outlier_indices = np.where(z_scores > 2.5)[0]`; report structure_id, metric, value, z_score.

**Derivation / references**

- Summary statistics, Pearson correlations, K-means clustering, and PCA follow standard definitions as implemented in NumPy/Scikit-learn.
- The outlier threshold |z| > 2.5 is a common heuristic in exploratory data analysis, roughly corresponding to tail probabilities ≲ 1% for Gaussian-like distributions.

---

## 7. crystal_packing_analyzer.py

This script orchestrates the above: it runs packing_metrics, interface_analyzer, contact_analyzer, channel_analyzer, and graph_analyzer on each structure, then (optionally) runs comparative_analyzer on the combined results. It does not define new formulas; all metrics are as in §§1–6.

**Function selection / usage**

- When you run `crystal_packing_analyzer.py` from the command line, you choose:
  - Which input structures (files or directories) to analyze.
  - Optional flags that enable/disable particular stages (packing, interfaces, contacts, channels, graphs, comparative summary).
- Internally, it:
  - Iterates over structures and calls the appropriate per-structure analyzer functions from §§1–5.
  - After all structures are processed, assembles a combined table and passes it to `comparative_analyzer.py` if requested.

---

## 8. structure_phylogeny.py

This script builds structure-based phylogenies from pairwise structural comparisons provided as an RMSD matrix (precomputed by your Coot LSQ/SSM pipeline, or optionally computed via `TM-align`).

It supports two related “distance matrix” modes for tree inference:

1. Default: use RMSD directly as the distance.
2. DALI-like “per-query pseudo Z-score” mode: convert RMSD to similarity, then to per-query z-scores, then to distances.

### 8.1 Default mode: RMSD distance matrix

**Mathematical definition**

Let `D` be the RMSD matrix with entries `D_ij >= 0`. The distance matrix for phylogeny construction is:

`d_ij = D_ij` for `i != j`, and `d_ii = 0`.

**LaTeX version**

d_{ij}=
\begin{cases}
D_{ij}, & i\neq j
0, & i=j
\end{cases}


**Tree construction**

The script feeds this distance matrix to neighbor-joining (NJ) when available (otherwise it uses an internal UPGMA fallback).

---

### 8.2 DALI-like per-query pseudo Z-scores (RMSD -> distance)

This mode is activated by `structure_phylogeny.py --pseudo-z per_query`.

Let `D_ij` be the RMSD distance between structures `i` and `j`. Define a similarity matrix `S`:

#### Step 1: RMSD -> similarity (`--similarity`)

1. `--similarity neg_rmsd`:

`S_ij = -D_ij`

1. `--similarity exp` (default):

`S_ij = exp(-D_ij / tau)`

where `tau` is provided by `--tau`.

**LaTeX version**

S_{ij}=
\begin{cases}
-D_{ij} & \text{(negrmsd)}
\exp\left(-D_{ij}/\tau\right) & \text{(exp)}
\end{cases}


#### Step 2: per-query z-score standardization

For each query `i`, compute the mean and population standard deviation over all `j != i`:

`mu_i = mean_{j != i}(S_ij)`

`sigma_i = sqrt( mean_{j != i}( (S_ij - mu_i)^2 ) )`

Then define the (asymmetric) query-centric z-score:

`z_ij = (S_ij - mu_i) / sigma_i`

(If `sigma_i` is numerically ~0, the code leaves z-scores at 0 for that row.)

**LaTeX version**

\mu_i=\frac{1}{m_i}\sum_{j\neq i} S_{ij},\quad
\sigma_i=\sqrt{\frac{1}{m_i}\sum_{j\neq i}\left(S_{ij}-\mu_i\right)^2},\quad
z_{ij}=\frac{S_{ij}-\mu_i}{\sigma_i}.


#### Step 3: symmetrize the z-score matrix

The code symmetrizes z-scores by averaging:

`z^sym_ij = (z_ij + z_ji) / 2`

so the downstream distance matrix is symmetric.

**LaTeX version**

z^{\mathrm{sym}}*{ij}=\frac{z*{ij}+z_{ji}}{2}


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

**LaTeX version**

z^{+}*{ij}=\max\left(0, z^{\mathrm{sym}}*{ij}\right),\qquad
d_{ij}=
\begin{cases}
\dfrac{1}{1+z^{+}*{ij}} & \text{(inv)}
\exp\left(-z^{+}*{ij}/z_{\mathrm{scale}}\right) & \text{(exp)}
\dfrac{z_{\max}-z^{+}*{ij}}{z*{\max}} & \text{(maxminus)}
\end{cases}


#### Query-centric ranking (optional)

When `--ranking-csv` is provided, the script also writes a ranking table per structure `i` computed from the *asymmetric* z-scores `z_ij` (before symmetrization):

`avg_z_i = mean_{j != i}( z_ij )`

`max_z_i = max_{j != i}( z_ij )`

and ranks by descending `(avg_z_i, max_z_i)`.

**LaTeX version**

\bar{z}*i=\frac{1}{m_i}\sum*{j\neq i} z_{ij},\qquad
z^{\max}*i=\max*{j\neq i} z_{ij}.


### 8.3 Derivation / references (why these formulas)

- The DALI query-centric Z-score framework is rooted in Holm & Sander’s distance-matrix alignment framework:
  - Holm, L.; Sander, C. *Protein structure comparison by alignment of distance matrices.* *J. Mol. Biol.* (1993). DOI: `10.1006/jmbi.1993.1489`
- The `exp(-D/tau)` similarity is a common distance-kernel/RBF-style choice:
  - Micchelli, C. *Interpolation of scattered data: Distance matrices and conditionally positive definite functions.* *Constructive Approx.* (1986). DOI: `10.1007/BF01893414`
- Neighbor-joining is the primary tree-building method:
  - Saitou, N.; Nei, M. *The neighbor-joining method: a new method for reconstructing phylogenetic trees.* *Mol. Biol. Evol.* (1987). DOI: `10.1093/oxfordjournals.molbev.a040454`
- The per-query z-score standardization uses conventional z-score definitions. The final `z^+ -> distance` mapping is a practical heuristic to produce non-negative distances usable by NJ/UPGMA.

---

## 9. dali_score.py

**Function selection / usage**

- Compute a Dali-like structural similarity score from two structures using residue equivalences.
- **Pairwise:** `python dali_score.py pdb_a pdb_b [-a alignment_file] [--chain-a ID] [--chain-b ID] [--dalilite-path DIR] [-o output.csv]`
- **All-vs-all:** `python dali_score.py --all-vs-all dir_or_file [dir_or_file ...] [--filter PATTERN] [--output-tree FILE] [--output-plot FILE] [--output-matrix FILE] [--output-ranking FILE]`
- Alignment sources (in order): (1) **DaliLite** (if `--dalilite-path` or `DALILITE_HOME` is set by user), (2) `--alignment` file, (3) biotite structural alignment (if installed), (4) sequence-order matching (same residue numbering). When DaliLite is used, its canonical Z-score and RMSD are reported.

**Libraries / functions used (dali_score)**

- **BioPython:** `PDBParser.get_structure()`, structure iteration for CA atoms; `atom.coord`.
- **NumPy:** `np.linalg.norm`, `np.exp`, `np.atleast_1d`; distance matrices for pairwise Cα–Cα distances.
- **biotite (optional):** `PDBFile.read()`, `atom_array[atom_array.atom_name == 'CA']`, `superimpose_structural_homologs()` returning `fixed_indices`, `mobile_indices` for residue equivalences.

### 9.1 Raw Dali score

**Reader-friendly formula**  
The *raw Dali score* measures structural similarity by comparing intramolecular Cα–Cα distances between two structures. For each pair of equivalent residues (i, j) in the aligned core, compute the distance between residues i and j within structure A and within structure B. The score φ(i,j) rewards similar distances and penalizes deviations; an envelope downweights long-range pairs. The total score is the sum over all residue pairs in the core.

**Mathematical definition**

S(A,B) = ∑_{i ∈ core} ∑_{j ∈ core} φ(i,j)

where

φ(i,j) = (θ − diff(i,j)) × exp(−(d*_{ij}/R₀)²)

θ = 0.2, R₀ = 20 Å, diff(i,j) = |d_{ij}^A − d_{ij}^B| / d*_{ij}, d*_{ij} = (d_{ij}^A + d_{ij}^B)/2.

Only terms with φ > 0 contribute. d_{ij}^A and d_{ij}^B are intramolecular Cα–Cα distances in structures A and B.

**Implementation (simplified):** For each pair (i,j) in the core, compute d_{ij}^A and d_{ij}^B from CA coordinates. d* = (d_A + d_B)/2; diff = |d_A − d_B|/d*; φ = max(0, (0.2 − diff) × exp(−(d*/20)²)). Sum all φ. Script: `compute_dali_score()`, `_residue_pair_score()`.

### 9.2 Z-score

**Reader-friendly formula**  
The *Z-score* normalizes the raw score for protein length so that scores are comparable across pairs of different sizes. It uses empirically fitted mean m(L) and standard deviation σ(L) for random structure pairs of effective length L = √(L_A × L_B).

**Mathematical definition**

Z(A,B) = (S(A,B) − m(L)) / σ(L)

L = √(L_A × L_B)

m(L) ≈ 7.95 + 0.71L − 2.59×10⁻⁴ L² − 1.92×10⁻⁶ L³   for L ≤ 400

m(L) = m(400) + (L − 400)   for L > 400

σ(L) = 0.5 × m(L)

**Implementation (simplified):** L = sqrt(L_A * L_B); m = polynomial or linear extrapolation; σ = 0.5*m; Z = (S − m) / σ. Script: `compute_z_score()`, `_mean_score_fit()`, `_std_score_fit()`.

### 9.3 DaliLite integration

When **DaliLite** is available (user specifies its installation directory), it is used as the primary alignment source. In this repository, **`dali_score.py` does not rely on bare `dali.pl --pdbfile1/--pdbfile2`** (which often yields empty internal data). It uses the same pattern as DaliLite v5 batch workflows:

1. **Short working directory:** comparisons run under a short temporary directory (typically a randomly named subdirectory of the system temp folder; DaliLite Fortran binaries enforce a maximum path length of about 80 characters).
2. **Staging:** PDB/mmCIF inputs may be copied or rewritten; PDB files without a `HEADER` record can get a dummy `HEADER` so `mkdssp` treats them as PDB, not mmCIF (`_stage_pdb_for_dalilite_import`).
3. **`import.pl`** builds `DAT/*.dat` per structure: `import.pl --pdbfile … --pdbid xxxx --dat DAT/`.
4. **`dali.pl`** pairwise comparison: `dali.pl --cd1 xxxxX --cd2 yyyyY --dat1 DAT` (with `--outfmt` for summary, equivalences, transrot as needed).

The script parses **`dali.pl` text output** for:

- **Z-score** (canonical Dali Z-score, used in preference to the empirical fit)
- **RMSD** of structurally equivalent Cα atoms where reported
- **Residue equivalences** from the structural equivalences block

**Configuration:** set `--dalilite-path` or **`DALILITE_HOME`** to the DaliLite installation root (directory containing `bin/dali.pl`).

**mkdssp (required for `import.pl`):** DaliLite calls **`mkdssp`** using the path in **`bin/mpidali.pm`** (`$DSSP_EXE`; upstream defaults often use `/usr/local/bin/mkdssp`). That file must exist and be executable—otherwise `DAT/*.dat` can be empty and comparisons fail silently. Fix by installing a compatible **`mkdssp`**, symlinking it to the path in `mpidali.pm`, or editing `$DSSP_EXE`. Perl modules under `bin/` (e.g. `FSSP.pm`) must be readable by the user. Optionally set **`MKDSSP`** to the full path of `mkdssp` so this script’s diagnostics match your install.

**Manual all-vs-all** (outside `dali_score.py`): batch-`import.pl` for each model into a shared `DAT/`, build `query.list` of five-character ids (`xxxxX` per `DAT/*.dat`), then `dali.pl --matrix --query query.list --dat1 DAT` from a short working directory path.

### 9.3.1 dalilite_superpose_scores.py

**`dalilite_superpose_scores.py`** runs the DaliLite pairwise path via `dali_score._dalilite_pair_via_dat`, then applies **translation–rotation** from `--outfmt transrot` to write a superposed target PDB (BioPython). It supports **all-vs-all**, CSV scores, Z-matrix, and Newick tree output.

Requirements and environment variables match **`dali_score.py`** (DaliLite path, mkdssp as above). Use **`--fallback-biotite`** when DaliLite reports no significant hit (often below Dali’s usual reporting threshold, Z ~ 2) but a structural alignment is still desired.

### 9.3.2 run_all_superpositions.py (Coot LSQ batch)

**`run_all_superpositions.py`** is a small driver that loops over **condition** subdirectories (e.g. `condition_a`) and **tags** (same idea as `--filter=tag_a` in the Coot LSQ/SSM scripts), calling **`superimpose_coot_LSQ.py`** with **`--pattern`**, **`--divider=LSQ_`**, per-tag reference globs, and **`CONDITION_PREFIX`**-derived filename tokens. Defaults such as **`REFERENCE_SUBDIR`** and **`REFERENCE_STEM`** are documented in the script and in the main **`README.md`** (superposition section).

### 9.4 All-vs-all and tree output

When `--all-vs-all` is used with one or more directories or PDB/CIF files, the script compares every pair of structures, collects Z-scores, and optionally generates:

- **Pairwise table** (`-o`): CSV with pdb_a, pdb_b, raw_score, z_score, n_core
- **Ranking** (`--output-ranking`, default `dali_ranking.csv`): structures ranked by average and max Z-score
- **Z-score matrix** (`--output-matrix`): symmetric matrix CSV for downstream use
- **Newick tree** (`--output-tree`, default `dali_tree.nwk`): phylogenetic tree from Z-score–derived distances (neighbor-joining or UPGMA)
- **Dendrogram plot** (`--output-plot`): PNG/SVG tree visualization (requires ete3 or biopython+matplotlib)

**Z-score → distance transform** (`--transform`): Higher Z = more similar. For trees, Z is converted to distance: `inv` (d = 1/(1+Z)), `maxminus` (normalized), or `exp` (d = exp(-Z/scale)). Tree building uses scikit-bio or Biopython for neighbor-joining when available, else a pure-Python UPGMA fallback. See Saitou & Nei (1987) for NJ.

**Filtering:** `--filter` limits which files are included: plain text matches as a substring on the basename; patterns with `*`, `?`, or `[` use `fnmatch` on basename or stem (see `_collect_pdb_files` in `dali_score.py`).

### 9.5 Residue equivalences

**Alignment file format (optional):** TSV or CSV with either:
- Two columns: `resnum_A`, `resnum_B` (for single chain per structure)
- Four columns: `chain_A`, `resnum_A`, `chain_B`, `resnum_B`

Lines starting with `#` are comments. Header row is auto-skipped.

**Automatic alignment:** When no alignment file is provided:
1. **biotite** (if installed): `superimpose_structural_homologs()` returns fixed/mobile indices of equivalent CA atoms; mapped to (chain, resseq, icode).
2. **sequence-order:** For structures sharing residue numbering (e.g. same chain IDs and residue numbers), use 1:1 correspondence by (chain, resseq).

**Derivation / references**

- Holm, L.; Sander, C. *Protein structure comparison by alignment of distance matrices.* *J. Mol. Biol.* **233**, 123–138 (1993). DOI: [10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489) — Original Dali method: raw score φ(i,j), Z-score normalization, empirical m(L) and σ(L).
- Holm, L. *Dali and the persistence of protein shape.* *Protein Sci.* **29**, 128–140 (2020; first published online 2019). DOI: [10.1002/pro.3749](https://doi.org/10.1002/pro.3749) — DaliLite standalone implementation.
- The Dali score uses intramolecular Cα–Cα distances only; superposition is not required. Equivalences can come from DaliLite, Coot/SSM (export alignment if available), biotite, or any structural alignment tool.

---

## Reference: Default parameters


| Script / class                       | Parameter                                             | Default                         |
| ------------------------------------ | ----------------------------------------------------- | ------------------------------- |
| interface_analyzer                   | contact_distance                                      | 5.0 Å                           |
| contact_analyzer                     | contact_distance                                      | 4.5 Å                           |
| channel_analyzer                     | probe_radius                                          | 1.4 Å                           |
| channel_analyzer                     | minimum channel radius (bottlenecks < this discarded) | 5.0 Å                           |
| graph_analyzer                       | contact_threshold                                     | 5.0 Å                           |
| packing_metrics                      | protein density (for solvent)                         | 0.81 Da/Å³                      |
| interface_analyzer                   | BSA (fallback when SASA unavailable)                  | 15 Å² per contact atom          |
| contact_analyzer                     | surface neighbor cutoff                               | 5.0 Å; < 12 neighbors → surface |
| contact_analyzer                     | estimated area per contact                            | 4 Å²                            |
| contact_analyzer                     | avg atomic volume (packing)                           | 20 Å³                           |
| comparative_analyzer                 | outlier z-score threshold                             | 2.5                             |
| interface_analyzer, contact_analyzer | H-bond max distance                                   | 3.5 Å                           |
| interface_analyzer, contact_analyzer | Electrostatic max distance                            | 5.0 Å                           |
| interface_analyzer, contact_analyzer | Hydrophobic range                                     | 3.5–4.5 Å                       |
| interface_analyzer, contact_analyzer | Van der Waals max distance                            | 5.0 Å                           |

---

## References (selected)

*Bibliographic metadata and DOIs below were checked against the [Crossref](https://www.crossref.org/) REST API (works endpoint). PubMed-indexed articles align with the same DOIs where applicable.*

- **Holm, L.; Sander, C.** (1993). Protein structure comparison by alignment of distance matrices. *J. Mol. Biol.* **233**, 123–138. DOI: [10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489)
- **Holm, L.** (2020). Dali and the persistence of protein shape. *Protein Sci.* **29**, 128–140 (print; first published online 2019). DOI: [10.1002/pro.3749](https://doi.org/10.1002/pro.3749)
- **Saitou, N.; Nei, M.** (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Mol. Biol. Evol.* **4**, 406–425. DOI: [10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)
- **Matthews, B.W.** (1968). Solvent content of protein crystals. *J. Mol. Biol.* **33**, 491–497. DOI: [10.1016/0022-2836(68)90205-2](https://doi.org/10.1016/0022-2836(68)90205-2)
- **Shrake, A.; Rupley, J.A.** (1973). Environment and exposure to solvent of protein atoms. *J. Mol. Biol.* **79**, 351–371. DOI: [10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)


