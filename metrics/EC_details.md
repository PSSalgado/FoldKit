## Electrostatic complementarity (EC) metrics

This document defines the electrostatic complementarity (EC) metrics reported by the
FoldKit lattice analysers. The EC implementation follows the definition of McCoy,
Epa & Colman (1997): EC is the correlation of electrostatic potential on facing
protein surfaces at an interface.

### References

- **McCoy, A. J.; Epa, V. C.; Colman, P. M.** (1997). *Electrostatic complementarity at protein/protein interfaces.* **Journal of Molecular Biology**, 268, 570–584. DOI: 10.1006/jmbi.1997.0987.

### Definitions and notation

Consider an interface between a reference chain *A* and a partner chain *B*.

- Let *P*_A = { *p*_i^A } for *i* = 1, …, *n* be a set of solvent-accessible surface sample points on chain *A*
  that are classified as interfacial.
- Let *P*_B = { *p*_j^B } for *j* = 1, …, *m* be the corresponding set on chain *B*.
- A discrete "facing" pairing is defined by mapping each *p*_i^A to its nearest neighbour in *P*_B:

```
π(i) = arg min_{1 ≤ j ≤ m} ‖ p_i^A − p_j^B ‖
```

This yields paired points (*p*_i^A, *p*_{π(i)}^B) for *i* = 1, …, *n*. The number of pairs is:

```
n_pairs = n
```

### Electrostatic potential at surface points

The electrostatic potential at a point *p* is computed as a Coulomb sum over partial charges:

```
φ(p) = Σ_k  q_k / ( ε_r × max(r_min, ‖ p − r_k ‖) )
```

where:
- *q*_k is the partial charge of atom *k* (in units of elementary charge *e*),
- *r*_k is the coordinate of atom *k*,
- ε_r is an effective relative dielectric constant,
- *r*_min prevents divergence at short distances.

For an interface, the potentials are evaluated in the *partner field* sense:

- φ_A(*i*) = φ_B(*p*_i^A): potential on *A*'s surface point due to charges on *B*,
- φ_B(*i*) = φ_A(*p*_{π(i)}^B): potential on *B*'s facing surface point due to charges on *A*.

### Per-interface EC (correlation score)

Electrostatic complementarity for the *A*–*B* interface is defined as the Pearson correlation:

```
EC(A,B) = r_AB = corr( φ_A, −φ_B )
```

with:

```
corr(x,y) = Σ_{i=1}^n (x_i − x̄)(y_i − ȳ)
            ─────────────────────────────────────────────────────────
            √(Σ_{i=1}^n (x_i − x̄)²) × √(Σ_{i=1}^n (y_i − ȳ)²)
```

Interpretation:
- *r*_AB > 0: potentials tend to be opposite-sign across the interface (electrostatically complementary).
- *r*_AB ≈ 0: weak or no relationship.
- *r*_AB < 0: potentials tend to be same-sign (electrostatically anti-complementary).

### Per-interface EC density

To normalise EC by interface size, an EC density is reported as:

```
ECdensity(A,B) = r_AB / BSA_AB
```

where BSA_AB is the buried surface area of the *A*–*B* interface (Å²).

Units: Å⁻².

### Lattice EC: reference chain against multiple partners

In lattice mode, chain *A* forms interfaces with multiple partner chains *B*_k.
The lattice summary reports two weighted EC scores computed from the set of per-partner
EC correlations { *r*_k }.

Because correlations are not additive, weighted aggregation is performed using Fisher’s
z-transform:

```
z_k = atanh(r_k)
```

For weights *w*_k, the weighted mean is:

```
z̄ = ( Σ_k w_k z_k ) / ( Σ_k w_k )
```

and the aggregated lattice EC is:

```
r_lattice = tanh(z̄)
```

Two weighting schemes are reported:

- **BSA-weighted**: *w*_k = BSA_{A B_k}
- **n_pairs-weighted**: *w*_k = n_{pairs,k}

### Lattice EC density

Lattice EC density is reported by normalising by the reference buried area:

```
ECdensity_lattice = r_lattice / A_buried
```

where:

```
A_buried = SASA_iso(A) − SASA_cluster(A)
```

Here SASA_iso(*A*) is the Shrake–Rupley solvent-accessible surface area of chain *A*
in isolation, and SASA_cluster(*A*) is the sum of per-residue SASA for chain *A* computed
in the full lattice model (all chains present, occluding the probe).

Units: Å⁻².

### Reported fields

Per interface:
- **EC (r)**: *r*_AB
- **n_pairs**: *n*_pairs
- **EC density**: *r*_AB / BSA_AB

Lattice summary:
- **EC (r, BSA-weighted Fisher-z)**: aggregated *r*_lattice with *w*_k = BSA
- **EC (r, n_pairs-weighted Fisher-z)**: aggregated *r*_lattice with *w*_k = *n*_pairs
- **EC density**: aggregated *r*_lattice / *A*_buried
- **EC by partner chain**: (*r*_k, *n*_{pairs,k}) for each partner.

### Operational note: two-phase lattice EC (same physics)

For **`interface_analyser_lattice_ec.py`**, optional **`--phase sasa`** then **`--phase ec`** splits SASA/BSA work from EC while preserving the same metric definitions; use when runtime or scheduler limits favour two jobs over one. Requirements (matching PDB fingerprint and CLI flags) and examples: **README.md** (Metrics, **`interface_analyser_lattice_ec.py`** subsection).

### Related: Caver tunnel “EC” (different geometry)

`metrics/caver_tunnel_analysis.py` reports a **tunnel**-specific analogue of electrostatic complementarity, built on a **Caver 3.0** centreline and a shell of **lining** residues, using the same class of **Coulomb-summed** potential described above for interfaces. The **geometric** construction is **not** the protein–protein, SASA-facing point pairing used in `interface_analyser_asu_ec.py`, `interface_analyser_lattice_ec.py`, and `electrostatic_complementarity.py`. **Do not** treat tunnel and interface *r* values as interchangeable without clear labelling. Full behaviour, options, and outputs are documented in **README.md** (Metrics) and **metrics/metrics_details.md** (Section 1.7).
