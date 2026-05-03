## Electrostatic complementarity (EC) metrics

This document defines the electrostatic complementarity (EC) metrics reported by the
FoldKit lattice analysers. The EC implementation follows the definition of McCoy,
Epa & Colman (1997): EC is the correlation of electrostatic potential on facing
protein surfaces at an interface.

### References

- **McCoy, A. J.; Epa, V. C.; Colman, P. M.** (1997). *Electrostatic complementarity at protein/protein interfaces.* **Journal of Molecular Biology**, 268, 570–584. DOI: 10.1006/jmbi.1997.0987.

### Definitions and notation

Consider an interface between a reference chain \(A\) and a partner chain \(B\).

- Let \(P_A = \{p_i^A\}_{i=1}^n\) be a set of solvent-accessible surface sample points on chain \(A\)
  that are classified as interfacial.
- Let \(P_B = \{p_j^B\}_{j=1}^m\) be the corresponding set on chain \(B\).
- A discrete "facing" pairing is defined by mapping each \(p_i^A\) to its nearest neighbour in \(P_B\):

\[
\pi(i) = \arg \min_{1 \le j \le m} \left\lVert p_i^A - p_j^B \right\rVert
\]

This yields paired points \((p_i^A,\ p_{\pi(i)}^B)\) for \(i = 1,\dots,n\). The number of pairs is:

\[
n_{\mathrm{pairs}} = n
\]

### Electrostatic potential at surface points

The electrostatic potential at a point \(p\) is computed as a Coulomb sum over partial charges:

\[
\phi(p) = \sum_{k} \frac{q_k}{\varepsilon_r\ \max(r_{\min},\ \lVert p - r_k \rVert)}
\]

where:
- \(q_k\) is the partial charge of atom \(k\) (in units of elementary charge \(e\)),
- \(r_k\) is the coordinate of atom \(k\),
- \(\varepsilon_r\) is an effective relative dielectric constant,
- \(r_{\min}\) prevents divergence at short distances.

For an interface, the potentials are evaluated in the *partner field* sense:

- \(\phi_A(i) = \phi_B(p_i^A)\): potential on \(A\)'s surface point due to charges on \(B\),
- \(\phi_B(i) = \phi_A(p_{\pi(i)}^B)\): potential on \(B\)'s facing surface point due to charges on \(A\).

### Per-interface EC (correlation score)

Electrostatic complementarity for the \(A\)–\(B\) interface is defined as the Pearson correlation:

\[
\mathrm{EC}(A,B) = r_{AB} = \mathrm{corr}\big(\phi_A,\ -\phi_B\big)
\]

with:

\[
\mathrm{corr}(x,y) =
\frac{\sum_{i=1}^n (x_i-\bar{x})(y_i-\bar{y})}
{\sqrt{\sum_{i=1}^n (x_i-\bar{x})^2}\ \sqrt{\sum_{i=1}^n (y_i-\bar{y})^2}}
\]

Interpretation:
- \(r_{AB} > 0\): potentials tend to be opposite-sign across the interface (electrostatically complementary).
- \(r_{AB} \approx 0\): weak or no relationship.
- \(r_{AB} < 0\): potentials tend to be same-sign (electrostatically anti-complementary).

### Per-interface EC density

To normalise EC by interface size, an EC density is reported as:

\[
\mathrm{ECdensity}(A,B) = \frac{r_{AB}}{\mathrm{BSA}_{AB}}
\]

where \(\mathrm{BSA}_{AB}\) is the buried surface area of the \(A\)–\(B\) interface (Å\(^2\)).

Units: Å\(^{-2}\).

### Lattice EC: reference chain against multiple partners

In lattice mode, chain \(A\) forms interfaces with multiple partner chains \(B_k\).
The lattice summary reports two weighted EC scores computed from the set of per-partner
EC correlations \(\{r_k\}\).

Because correlations are not additive, weighted aggregation is performed using Fisher’s
z-transform:

\[
z_k = \operatorname{atanh}(r_k)
\]

For weights \(w_k\), the weighted mean is:

\[
\bar{z} = \frac{\sum_k w_k z_k}{\sum_k w_k}
\]

and the aggregated lattice EC is:

\[
r_{\mathrm{lattice}} = \tanh(\bar{z})
\]

Two weighting schemes are reported:

- **BSA-weighted**: \(w_k = \mathrm{BSA}_{A B_k}\)
- **n\_pairs-weighted**: \(w_k = n_{\mathrm{pairs},k}\)

### Lattice EC density

Lattice EC density is reported by normalising by the reference buried area:

\[
\mathrm{ECdensity}_{\mathrm{lattice}} = \frac{r_{\mathrm{lattice}}}{A_{\mathrm{buried}}}
\]

where:

\[
A_{\mathrm{buried}} = \mathrm{SASA}_{\mathrm{iso}}(A) - \mathrm{SASA}_{\mathrm{cluster}}(A)
\]

Here \(\mathrm{SASA}_{\mathrm{iso}}(A)\) is the Shrake–Rupley solvent-accessible surface area of chain \(A\)
in isolation, and \(\mathrm{SASA}_{\mathrm{cluster}}(A)\) is the sum of per-residue SASA for chain \(A\) computed
in the full lattice model (all chains present, occluding the probe).

Units: Å\(^{-2}\).

### Reported fields

Per interface:
- **EC (r)**: \(r_{AB}\)
- **n\_pairs**: \(n_{\mathrm{pairs}}\)
- **EC density**: \(r_{AB}/\mathrm{BSA}_{AB}\)

Lattice summary:
- **EC (r, BSA-weighted Fisher-z)**: aggregated \(r_{\mathrm{lattice}}\) with \(w_k=\mathrm{BSA}\)
- **EC (r, n\_pairs-weighted Fisher-z)**: aggregated \(r_{\mathrm{lattice}}\) with \(w_k=n_{\mathrm{pairs}}\)
- **EC density**: aggregated \(r_{\mathrm{lattice}}/A_{\mathrm{buried}}\)
- **EC by partner chain**: \((r_k,\ n_{\mathrm{pairs},k})\) for each partner.

### Operational note: two-phase lattice EC (same physics)

For **`interface_analyser_lattice_ec.py`**, optional **`--phase sasa`** then **`--phase ec`** splits SASA/BSA work from EC while preserving the same metric definitions; use when runtime or scheduler limits favour two jobs over one. Requirements (matching PDB fingerprint and CLI flags) and examples: **README.md** (Metrics, **`interface_analyser_lattice_ec.py`** subsection).

### Related: Caver tunnel “EC” (different geometry)

`metrics/caver_tunnel_analysis.py` reports a **tunnel**-specific analogue of electrostatic complementarity, built on a **Caver 3.0** centreline and a shell of **lining** residues, using the same class of **Coulomb-summed** potential described above for interfaces. The **geometric** construction is **not** the protein–protein, SASA-facing point pairing used in `interface_analyser_asu_ec.py`, `interface_analyser_lattice_ec.py`, and `electrostatic_complementarity.py`. **Do not** treat tunnel and interface \(r\) values as interchangeable without clear labelling. Full behaviour, options, and outputs are documented in **README.md** (Metrics) and **metrics/metrics_details.md** (Section 1.7).

