#!/usr/bin/env python3
"""
Cross-structure interface comparison: per-interface table, metric matrices, and heatmaps
from FoldKit interface analyser text reports (EC or charge lattice/ASU).

Chain pairs are canonised as ``min(chain_id)-max(chain_id)`` (string sort). Matrix columns
are those canonical labels sorted **lexicographically** (e.g. ``A-B``, ``A-C``, ``B-C``).

The per-interface CSV starts with the low-detail **all-chains** summary columns (when present
in the report), then structure / pair identifiers, **BSA**, **contact counts** (total and,
for EC reports, H-bond / electrostatic / hydrophobic / van der Waals breakdown), **charged**
contact counts for charge-only reports, **EC (r)**, parser **EC density** (r/Å²), and
**EC per 1000 Å²** (rounded to two decimals). Matrices labelled ``ec_density`` use the latter.

**`--heatmap-annotate-metrics`** selects which heatmaps get **one line** of numeric text per cell.
With no **`--heatmap-cell-text`** override, that text matches the colour values. Use
**`--heatmap-cell-text METRIC=FIELD`** to print a different CSV column (e.g. **`ec=ec_density_per_1000_A2`**
for **EC (r)** colouring with **EC per 1000 Å²** labels).

Outputs (stem from ``--prefix`` or ``--output-dir/interface_analysis``):

  - ``*_per_interface.csv``
  - ``*_matrix_<metric>.csv`` + ``*_heatmap_<metric>.<fmt>`` per selected ``--metrics`` name

Default matrices: ``bsa``, ``ec``, ``ec_density``, ``contacts`` (atom contact count). Optional:
``charged_opposite``, ``charged_same`` (charge complementarity reports only; omit from
``--metrics`` when using EC-only inputs to avoid empty figures).

Heatmap styling mirrors ``utils/foldkit_heatmap.py`` / ``ranking/rmsd_to_csv.py``
(``--heatmap-cmap``, ``--heatmap-vmin`` / ``--heatmap-vmax``, ``--heatmap-diverging-center``,
``--heatmap-colorbar-orientation``, ``--heatmap-short-labels``, colour-bar ticks, etc.).
Use **`--heatmap-boundaries`** / **`--heatmap-boundaries-by-metric`** for discrete BSA-style buckets:
comma-separated ascending edges define intervals; the colour-bar labels show each numeric range.
Heatmap **file format** is set with ``--heatmap-format`` (**png**, **svg**, or **pdf**).

When one interface dominates the BSA scale, use **`--heatmap-bsa-robust`** to split the colour
bar between the bulk and the dominant interface(s). A canonical chain pair is auto-flagged as
an outlier when its **minimum BSA across structures** exceeds **`--heatmap-bsa-outlier-factor`**
(default ``3``) times the median of the remaining cells; the cap is the rounded maximum of
non-outlier cells. Override the auto rule with **`--heatmap-bsa-split-at VALUE`** (Å²). Cells
below the cap use ``--heatmap-cmap`` and cells above use a contrasting gradient
(**`--heatmap-bsa-above-cmap`**, default ``Reds``), so above-cap interfaces remain comparable
to each other. Above-cap pairs are listed with their raw Å² values in a figure note (disable
with **`--heatmap-bsa-no-outlier-note`**). For binned readouts use ``--heatmap-boundaries-by-metric
bsa=…``; for log-compressed scales use ``--heatmap-scale-by-metric bsa=log1p``. Explicit
``--heatmap-vmax-by-metric bsa=…`` or ``--heatmap-boundaries-by-metric bsa=…`` settings take
precedence over ``--heatmap-bsa-robust`` (a stderr warning is emitted).

Examples (repository root)::

  PYTHONPATH=. python metrics/interface_analysis_matrix.py \\
    lattice_ec_*.txt --output-dir ./compare/

  PYTHONPATH=. python metrics/interface_analysis_matrix.py \\
    reports/*.txt --prefix ./out/run1 --metrics bsa ec_density contacts
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from typing import Any

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from metrics.interface_mol_report_charge_csv import parse_charge_report_text  # noqa: E402
from metrics.interface_mol_report_ec_csv import parse_ec_report_text  # noqa: E402


SUMMARY_KEYS = (
    "all_total_interfaces",
    "all_total_buried_surface_area_A2",
    "all_average_buried_area_per_interface_A2",
    "all_sasa_isolated_sum_A2",
)

PER_IFACE_COLUMNS = (
    *SUMMARY_KEYS,
    "report_txt",
    "structure_basename",
    "chain_pair_canonical",
    "chain_pair_raw",
    "interface_number",
    "buried_surface_area_A2",
    "contact_count",
    "contact_area_A2",
    "hbonds",
    "electrostatic",
    "hydrophobic",
    "van_der_waals",
    "charged_contacts_opposite",
    "charged_contacts_same",
    "ec_r",
    "ec_density",
    "ec_density_denominator",
    "ec_density_per_1000_A2",
)

METRIC_CHOICES: tuple[str, ...] = (
    "bsa",
    "ec",
    "ec_density",
    "contacts",
    "charged_opposite",
    "charged_same",
)

DEFAULT_METRICS: tuple[str, ...] = ("bsa", "ec", "ec_density", "contacts")

_METRIC_FIELD = {
    "bsa": "buried_surface_area_A2",
    "ec": "ec_r",
    "ec_density": "ec_density_per_1000_A2",
    "contacts": "contact_count",
    "charged_opposite": "charged_contacts_opposite",
    "charged_same": "charged_contacts_same",
}

_METRIC_STEM_SUFFIX = {
    "bsa": "bsa",
    "ec": "ec",
    "ec_density": "ec_density",
    "contacts": "contacts",
    "charged_opposite": "charged_contacts_opposite",
    "charged_same": "charged_contacts_same",
}

_METRIC_HEATMAP_LABEL = {
    "bsa": "Buried surface area (Å²)",
    "ec": "EC (r)",
    "ec_density": "EC per 1000 Å² (r)",
    "contacts": "Contact count",
    "charged_opposite": "Opposite-sign charged contacts",
    "charged_same": "Same-sign charged contacts",
}

_METRIC_TITLE = {
    "bsa": "Pairwise interface BSA across structures",
    "ec": "Pairwise interface EC (r) across structures",
    "ec_density": "Pairwise interface EC per 1000 Å² across structures",
    "contacts": "Pairwise interface contact counts across structures",
    "charged_opposite": "Pairwise opposite-sign charged contacts across structures",
    "charged_same": "Pairwise same-sign charged contacts across structures",
}


def _ec_density_per_1000_A2(raw: Any) -> Any:
    """EC per Å² × 1000 → EC per 1000 Å²; rounded to 2 decimals. Empty if not parseable."""
    if raw is None or raw == "":
        return ""
    try:
        return round(float(raw) * 1000.0, 2)
    except (TypeError, ValueError):
        return ""


def _parse_csv_floats(s: str | None) -> list[float] | None:
    if s is None or not str(s).strip():
        return None
    return [float(x.strip()) for x in str(s).split(",") if x.strip()]


def _structure_axis_label(sb: str, *, short_fn, use_short: bool) -> str:
    """Structure tick labels for heatmaps.

    When ``--heatmap-short-labels`` is set, use the structure basename up to the first underscore,
    matching the common ``sXXX_*`` naming used by FoldKit pipelines.
    """
    if not use_short:
        return sb
    raw = os.path.splitext(os.path.basename(sb))[0]
    head = raw.split("_", 1)[0].strip()
    return head or raw or sb


def _heatmap_format_choice(s: str) -> str:
    """Backward-compatible format normaliser (foldkit_heatmap owns the canonical CLI)."""
    t = str(s).strip().lower().lstrip(".")
    if t in ("png", "pdf", "svg"):
        return t
    raise argparse.ArgumentTypeError(f"heatmap format must be one of: png, svg, pdf (got {s!r})")


def _import_foldkit_heatmap_helpers():
    """Heatmap helpers (same module as ``rmsd_to_csv`` / ``foldkit_heatmap`` CLI)."""
    from utils.foldkit_heatmap import (  # noqa: E402
        _add_raster_colorbar,
        _apply_heatmap_y_axis_right,
        _draw_heatmap_vector_cells,
        _heatmap_vector_format,
        _savefig_vector_heatmap,
        add_generic_heatmap_args,
        apply_heatmap_scale,
        boundary_bin_centres_and_labels,
        heatmap_should_rotate_xticklabels,
        heatmap_ticks_for_norm,
        make_boundary_norm_cmap,
        make_heatmap_norm,
        parse_boundaries_csv,
        short_heatmap_label,
    )

    return (
        make_heatmap_norm,
        short_heatmap_label,
        _apply_heatmap_y_axis_right,
        _heatmap_vector_format,
        _draw_heatmap_vector_cells,
        _add_raster_colorbar,
        _savefig_vector_heatmap,
        add_generic_heatmap_args,
        apply_heatmap_scale,
        boundary_bin_centres_and_labels,
        heatmap_should_rotate_xticklabels,
        heatmap_ticks_for_norm,
        make_boundary_norm_cmap,
        parse_boundaries_csv,
    )


def _median_finite_2d(arr) -> float | None:
    import numpy as np

    v = arr[np.isfinite(arr)]
    if v.size == 0:
        return None
    return float(np.median(v))


def _data_range_finite(
    arr,
    *,
    positive_only: bool,
) -> tuple[float | None, float | None]:
    import numpy as np

    v = arr[np.isfinite(arr)]
    if v.size == 0:
        return None, None
    if positive_only:
        vp = v[v > 0]
        if vp.size > 0:
            v = vp
    return float(np.min(v)), float(np.max(v))


def _detect_bsa_outlier_pairs(
    pair_to_values: dict[str, list[float]],
    *,
    factor: float,
) -> list[str]:
    """
    Iteratively flag canonical chain pairs whose minimum positive BSA across structures exceeds
    ``factor * median(remaining positive cells)``. Returns outlier pair labels in detection order.

    Each iteration the candidate with the largest qualifying minimum is removed, then the median
    is recomputed over the cells that remain. This handles the case where multiple interfaces are
    each individually dominant: once one is removed, the median drops and the next interface can
    qualify as an outlier in turn.
    """
    import math

    k = float(factor)
    if not math.isfinite(k) or k <= 1.0:
        return []
    pair_pos: dict[str, list[float]] = {}
    for p, vs in pair_to_values.items():
        clean = [float(v) for v in vs if v is not None and math.isfinite(float(v)) and float(v) > 0]
        if clean:
            pair_pos[str(p)] = clean
    if not pair_pos:
        return []

    remaining = set(pair_pos.keys())
    outliers: list[str] = []
    while remaining:
        cells: list[float] = []
        for p in remaining:
            cells.extend(pair_pos[p])
        if not cells:
            break
        cells_sorted = sorted(cells)
        n = len(cells_sorted)
        m = (
            cells_sorted[n // 2]
            if n % 2 == 1
            else 0.5 * (cells_sorted[n // 2 - 1] + cells_sorted[n // 2])
        )
        threshold = k * float(m)
        candidate: str | None = None
        best_min = float("-inf")
        for p in remaining:
            mp = min(pair_pos[p])
            if mp > threshold and mp > best_min:
                candidate, best_min = p, mp
        if candidate is None:
            break
        outliers.append(candidate)
        remaining.discard(candidate)
    return outliers


def _cbar_ticks_from_step(vmin: float, vmax: float, step: float) -> list[float] | None:
    import numpy as np

    if step <= 0:
        return None
    lo, hi = float(vmin), float(vmax)
    if hi < lo:
        lo, hi = hi, lo
    start = np.floor(lo / step) * step
    vals = np.arange(start, hi + 0.5 * step, step, dtype=float)
    vals = vals[(vals >= lo - 1e-9) & (vals <= hi + 1e-9)]
    return [float(x) for x in vals]


def _apply_custom_order(items: list[str], wanted: list[str]) -> list[str]:
    """
    Return a new list where items mentioned in wanted come first (in that order),
    followed by any remaining items in their original order.
    """
    if not wanted:
        return list(items)
    s = set(items)
    bad = [x for x in wanted if x not in s]
    if bad:
        raise ValueError(f"Unknown x-axis item(s): {', '.join(bad)}")
    out: list[str] = []
    seen: set[str] = set()
    for x in wanted:
        if x not in seen:
            out.append(x)
            seen.add(x)
    for x in items:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out


def _parse_csv_list(s: str | None) -> list[str]:
    if s is None or not str(s).strip():
        return []
    return [x.strip() for x in str(s).split(",") if x.strip()]


def _apply_custom_order_by_key(items: list[str], wanted: list[str], *, key_fn) -> list[str]:
    """
    Like _apply_custom_order, but match wanted against key_fn(item).

    Requires that key_fn(item) is unique across items.
    """
    if not wanted:
        return list(items)
    key_to_item: dict[str, str] = {}
    dup: dict[str, list[str]] = {}
    for it in items:
        k = str(key_fn(it))
        if k in key_to_item and key_to_item[k] != it:
            dup.setdefault(k, [key_to_item[k]]).append(it)
        key_to_item[k] = it
    if dup:
        ks = ", ".join(sorted(dup.keys()))
        raise ValueError(f"Cannot use --heatmap-x-order with non-unique short labels: {ks}")
    translated: list[str] = []
    missing: list[str] = []
    for w in wanted:
        if w in key_to_item:
            translated.append(key_to_item[w])
        else:
            missing.append(w)
    if missing:
        raise ValueError(f"Unknown x-axis item(s): {', '.join(missing)}")
    return _apply_custom_order(items, translated)


def _canon_chain_pair(c1: str, c2: str) -> str:
    a = str(c1 or "").strip()
    b = str(c2 or "").strip()
    if not a or not b:
        return ""
    x, y = sorted((a, b))
    return f"{x}-{y}"


def _expand_inputs(raw: list[str]) -> list[str]:
    out: list[str] = []
    for item in raw:
        item = os.path.expanduser(item.strip())
        if not item:
            continue
        if any(ch in item for ch in "*?["):
            found = sorted(glob.glob(item))
            if not found:
                print(f"Warning: no files matched glob {item!r}", file=sys.stderr)
            out.extend(found)
        elif os.path.isfile(item):
            out.append(item)
        else:
            print(f"Warning: not a file, skipping {item!r}", file=sys.stderr)
    seen: set[str] = set()
    uniq: list[str] = []
    for p in out:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            uniq.append(ap)
    return uniq


def _parse_report_text(text: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ifaces, sums = parse_ec_report_text(text)
    if not ifaces:
        ifaces, sums = parse_charge_report_text(text)
    return ifaces, sums


def _summary_by_structure(summaries: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    m: dict[str, dict[str, Any]] = {}
    for s in summaries:
        key = str(s.get("structure_basename") or "").strip()
        if not key:
            continue
        m[key] = dict(s)
    return m


def _collect_rows(
    report_paths: list[str],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    global_summary: dict[str, dict[str, Any]] = {}

    for path in report_paths:
        rtxt = os.path.basename(path)
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
        ifaces, sums = _parse_report_text(text)
        local_sum = _summary_by_structure(sums)
        for k, v in local_sum.items():
            global_summary[k] = v

        for iface in ifaces:
            sb = str(iface.get("structure_basename") or "").strip()
            sum_rec = local_sum.get(sb, {})
            pair_raw = str(iface.get("chain_pair") or "").strip()
            c1 = str(iface.get("chain1_id") or "").strip()
            c2 = str(iface.get("chain2_id") or "").strip()
            canon = _canon_chain_pair(c1, c2)
            row: dict[str, Any] = {
                "report_txt": rtxt,
                "structure_basename": sb,
                "chain_pair_raw": pair_raw,
                "chain_pair_canonical": canon,
                "interface_number": iface.get("interface_number", ""),
                "buried_surface_area_A2": iface.get("buried_surface_area_A2", ""),
                "contact_count": iface.get("contact_count", ""),
                "contact_area_A2": iface.get("contact_area_A2", ""),
                "hbonds": iface.get("hbonds", ""),
                "electrostatic": iface.get("electrostatic", ""),
                "hydrophobic": iface.get("hydrophobic", ""),
                "van_der_waals": iface.get("van_der_waals", ""),
                "charged_contacts_opposite": iface.get("charged_contacts_opposite", ""),
                "charged_contacts_same": iface.get("charged_contacts_same", ""),
                "ec_r": iface.get("ec_r", ""),
                "ec_density": iface.get("ec_density", ""),
                "ec_density_denominator": iface.get("ec_density_denominator", ""),
                "ec_density_per_1000_A2": _ec_density_per_1000_A2(iface.get("ec_density")),
            }
            for sk in SUMMARY_KEYS:
                row[sk] = sum_rec.get(sk, "")
            rows.append(row)

    return rows, global_summary


def _dedupe_iface_per_structure(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    key_to_row: dict[tuple[str, str], dict[str, Any]] = {}
    order: list[tuple[str, str]] = []
    for r in rows:
        sb = str(r.get("structure_basename") or "").strip()
        cp = str(r.get("chain_pair_canonical") or "").strip()
        if not sb or not cp:
            continue
        k = (sb, cp)
        if k not in key_to_row:
            order.append(k)
        key_to_row[k] = r
    return [key_to_row[k] for k in order]


def _ordered_canonical_pairs(rows: list[dict[str, Any]]) -> list[str]:
    pairs_set = {
        str(r["chain_pair_canonical"]) for r in rows if str(r.get("chain_pair_canonical") or "").strip()
    }
    return sorted(pairs_set)


def _value_by_struct_pair(
    rows: list[dict[str, Any]],
    field: str,
) -> dict[tuple[str, str], Any]:
    out: dict[tuple[str, str], Any] = {}
    for r in rows:
        sb = str(r.get("structure_basename") or "").strip()
        cp = str(r.get("chain_pair_canonical") or "").strip()
        if not sb or not cp:
            continue
        out[(sb, cp)] = r.get(field, "")
    return out


def _pair_matrix_numpy(
    structures: list[str],
    ordered_pairs: list[str],
    transpose: bool,
    value_map: dict[tuple[str, str], Any],
):
    """Float matrix [n_row x n_col] aligned with transpose layout (same indexing as heatmap colour)."""
    import numpy as np

    n_row = len(ordered_pairs) if transpose else len(structures)
    n_col = len(structures) if transpose else len(ordered_pairs)
    arr = np.full((n_row, n_col), np.nan, dtype=float)
    if not transpose:
        for i, sb in enumerate(structures):
            for j, p in enumerate(ordered_pairs):
                v = value_map.get((sb, p), None)
                if v is None or v == "":
                    continue
                try:
                    arr[i, j] = float(v)
                except (TypeError, ValueError):
                    arr[i, j] = np.nan
    else:
        for i, p in enumerate(ordered_pairs):
            for j, sb in enumerate(structures):
                v = value_map.get((sb, p), None)
                if v is None or v == "":
                    continue
                try:
                    arr[i, j] = float(v)
                except (TypeError, ValueError):
                    arr[i, j] = np.nan
    return arr


_CELL_TEXT_FIELD_ALLOWED = frozenset((*PER_IFACE_COLUMNS, "same"))


def _heatmap_cell_text_map_from_cli(items: list[str]) -> dict[str, str]:
    """
    Parse METRIC=FIELD tokens. FIELD must be 'same' or a per-interface CSV column name.
    """
    out: dict[str, str] = {}
    for raw in items:
        item = raw.strip()
        if "=" not in item:
            raise ValueError(f"Invalid --heatmap-cell-text {raw!r}; expected METRIC=FIELD.")
        mk, fk = item.split("=", 1)
        mk, fk = mk.strip(), fk.strip()
        if mk not in METRIC_CHOICES:
            raise ValueError(
                f"Unknown heatmap metric {mk!r} in --heatmap-cell-text (choices: {', '.join(METRIC_CHOICES)})."
            )
        if fk not in _CELL_TEXT_FIELD_ALLOWED:
            raise ValueError(
                f"Unknown annotation field {fk!r} in --heatmap-cell-text "
                f"(use 'same' or a per-interface CSV column name)."
            )
        out[mk] = fk
    return out


def _metric_float_map_from_cli(items: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for raw in items:
        item = (raw or "").strip()
        if "=" not in item:
            raise ValueError(f"Invalid token {raw!r}; expected METRIC=NUMBER.")
        mk, vs = item.split("=", 1)
        mk, vs = mk.strip(), vs.strip()
        if mk not in METRIC_CHOICES:
            raise ValueError(f"Unknown heatmap metric {mk!r} (choices: {', '.join(METRIC_CHOICES)}).")
        try:
            out[mk] = float(vs)
        except ValueError:
            raise ValueError(f"Invalid numeric value in {raw!r}; expected METRIC=NUMBER.") from None
    return out


def _metric_str_map_from_cli(items: list[str], *, allowed: set[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in items:
        item = (raw or "").strip()
        if "=" not in item:
            raise ValueError(f"Invalid token {raw!r}; expected METRIC=VALUE.")
        mk, vs = item.split("=", 1)
        mk, vs = mk.strip(), vs.strip()
        if mk not in METRIC_CHOICES:
            raise ValueError(f"Unknown heatmap metric {mk!r} (choices: {', '.join(METRIC_CHOICES)}).")
        if vs not in allowed:
            raise ValueError(f"Invalid value {vs!r} for {mk!r} (allowed: {', '.join(sorted(allowed))}).")
        out[mk] = vs
    return out


def _metric_edges_map_from_cli(items: list[str]) -> dict[str, list[float]]:
    from utils.foldkit_heatmap import parse_boundaries_csv  # noqa: E402

    out: dict[str, list[float]] = {}
    for raw in items:
        item = (raw or "").strip()
        if "=" not in item:
            raise ValueError(f"Invalid token {raw!r}; expected METRIC=V0,V1,…")
        mk, rest = item.split("=", 1)
        mk, rest = mk.strip(), rest.strip()
        if mk not in METRIC_CHOICES:
            raise ValueError(f"Unknown heatmap metric {mk!r} (choices: {', '.join(METRIC_CHOICES)}).")
        edges = parse_boundaries_csv(rest)
        if len(edges) < 2:
            raise ValueError(f"Need at least two strictly distinct boundary edges for {mk!r} (got {rest!r}).")
        out[mk] = edges
    return out


def _write_per_interface_csv(path: str, rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    fieldnames = list(PER_IFACE_COLUMNS)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in sorted(
            rows,
            key=lambda x: (
                str(x.get("structure_basename") or ""),
                str(x.get("chain_pair_canonical") or ""),
            ),
        ):
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _orient_plot_matrix(a, *, transpose: bool, vector_cells: bool):
    """Return matrix for plotting with a consistent origin='lower' layout.

    Both the raster (PNG) and vector (SVG/PDF) paths are rendered with row 0 at the bottom
    and column 0 at the left so axis labels always run bottom→top / left→right in list order.
    """
    import numpy as np

    x = np.array(a, dtype=float, copy=True)
    return x


def _heatmap_row_col_labels(
    structures: list[str],
    ordered_pairs: list[str],
    *,
    transpose: bool,
    vector_cells: bool,
    short_heatmap_label,
    args: argparse.Namespace,
) -> tuple[list[str], list[str]]:
    short = bool(getattr(args, "heatmap_short_labels", False))
    if not transpose:
        row_labels = list(structures)
        col_labels = list(ordered_pairs)
        if short:
            row_labels = [
                _structure_axis_label(str(x), short_fn=short_heatmap_label, use_short=True)
                for x in row_labels
            ]
    else:
        row_labels = list(ordered_pairs)
        col_labels = list(structures)
        if short:
            col_labels = [
                _structure_axis_label(str(x), short_fn=short_heatmap_label, use_short=True)
                for x in col_labels
            ]
    return row_labels, col_labels


def _write_matrix_csv(
    path: str,
    structures: list[str],
    ordered_pairs: list[str],
    values: dict[tuple[str, str], Any],
    *,
    transpose: bool,
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    # Default: rows=structures, cols=interfaces (canonical pairs).
    # Transpose: rows=interfaces, cols=structures.
    if not transpose:
        headers = ["structure_basename", *ordered_pairs]
    else:
        headers = ["chain_pair_canonical", *structures]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        if not transpose:
            for sb in structures:
                row_out: list[Any] = [sb]
                for p in ordered_pairs:
                    row_out.append(values.get((sb, p), ""))
                w.writerow(row_out)
        else:
            for p in ordered_pairs:
                row_out = [p]
                for sb in structures:
                    row_out.append(values.get((sb, p), ""))
                w.writerow(row_out)


def _write_heatmap(
    path: str,
    structures: list[str],
    ordered_pairs: list[str],
    values: dict[tuple[str, str], Any],
    *,
    rows: list[dict[str, Any]],
    metric_key: str,
    args: argparse.Namespace,
    cbar_label: str,
    title: str,
    transpose: bool,
) -> str | None:
    try:
        import numpy as np
    except ImportError:
        print(f"Warning: NumPy not installed; skipping heatmap {path!r}.", file=sys.stderr)
        return None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"Warning: Matplotlib not installed; skipping heatmap {path!r}.", file=sys.stderr)
        return None

    try:
        (
            make_heatmap_norm,
            short_heatmap_label,
            apply_y_right,
            is_vector_fmt,
            draw_vector_cells,
            add_raster_cbar,
            savefig_vector,
            _add_generic_heatmap_args,
            apply_heatmap_scale,
            boundary_bin_centres_and_labels,
            should_rotate_xlabels,
            ticks_for_norm,
            make_boundary_norm_cmap,
            parse_boundaries_csv,
        ) = _import_foldkit_heatmap_helpers()
    except ImportError as e:
        print(
            f"Warning: could not import utils.foldkit_heatmap helpers ({e}); skipping heatmap {path!r}.",
            file=sys.stderr,
        )
        return None

    n_row = len(ordered_pairs) if transpose else len(structures)
    n_col = len(structures) if transpose else len(ordered_pairs)
    if n_row == 0 or n_col == 0:
        print(f"Warning: empty matrix; skipping heatmap {path!r}.", file=sys.stderr)
        return None

    arr = _pair_matrix_numpy(structures, ordered_pairs, transpose, values)
    scale_by: dict[str, str] = getattr(args, "heatmap_scale_by_metric", None) or {}
    scale_mode = str(scale_by.get(metric_key, getattr(args, "heatmap_scale", "linear")) or "linear")
    try:
        arr, lab_pre, lab_suf = apply_heatmap_scale(arr, scale_mode)
    except ValueError as exc:
        print(f"Warning: {exc}; using linear scale for {path!r}.", file=sys.stderr)
        lab_pre, lab_suf = "", ""

    cmap_name = str(getattr(args, "heatmap_cmap", "viridis_r") or "viridis_r").strip()

    positive_only = metric_key == "bsa" and bool(getattr(args, "heatmap_bsa_positive_range", False))
    data_min, data_max = _data_range_finite(arr, positive_only=positive_only)
    if data_min is None or data_max is None:
        print(f"Warning: no finite values for heatmap {path!r}; skipping.", file=sys.stderr)
        return None

    vmin_by: dict[str, float] = getattr(args, "heatmap_vmin_by_metric", None) or {}
    vmax_by: dict[str, float] = getattr(args, "heatmap_vmax_by_metric", None) or {}
    vmin_u = vmin_by.get(metric_key, getattr(args, "heatmap_vmin", None))
    vmax_u = vmax_by.get(metric_key, getattr(args, "heatmap_vmax", None))

    bd_global = parse_boundaries_csv(getattr(args, "heatmap_boundaries", None))
    bd_by_metric: dict[str, list[float]] = getattr(args, "heatmap_boundaries_by_metric", None) or {}
    edges_metric = bd_by_metric.get(metric_key)
    edges_use = edges_metric if edges_metric else bd_global
    use_boundary = bool(edges_use and len(edges_use) >= 2)

    bsa_robust_requested = bool(metric_key == "bsa" and getattr(args, "heatmap_bsa_robust", False))
    bsa_robust_active = False
    bsa_robust_cap: float | None = None
    bsa_robust_outliers: list[str] = []
    bsa_robust_outlier_factor: float = float(getattr(args, "heatmap_bsa_outlier_factor", 3.0))
    bsa_robust_split_at: float | None = getattr(args, "heatmap_bsa_split_at", None)
    if bsa_robust_split_at is not None:
        try:
            bsa_robust_split_at = float(bsa_robust_split_at)
        except (TypeError, ValueError):
            bsa_robust_split_at = None
    if bsa_robust_requested:
        if use_boundary:
            print(
                f"Warning: --heatmap-bsa-robust ignored because boundary buckets "
                f"(--heatmap-boundaries* for {metric_key!r}) are set; the boundary settings win.",
                file=sys.stderr,
            )
        elif vmax_u is not None:
            print(
                f"Warning: --heatmap-bsa-robust ignored because an explicit "
                f"--heatmap-vmax-by-metric bsa={vmax_u} was provided; the explicit value wins.",
                file=sys.stderr,
            )
        elif bsa_robust_split_at is None and not (bsa_robust_outlier_factor > 1.0):
            print(
                f"Warning: --heatmap-bsa-outlier-factor must be > 1 "
                f"(got {bsa_robust_outlier_factor}); using 3.0.",
                file=sys.stderr,
            )
            bsa_robust_outlier_factor = 3.0

    ticks_boundary: list[float] | None = None
    tick_labels_boundary: list[str] | None = None
    cmap_obj = None

    if use_boundary:
        if vmin_u is not None or vmax_u is not None:
            print(
                f"Warning: vmin/vmax ignored for heatmap {metric_key} ({path!r}) "
                "because boundary buckets (--heatmap-boundaries*) are set.",
                file=sys.stderr,
            )
        try:
            norm, cmap_obj = make_boundary_norm_cmap(edges_use, cmap_name)
        except ValueError as exc:
            print(f"Warning: invalid boundaries for heatmap {path!r}: {exc}; skipping.", file=sys.stderr)
            return None
        try:
            ticks_boundary, tick_labels_boundary = boundary_bin_centres_and_labels(edges_use)
        except ValueError as exc:
            print(f"Warning: {exc}; skipping heatmap {path!r}.", file=sys.stderr)
            return None
    else:
        disp_vmin = float(vmin_u) if vmin_u is not None else data_min
        disp_vmax = float(vmax_u) if vmax_u is not None else data_max
        if disp_vmin > disp_vmax:
            print(
                f"Error: heatmap scale requires vmin < vmax (got {disp_vmin} {disp_vmax}) for {path!r}.",
                file=sys.stderr,
            )
            return None
        if disp_vmax - disp_vmin <= 1e-15:
            mid = 0.5 * (disp_vmin + disp_vmax)
            eps = 1e-6 * max(abs(mid), 1.0)
            disp_vmin, disp_vmax = mid - eps, mid + eps

        # BSA robust mode: split the colour bar between the bulk and dominant interface(s).
        # Done before the diverging-centre handling so the cap takes precedence over a median
        # centre. Manual --heatmap-bsa-split-at wins over auto detection.
        if bsa_robust_requested and not use_boundary and vmax_u is None:
            from utils.foldkit_heatmap import nice_round_up  # noqa: E402

            try:
                v_pos = arr[np.isfinite(arr) & (arr > 0)]
            except Exception:
                v_pos = None
            if v_pos is not None and v_pos.size:
                v_min_pos = float(np.nanmin(v_pos))
                v_max_pos = float(np.nanmax(v_pos))

                # Group positive BSA values by canonical chain pair (uses raw rows for the
                # parsed BSA, which is more robust than reverse-engineering the matrix).
                pair_to_values: dict[str, list[float]] = {}
                for r in rows:
                    cp = str(r.get("chain_pair_canonical") or "").strip()
                    if not cp:
                        continue
                    try:
                        bv = float(r.get("buried_surface_area_A2") or "")
                    except (TypeError, ValueError):
                        continue
                    if not (np.isfinite(bv) and bv > 0):
                        continue
                    pair_to_values.setdefault(cp, []).append(bv)

                if bsa_robust_split_at is not None:
                    cap_val = float(bsa_robust_split_at)
                    if np.isfinite(cap_val) and v_min_pos < cap_val < v_max_pos:
                        bsa_robust_active = True
                        bsa_robust_cap = cap_val
                        disp_vmin = v_min_pos
                        disp_vmax = v_max_pos
                        bsa_robust_outliers = sorted(
                            p for p, vs in pair_to_values.items() if min(vs) > cap_val
                        )
                        print(
                            f"Heatmap (bsa): robust split at cap = {cap_val:g} "
                            f"(--heatmap-bsa-split-at).",
                            file=sys.stderr,
                        )
                    else:
                        print(
                            f"Heatmap (bsa): --heatmap-bsa-split-at={cap_val:g} outside data "
                            f"range ({v_min_pos:g}, {v_max_pos:g}); robust mode not applied.",
                            file=sys.stderr,
                        )
                else:
                    outliers = _detect_bsa_outlier_pairs(
                        pair_to_values, factor=bsa_robust_outlier_factor
                    )
                    if outliers:
                        non_outlier_cells: list[float] = []
                        for p, vs in pair_to_values.items():
                            if p not in outliers:
                                non_outlier_cells.extend(vs)
                        if non_outlier_cells:
                            cap_raw = max(non_outlier_cells)
                            cap_nice = nice_round_up(cap_raw)
                            if np.isfinite(cap_nice) and v_min_pos < cap_nice < v_max_pos:
                                bsa_robust_active = True
                                bsa_robust_cap = float(cap_nice)
                                bsa_robust_outliers = list(outliers)
                                disp_vmin = v_min_pos
                                disp_vmax = v_max_pos
                                print(
                                    f"Heatmap (bsa): robust split at cap = {cap_nice:g} "
                                    f"(max non-outlier = {cap_raw:g}); outlier interfaces: "
                                    f"{', '.join(outliers)} (k = {bsa_robust_outlier_factor:g}).",
                                    file=sys.stderr,
                                )
                            else:
                                print(
                                    f"Heatmap (bsa): auto cap {cap_nice:g} outside data range "
                                    f"({v_min_pos:g}, {v_max_pos:g}); robust mode not applied.",
                                    file=sys.stderr,
                                )
                        else:
                            print(
                                "Heatmap (bsa): all interfaces flagged as outliers; robust mode "
                                "not applied. Try a larger --heatmap-bsa-outlier-factor or use "
                                "--heatmap-bsa-split-at.",
                                file=sys.stderr,
                            )
                    else:
                        print(
                            f"Heatmap (bsa): no interface qualifies as outlier "
                            f"(k = {bsa_robust_outlier_factor:g}); using linear scale. Pass "
                            f"--heatmap-bsa-split-at VALUE or lower "
                            f"--heatmap-bsa-outlier-factor to force a split.",
                            file=sys.stderr,
                        )

        if bsa_robust_active and bsa_robust_cap is not None:
            norm = make_heatmap_norm(disp_vmin, disp_vmax, "median", bsa_robust_cap)
        else:
            dc_raw = str(getattr(args, "heatmap_diverging_center", "none") or "none").strip().lower()
            vc_for_norm: float | None = None
            if dc_raw == "median":
                vc_for_norm = _median_finite_2d(arr)
                if vc_for_norm is None:
                    print(f"Warning: could not compute median for heatmap {path!r}; using linear scale.", file=sys.stderr)
                    dc_use = "none"
                else:
                    dc_use = "median"
                    print(f"Heatmap ({metric_key}): median centre = {vc_for_norm:.4g}", file=sys.stderr)
            else:
                dc_use = "none"

            norm = make_heatmap_norm(disp_vmin, disp_vmax, dc_use, vc_for_norm if dc_use == "median" else None)

    cb_orient = str(getattr(args, "heatmap_colorbar_orientation", "vertical") or "vertical").lower()
    figscale = float(getattr(args, "heatmap_figscale", 1.0))
    cell_in = 0.52 * figscale
    ml, mr, mt, mb = 1.55, 1.45, 1.15, 1.25
    if cb_orient == "vertical":
        mr += 1.15
    else:
        mb += 1.35
    fig_w = ml + mr + cell_in * max(n_col, 1)
    fig_h = mt + mb + cell_in * max(n_row, 1)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    if not use_boundary:
        if bsa_robust_active:
            from utils.foldkit_heatmap import compose_split_cmap  # noqa: E402

            above_name = str(getattr(args, "heatmap_bsa_above_cmap", "Reds") or "Reds")
            cmap_obj = compose_split_cmap(below_name=cmap_name, above_name=above_name)
            cmap_obj.set_bad(color="#e0e0e0")
        else:
            try:
                cmap_obj = plt.get_cmap(cmap_name).copy()
            except ValueError:
                print(
                    f"Warning: unknown colour map {cmap_name!r}; using viridis_r.",
                    file=sys.stderr,
                )
                cmap_obj = plt.get_cmap("viridis_r").copy()
            cmap_obj.set_bad(color="#e0e0e0")

    fmt_out = str(getattr(args, "heatmap_format", "png") or "png").lower().lstrip(".")
    vector_cells = bool(is_vector_fmt(fmt_out))

    plot_arr = _orient_plot_matrix(arr, transpose=transpose, vector_cells=vector_cells)

    ann_set = frozenset(getattr(args, "heatmap_annotate_metrics", ()) or ())
    want_cell_ann = metric_key in ann_set
    cell_text_map: dict[str, str] = getattr(args, "heatmap_cell_text_map", None) or {}
    text_spec = cell_text_map.get(metric_key)

    single_ann_plot = plot_arr
    single_cell_ann = False
    if want_cell_ann:
        single_cell_ann = True
        if text_spec is None or text_spec == "same":
            single_ann_plot = plot_arr
        else:
            vm = _value_by_struct_pair(rows, text_spec)
            ann_raw = _pair_matrix_numpy(structures, ordered_pairs, transpose, vm)
            single_ann_plot = _orient_plot_matrix(ann_raw, transpose=transpose, vector_cells=vector_cells)

    row_labels, col_labels = _heatmap_row_col_labels(
        structures,
        ordered_pairs,
        transpose=transpose,
        vector_cells=vector_cells,
        short_heatmap_label=short_heatmap_label,
        args=args,
    )

    annotate_fs = float(getattr(args, "heatmap_annotate_fontsize", 7.0))
    single_ann_fmt = str(getattr(args, "heatmap_annotate_fmt", "{:.2f}") or "{:.2f}")
    _INT_ANN_FIELDS = {
        "contact_count",
        "hbonds",
        "electrostatic",
        "hydrophobic",
        "van_der_waals",
        "charged_contacts_opposite",
        "charged_contacts_same",
    }
    if single_cell_ann:
        # Ensure contact-like counts display as integers even if the global annotation format is decimal.
        if (text_spec in _INT_ANN_FIELDS) or (text_spec is None and metric_key in ("contacts", "charged_opposite", "charged_same")):
            single_ann_fmt = "{:.0f}"

    if vector_cells:
        if single_cell_ann:
            draw_vector_cells(
                ax,
                plot_arr,
                cmap_obj,
                norm,
                annotate_value_arr=single_ann_plot,
                annotate_fmt=single_ann_fmt,
                annotate_fontsize=annotate_fs,
                annotate_offdiag_only=False,
            )
        else:
            draw_vector_cells(
                ax,
                plot_arr,
                cmap_obj,
                norm,
                annotate_value_arr=None,
            )
        ax.set_xlim(-0.5, n_col - 0.5)
        ax.set_ylim(-0.5, n_row - 0.5)
    else:
        # PNG (and fallback): imshow is fine for raster formats.
        from utils.foldkit_heatmap import (  # noqa: E402
            _CELL_BORDER_COLOR,
            _CELL_BORDER_LINEWIDTH_PT,
        )

        ax.imshow(
            np.ma.masked_invalid(plot_arr),
            aspect="auto",
            cmap=cmap_obj,
            norm=norm,
            interpolation="nearest",
            origin="lower",
        )
        ax.set_xticks(np.arange(-0.5, n_col, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_row, 1), minor=True)
        ax.grid(
            which="minor",
            color=_CELL_BORDER_COLOR,
            linestyle="-",
            linewidth=float(_CELL_BORDER_LINEWIDTH_PT),
        )
        ax.tick_params(which="minor", bottom=False, left=False)
        if single_cell_ann:
            from utils.foldkit_heatmap import _format_cell_value, _relative_luminance_srgb  # noqa: E402

            for i in range(n_row):
                for j in range(n_col):
                    tv = single_ann_plot[i, j]
                    if not np.isfinite(tv):
                        continue
                    cv = plot_arr[i, j]
                    if not np.isfinite(cv):
                        continue
                    fc = cmap_obj(norm(float(cv)))
                    luma = _relative_luminance_srgb(fc)
                    tc = "white" if luma <= 0.45 else "0.15"
                    txt = _format_cell_value(float(tv), single_ann_fmt)
                    ax.text(
                        j,
                        i,
                        txt,
                        ha="center",
                        va="center",
                        fontsize=annotate_fs,
                        color=tc,
                        clip_on=True,
                    )

    # Square cells (equal x/y scale), matching utils/foldkit_heatmap.
    ax.set_aspect("equal", adjustable="box")

    axis_fs = float(getattr(args, "heatmap_axis_labelsize", 9.0))
    ax.set_xticks(range(n_col))
    rotate_x = bool(should_rotate_xlabels(list(col_labels), ncols=n_col))
    ax.set_xticklabels(col_labels, rotation=(45 if rotate_x else 0), ha=("right" if rotate_x else "center"))
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(row_labels)
    ax.tick_params(axis="both", which="major", labelsize=axis_fs)
    if transpose:
        ax.set_xlabel("Structure")
        ax.set_ylabel("Interface (canonical chain pair)")
    else:
        ax.set_xlabel("Interface (canonical chain pair)")
        ax.set_ylabel("Structure")
    ax.set_title(title)

    if bool(getattr(args, "heatmap_y_axis_right", False)):
        apply_y_right(ax)

    y_right = bool(getattr(args, "heatmap_y_axis_right", False))
    if cb_orient == "horizontal":
        r_edge = 0.92 if y_right else 0.96
        plt.tight_layout(rect=(0.06, 0.22, r_edge, 0.94))
    elif y_right and cb_orient == "vertical":
        plt.tight_layout(rect=(0.14, 0.06, 0.92, 0.92))
    else:
        plt.tight_layout(rect=(0.08, 0.06, 0.92, 0.92))

    # Raster colour bar (even in SVG/PDF), sized to match the heatmap grid extent.
    user_ticks = _parse_csv_floats(getattr(args, "heatmap_cbar_ticks", None))
    tick_step = getattr(args, "heatmap_cbar_tick_step", None)
    ticks_use: list[float] | None = None
    tick_labels_use: list[str] | None = None
    if user_ticks:
        ticks_use = user_ticks
    elif use_boundary and ticks_boundary is not None:
        ticks_use = ticks_boundary
        tick_labels_use = tick_labels_boundary
    elif tick_step is not None:
        bticks = _cbar_ticks_from_step(float(norm.vmin), float(norm.vmax), float(tick_step))
        if bticks:
            ticks_use = bticks
    else:
        tt = ticks_for_norm(norm, max_ticks=6)
        if tt:
            ticks_use = tt
    if bsa_robust_active and bsa_robust_cap is not None:
        # Force the cap value onto the colour bar (and dedupe nearby ticks) so the split
        # point is unambiguously labelled.
        cap_val = float(bsa_robust_cap)
        existing = list(ticks_use) if ticks_use else []
        merged: list[float] = []
        eps = 1e-6 * max(abs(cap_val), 1.0)
        for t in existing:
            if abs(float(t) - cap_val) > eps:
                merged.append(float(t))
        merged.append(cap_val)
        ticks_use = sorted(merged)
        tick_labels_use = None
    add_raster_cbar(
        ax,
        cmap_obj,
        norm,
        nrows=n_row,
        ncols=n_col,
        label=f"{lab_pre}{cbar_label}{lab_suf}",
        orientation=cb_orient,
        vertical_colorbar_on_left=bool(y_right and cb_orient == "vertical"),
        ticks=ticks_use,
        tick_labels=tick_labels_use,
        tick_labelsize=float(getattr(args, "heatmap_cbar_tick_fontsize", 10.0)),
        label_fontsize=float(getattr(args, "heatmap_cbar_label_fontsize", 12.0)),
        thickness_scale=float(getattr(args, "heatmap_cbar_thickness_scale", 1.65)),
    )

    if (
        bsa_robust_active
        and bsa_robust_cap is not None
        and not bool(getattr(args, "heatmap_bsa_no_outlier_note", False))
    ):
        from utils.foldkit_heatmap import format_outlier_note  # noqa: E402

        use_short = bool(getattr(args, "heatmap_short_labels", False))
        items: list[tuple[str, str, float]] = []
        seen: set[tuple[str, str]] = set()
        for r in rows:
            try:
                bsa = float(r.get("buried_surface_area_A2") or "")
            except (TypeError, ValueError):
                continue
            if not (np.isfinite(bsa) and bsa > bsa_robust_cap):
                continue
            sb = str(r.get("structure_basename") or "").strip()
            cp = str(r.get("chain_pair_canonical") or "").strip()
            if not sb or not cp:
                continue
            key = (sb, cp)
            if key in seen:
                continue
            seen.add(key)
            sb_label = _structure_axis_label(sb, short_fn=short_heatmap_label, use_short=use_short)
            items.append((sb_label, cp, bsa))
        note_text = format_outlier_note(items, cap=float(bsa_robust_cap), units="\u00c5\u00b2")
        if note_text:
            note_fs = max(6.0, float(getattr(args, "heatmap_axis_labelsize", 9.0)) * 0.85)
            n_note_lines = note_text.count("\n") + 1
            extra_bottom = 0.045 + 0.018 * n_note_lines
            cur_bottom = float(fig.subplotpars.bottom)
            if cur_bottom < extra_bottom:
                fig.subplots_adjust(bottom=min(0.45, cur_bottom + extra_bottom))
            fig.text(
                0.02,
                0.01,
                note_text,
                ha="left",
                va="bottom",
                fontsize=note_fs,
                family="monospace",
            )
    root, _ext = os.path.splitext(path)
    path_out = root + "." + fmt_out

    os.makedirs(os.path.dirname(os.path.abspath(path_out)) or ".", exist_ok=True)
    dpi = int(getattr(args, "heatmap_dpi", 150))
    if vector_cells and fmt_out in ("svg", "pdf"):
        savefig_vector(fig, path_out, fmt_out, dpi=dpi)
    else:
        fig.savefig(path_out, format=fmt_out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path_out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "From interface analyser text reports, write per-interface CSV (summary, BSA, contacts, EC …), "
            "metric matrices with canonical pair columns, and optional heatmaps."
        )
    )
    ap.add_argument(
        "reports",
        nargs="+",
        metavar="PATH",
        help="One or more analyser .txt paths (shell globs ok).",
    )
    ap.add_argument(
        "-d",
        "--output-dir",
        metavar="DIR",
        help="Write outputs under DIR with basename interface_analysis_*.",
    )
    ap.add_argument(
        "--prefix",
        metavar="STEM",
        help="Output file stem (no extension), e.g. ./out/run1 → ./out/run1_per_interface.csv …",
    )
    ap.add_argument(
        "--metrics",
        nargs="+",
        choices=list(METRIC_CHOICES),
        metavar="NAME",
        default=list(DEFAULT_METRICS),
        help=(
            "Matrices and heatmaps to write: bsa, ec, ec_density, contacts (default set); optional "
            "charged_opposite, charged_same for charge complementarity reports."
        ),
    )

    hg = ap.add_argument_group(
        "heatmap",
        "Figure options (same spirit as utils/foldkit_heatmap.py and ranking/rmsd_to_csv.py).",
    )
    try:
        from utils.foldkit_heatmap import add_generic_heatmap_args  # noqa: E402
    except ImportError:
        add_generic_heatmap_args = None
    hg.add_argument(
        "--heatmap-structures-x-axis",
        action="store_true",
        help="Transpose matrix plots/CSVs so structures are on x-axis and interfaces on y-axis.",
    )
    hg.add_argument(
        "--heatmap-x-order",
        default=None,
        metavar="A,B,C,…",
        help=(
            "Comma-separated order for the x-axis labels (partial lists allowed). "
            "Applies to structures when using --heatmap-structures-x-axis; otherwise applies to "
            "interface canonical pairs."
        ),
    )
    hg.add_argument(
        "--heatmap-vmin-by-metric",
        action="append",
        default=None,
        metavar="METRIC=V",
        help="Per-heatmap vmin override (repeatable), e.g. --heatmap-vmin-by-metric bsa=0.",
    )
    hg.add_argument(
        "--heatmap-vmax-by-metric",
        action="append",
        default=None,
        metavar="METRIC=V",
        help="Per-heatmap vmax override (repeatable), e.g. --heatmap-vmax-by-metric bsa=800.",
    )
    hg.add_argument(
        "--heatmap-scale-by-metric",
        action="append",
        default=None,
        metavar="METRIC=SCALE",
        help=(
            "Per-heatmap scale override (repeatable), e.g. --heatmap-scale-by-metric contacts=log1p. "
            "SCALE choices match --heatmap-scale."
        ),
    )
    hg.add_argument(
        "--heatmap-boundaries-by-metric",
        action="append",
        default=None,
        metavar="METRIC=V0,V1,…",
        help=(
            "Per-heatmap discrete colour buckets (repeatable): ascending edges V0,V1,… produce "
            "len−1 intervals; colour-bar labels show each range. Overrides global --heatmap-boundaries "
            "for that metric."
        ),
    )
    if add_generic_heatmap_args is not None:
        add_generic_heatmap_args(hg, prefix="heatmap-")
    hg.add_argument(
        "--heatmap-annotate-metrics",
        nargs="*",
        choices=list(METRIC_CHOICES),
        default=[],
        metavar="NAME",
        help=(
            "Draw one numeric label per cell on these heatmaps (repeatable metric names). "
            "Defaults to the colour value; use --heatmap-cell-text to show a different CSV column."
        ),
    )
    hg.add_argument(
        "--heatmap-cell-text",
        action="append",
        default=None,
        metavar="METRIC=FIELD",
        help=(
            "For heatmap METRIC, print FIELD in each cell (colour still follows that heatmap's metric). "
            "FIELD is a per-interface CSV column name, or 'same' for the colour values (default when omitted). "
            "Example: --heatmap-cell-text ec=ec_density_per_1000_A2 colours by EC (r) with EC per 1000 Å² labels. "
            "Repeat the flag for multiple metrics."
        ),
    )
    hg.add_argument(
        "--heatmap-annotate-fmt",
        default="{:.2f}",
        metavar="FMT",
        help="Format string for single-metric cell annotations (default: '{:.2f}').",
    )
    hg.add_argument(
        "--heatmap-annotate-fontsize",
        type=float,
        default=7.0,
        metavar="PT",
        help="Font size for annotated heatmap cells (default: 7).",
    )
    hg.add_argument(
        "--heatmap-bsa-positive-range",
        action="store_true",
        help=(
            "For BSA heatmaps only: set default vmin/vmax from positive finite cells "
            "(RMSD-style autoscale)."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-annotate-contacts",
        action="store_true",
        help=(
            "Convenience: for the BSA heatmap, annotate cells with contact_count values "
            "(colour still uses BSA). Equivalent to: --heatmap-annotate-metrics bsa "
            "--heatmap-cell-text bsa=contact_count."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-robust",
        action="store_true",
        help=(
            "BSA only: split the colour bar between the bulk and any dominant interface(s). A "
            "canonical chain pair is auto-flagged as an outlier when its minimum BSA across "
            "structures exceeds --heatmap-bsa-outlier-factor (default 3) times the median of "
            "the remaining cells; the cap is the rounded max of non-outlier cells. Cells below "
            "the cap use --heatmap-cmap; above use --heatmap-bsa-above-cmap (default Reds). "
            "Above-cap pairs are listed in a figure note."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-outlier-factor",
        type=float,
        default=3.0,
        metavar="K",
        help=(
            "Multiplier (> 1) for the auto outlier rule: a canonical chain pair is an outlier "
            "when its min positive BSA across structures > K x median(remaining cells). "
            "Default: 3.0."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-split-at",
        type=float,
        default=None,
        metavar="VALUE",
        help=(
            "BSA only: explicit colour-bar split value (Å²). Wins over auto detection. Use "
            "when the auto rule does not isolate the desired interface(s)."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-above-cmap",
        default="Reds",
        metavar="NAME",
        help="Matplotlib colour map for cells above the BSA cap (default: Reds).",
    )
    hg.add_argument(
        "--heatmap-bsa-no-outlier-note",
        action="store_true",
        help="Disable the figure-level note listing above-cap BSA pairs.",
    )
    args = ap.parse_args(argv)
    try:
        args.heatmap_cell_text_map = _heatmap_cell_text_map_from_cli(list(args.heatmap_cell_text or []))
    except ValueError as exc:
        ap.error(str(exc))
    try:
        args.heatmap_vmin_by_metric = _metric_float_map_from_cli(list(args.heatmap_vmin_by_metric or []))
        args.heatmap_vmax_by_metric = _metric_float_map_from_cli(list(args.heatmap_vmax_by_metric or []))
        args.heatmap_scale_by_metric = _metric_str_map_from_cli(
            list(args.heatmap_scale_by_metric or []),
            allowed={"linear", "log10", "log1p", "clip_p95", "clip_p98"},
        )
        args.heatmap_boundaries_by_metric = _metric_edges_map_from_cli(
            list(args.heatmap_boundaries_by_metric or [])
        )
    except ValueError as exc:
        ap.error(str(exc))

    from utils.foldkit_heatmap import parse_boundaries_csv as _parse_boundaries_csv_main  # noqa: E402

    _gb_edges = _parse_boundaries_csv_main(getattr(args, "heatmap_boundaries", None))
    if len(_gb_edges) == 1:
        ap.error("--heatmap-boundaries needs at least two distinct comma-separated edge values.")

    if bool(getattr(args, "heatmap_bsa_annotate_contacts", False)):
        # Make sure BSA is in the annotation set.
        ann = list(getattr(args, "heatmap_annotate_metrics", None) or [])
        if "bsa" not in ann:
            ann.append("bsa")
        args.heatmap_annotate_metrics = ann
        # Only set default if user didn't explicitly define bsa cell-text.
        m = dict(getattr(args, "heatmap_cell_text_map", None) or {})
        m.setdefault("bsa", "contact_count")
        args.heatmap_cell_text_map = m

    paths = _expand_inputs(list(args.reports))
    if not paths:
        ap.error("No report files found.")

    if args.prefix and args.output_dir:
        ap.error("Use either --prefix or --output-dir, not both.")
    if args.prefix:
        prefix = os.path.abspath(args.prefix)
    elif args.output_dir:
        od = os.path.abspath(args.output_dir)
        os.makedirs(od, exist_ok=True)
        prefix = os.path.join(od, "interface_analysis")
    else:
        prefix = os.path.join(os.getcwd(), "interface_analysis")

    metrics_sel = list(dict.fromkeys(args.metrics))  # preserve order, dedupe
    transpose = bool(getattr(args, "heatmap_structures_x_axis", False))

    rows, _ = _collect_rows(paths)
    rows = _dedupe_iface_per_structure(rows)

    structures = sorted(
        {str(r["structure_basename"]) for r in rows if str(r.get("structure_basename") or "").strip()}
    )
    if not structures:
        ap.error("No interface rows parsed (check report format: EC or charge interface analyser output).")

    ordered_pairs = _ordered_canonical_pairs(rows)

    # Optional user-defined x-axis ordering.
    wanted_x = _parse_csv_list(getattr(args, "heatmap_x_order", None))
    try:
        if wanted_x:
            if bool(getattr(args, "heatmap_structures_x_axis", False)):
                # Structures are on x-axis.
                if bool(getattr(args, "heatmap_short_labels", False)):
                    structures = _apply_custom_order_by_key(
                        structures,
                        wanted_x,
                        key_fn=lambda sb: _structure_axis_label(sb, short_fn=lambda x: x, use_short=True),
                    )
                else:
                    structures = _apply_custom_order(structures, wanted_x)
            else:
                # Interfaces are on x-axis.
                ordered_pairs = _apply_custom_order(ordered_pairs, wanted_x)
    except ValueError as exc:
        ap.error(str(exc))

    p_iface = f"{prefix}_per_interface.csv"
    _write_per_interface_csv(p_iface, rows)

    written: list[str] = [p_iface]

    value_maps = {k: _value_by_struct_pair(rows, _METRIC_FIELD[k]) for k in METRIC_CHOICES}

    for key in metrics_sel:
        if key not in METRIC_CHOICES:
            continue
        stem_suffix = _METRIC_STEM_SUFFIX[key]
        p_mat = f"{prefix}_matrix_{stem_suffix}.csv"
        p_hm_base = f"{prefix}_heatmap_{stem_suffix}"
        _write_matrix_csv(
            p_mat,
            structures,
            ordered_pairs,
            value_maps[key],
            transpose=transpose,
        )
        written.append(p_mat)
        out_hm = _write_heatmap(
            p_hm_base,
            structures,
            ordered_pairs,
            value_maps[key],
            rows=rows,
            metric_key=key,
            args=args,
            cbar_label=_METRIC_HEATMAP_LABEL[key],
            title=_METRIC_TITLE[key],
            transpose=transpose,
        )
        if out_hm:
            written.append(out_hm)

    print(
        "Wrote:\n  "
        + "\n  ".join(written)
        + f"\nStructures: {len(structures)}, canonical interfaces: {len(ordered_pairs)}.\n"
        f"Matrices / heatmaps: {', '.join(metrics_sel)}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
