#!/usr/bin/env python3
"""
Join lattice packing JSON outputs with lattice EC text reports (matched by structure stem).

Writes:
  - combined_lattice_vs_ec.csv (scaled packing scalars + reference-chain EC summary)
  - lattice_compare_matrix_<metric>.csv + heatmaps (reuse interface_analysis_matrix._write_heatmap)
  - per-structure {stem}_full.csv (packing + EC summary + interfaces in one table, ``section`` column)

Heatmap defaults for this script: structures on the x-axis, horizontal colour bar, short structure
labels. Override with ``--no-heatmap-structures-x-axis``, ``--heatmap-colorbar-orientation vertical``,
and ``--no-heatmap-short-labels`` as needed.

Row order for ``combined_lattice_vs_ec.csv``: default is alphabetical by structure stem; use
``--structure-order`` or ``--structure-order-file``. Matrix/heatmap column order: ``--heatmap-x-order``.

Requires PYTHONPATH pointing at the FoldKit repo root (same as other metrics scripts).
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

from metrics.interface_analysis_matrix import (  # noqa: E402
    METRIC_CHOICES,
    PER_IFACE_COLUMNS,
    SUMMARY_KEYS,
    DEFAULT_METRICS,
    _METRIC_FIELD,
    _METRIC_HEATMAP_LABEL,
    _METRIC_STEM_SUFFIX,
    _METRIC_TITLE,
    _apply_custom_order,
    _apply_custom_order_by_key,
    _dedupe_iface_per_structure,
    _parse_csv_list,
    _structure_axis_label,
    _value_by_struct_pair,
    _write_heatmap,
    _write_matrix_csv,
)
from metrics.interface_matrix_heatmap_cli import (  # noqa: E402
    add_interface_matrix_heatmap_argument_group,
    apply_default_lattice_compare_heatmap_annotations,
    finalize_interface_matrix_heatmap_args,
)
from metrics.interface_mol_report_ec_csv import (  # noqa: E402
    _CSV_DETAILS_FIELDS,
    parse_ec_report_text,
)
from metrics.lattice_packing_report_csv import (  # noqa: E402
    _discover_from_dirs,
    _ensure_percent_fields,
    _iter_paths,
    _load_json_records,
)


def _stem_from_any(path_or_name: str) -> str:
    return os.path.splitext(os.path.basename(path_or_name.strip()))[0]


def _normalize_packing_stem(stem: str) -> str:
    """Map packing artefact names (packing_<stem>.json, <stem>_lattice_packing.json) to structure stem."""
    s = str(stem or "").strip()
    if s.startswith("packing_"):
        s = s[len("packing_") :]
    if s.endswith("_lattice_packing"):
        s = s[: -len("_lattice_packing")]
    return s


def _packing_row_stem(row: dict[str, Any]) -> str:
    candidates: list[str] = []
    if row.get("structure_stem"):
        candidates.append(str(row["structure_stem"]))
    if row.get("structure_basename"):
        candidates.append(_stem_from_any(str(row["structure_basename"])))
    inp = row.get("input") or ""
    if inp:
        candidates.append(_stem_from_any(str(inp)))
    src = row.get("source_json") or ""
    if src:
        candidates.append(_stem_from_any(str(src)))
    for raw in candidates:
        stem = _normalize_packing_stem(raw)
        if stem:
            return stem
    return ""


def _load_structure_order(order_csv: str | None, order_file: str | None) -> list[str]:
    """Stems from comma-separated CLI and/or one stem per line in a file (# comments allowed)."""
    wanted: list[str] = []
    if order_csv:
        wanted.extend(_parse_csv_list(order_csv))
    if order_file:
        path = os.path.expanduser(str(order_file).strip())
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Structure order file not found: {path}")
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                raw = line.strip()
                if not raw or raw.startswith("#"):
                    continue
                wanted.append(raw.split(",")[0].strip())
    # De-duplicate while preserving first occurrence.
    seen: set[str] = set()
    uniq: list[str] = []
    for w in wanted:
        if w and w not in seen:
            seen.add(w)
            uniq.append(w)
    return uniq


def _expand_ec_reports(inputs: list[str]) -> list[str]:
    out: list[str] = []
    for raw in inputs:
        raw = os.path.expanduser(raw.strip())
        if not raw:
            continue
        if any(ch in raw for ch in "*?["):
            matches = sorted(glob.glob(raw))
            if matches:
                out.extend(matches)
            else:
                print(f"Warning: skipping unmatched glob {raw!r}", file=sys.stderr)
        elif os.path.isfile(raw):
            out.append(raw)
        else:
            print(f"Warning: skipping missing path {raw!r}", file=sys.stderr)
    seen: set[str] = set()
    uniq: list[str] = []
    for p in out:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            uniq.append(ap)
    return uniq


def _collect_ec_rows_and_summaries(
    paths: list[str],
    *,
    reference_chain: str,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], dict[str, str]]:
    """
    Parse EC reports; build synthetic per-interface rows with reference-first chain_pair_canonical.
    Returns (iface_rows_for_heatmap, summary_by_stem, stem_to_report_basename).
    """
    ref = reference_chain.strip()
    iface_accum: list[dict[str, Any]] = []
    summary_by_stem: dict[str, dict[str, Any]] = {}
    stem_report: dict[str, str] = {}

    for path in paths:
        rtxt = os.path.basename(path)
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
        ifaces, sums = parse_ec_report_text(text)

        # Latest summary wins per structure basename stem.
        for s in sums:
            sb = str(s.get("structure_basename") or "").strip()
            if not sb:
                continue
            stem = _stem_from_any(sb)
            summary_by_stem[stem] = dict(s)
            stem_report[stem] = rtxt

        for iface in ifaces:
            sb = str(iface.get("structure_basename") or "").strip()
            if not sb:
                continue
            c1 = str(iface.get("chain1_id") or "").strip()
            c2 = str(iface.get("chain2_id") or "").strip()
            if ref not in (c1, c2):
                continue
            partner = c2 if c1 == ref else c1
            if not partner:
                continue
            iface_col = f"{ref}-{partner}"

            row = dict(iface)
            row["report_txt"] = rtxt
            row["structure_basename"] = sb
            row["chain_pair"] = iface_col
            row["chain_pair_canonical"] = iface_col
            row["ec_density_per_1000_A2"] = _ec_density_per_1000_A2_val(row.get("ec_density"))
            stem_i = _stem_from_any(sb)
            sum_rec = summary_by_stem.get(stem_i, {})
            for sk in SUMMARY_KEYS:
                row[sk] = sum_rec.get(sk, "")
            iface_accum.append(row)

    iface_accum = _dedupe_iface_per_structure(iface_accum)
    return iface_accum, summary_by_stem, stem_report


def _ordered_ref_pairs(rows: list[dict[str, Any]]) -> list[str]:
    pairs_set = {
        str(r["chain_pair_canonical"]) for r in rows if str(r.get("chain_pair_canonical") or "").strip()
    }
    return sorted(pairs_set)


COMBINED_EC_KEYS_REF = (
    "total_interfaces",
    "total_buried_surface_area_A2",
    "average_buried_area_per_interface_A2",
    "reference_chain_BSA_A2",
    "reference_chain_BSA_per_residue_reference_chain_A2",
    "reference_chain_BSA_per_atom_reference_chain_A2",
    "reference_chain_BSA_per_kDa_reference_chain_A2",
    "lattice_burial_fraction_percent",
    "lattice_contact_residue_fraction_percent",
    "lattice_ec_r_bsa_weighted",
    "lattice_ec_r_npairs_weighted",
)

COMBINED_ALL_KEYS_EXTRA = (
    "all_total_interfaces",
    "all_total_buried_surface_area_A2",
    "all_average_buried_area_per_interface_A2",
)

PACKING_TAIL_KEYS = (
    "volume_volume_a3",
    "volume_source",
    "lattice_reference_chain",
    "atom_density_per_1000_A3",
    "mass_density_Da_per_1000_A3",
    "packing_density_percent",
    "matthews_a3_per_Da",
    "estimated_solvent_percent",
)

# Omitted from per-stem merged CSV (lattice EC does not populate these in practice).
_MERGED_EC_OMIT_COLS: frozenset[str] = frozenset(
    {
        "distance_min_A",
        "distance_max_A",
        "distance_avg_A",
        "polarity_1_chain",
        "polarity_1_charged",
        "polarity_1_polar",
        "polarity_1_apolar",
        "polarity_1_other",
        "polarity_2_chain",
        "polarity_2_charged",
        "polarity_2_polar",
        "polarity_2_apolar",
        "polarity_2_other",
        "accessibility_1_chain",
        "accessibility_1_avg_sasa_A2",
        "accessibility_1_accessible_fraction",
        "accessibility_2_chain",
        "accessibility_2_avg_sasa_A2",
        "accessibility_2_accessible_fraction",
    }
)

_MERGED_EC_FIELDS_ORDER: tuple[str, ...] = tuple(
    f for f in _CSV_DETAILS_FIELDS if f not in _MERGED_EC_OMIT_COLS
)


def combined_csv_fieldnames(*, include_all_chains: bool) -> list[str]:
    """Fixed column order for combined_lattice_vs_ec.csv."""
    names: list[str] = ["structure_stem", *COMBINED_EC_KEYS_REF, "lattice_ec_density_bsa_weighted_per_1000_A2"]
    if include_all_chains:
        names.extend(COMBINED_ALL_KEYS_EXTRA)
    names.extend(PACKING_TAIL_KEYS)
    return names


def _ec_density_per_1000_A2_val(raw: Any) -> Any:
    if raw is None or raw == "":
        return ""
    try:
        return round(float(raw) * 1000.0, 2)
    except (TypeError, ValueError):
        return ""


def _build_combined_row(
    *,
    stem: str,
    flat_packing: dict[str, Any],
    ec_summary: dict[str, Any],
    include_all_chains: bool,
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "structure_stem": stem,
        "volume_volume_a3": flat_packing.get("volume_volume_a3", ""),
        "volume_source": flat_packing.get("volume_source", ""),
        "atom_density_per_1000_A3": "",
        "mass_density_Da_per_1000_A3": "",
        "packing_density_percent": flat_packing.get("lattice_packing_density_percent", ""),
        "matthews_a3_per_Da": flat_packing.get("lattice_matthews_a3_per_da", ""),
        "estimated_solvent_percent": flat_packing.get("estimated_solvent_content_percent", ""),
        "lattice_reference_chain": ec_summary.get("lattice_reference_chain", ""),
    }
    ad = flat_packing.get("lattice_atom_density_atoms_per_a3")
    md = flat_packing.get("lattice_mass_density_da_per_a3")
    try:
        if ad not in (None, ""):
            out["atom_density_per_1000_A3"] = float(ad) * 1000.0
    except (TypeError, ValueError):
        pass
    try:
        if md not in (None, ""):
            out["mass_density_Da_per_1000_A3"] = float(md) * 1000.0
    except (TypeError, ValueError):
        pass

    dens = ec_summary.get("lattice_ec_density_bsa_weighted")
    try:
        if dens not in (None, ""):
            out["lattice_ec_density_bsa_weighted_per_1000_A2"] = float(dens) * 1000.0
        else:
            out["lattice_ec_density_bsa_weighted_per_1000_A2"] = ""
    except (TypeError, ValueError):
        out["lattice_ec_density_bsa_weighted_per_1000_A2"] = ""

    for k in COMBINED_EC_KEYS_REF:
        if k == "reference_chain_BSA_A2":
            v = ec_summary.get("reference_chain_BSA_A2")
            if v in ("", None):
                iso = ec_summary.get("sasa_reference_isolated_A2")
                clu = ec_summary.get("sasa_reference_in_cluster_A2")
                try:
                    if iso not in ("", None) and clu not in ("", None):
                        v = float(iso) - float(clu)
                except (TypeError, ValueError):
                    v = ""
            out[k] = v
        else:
            out[k] = ec_summary.get(k, "")

    if include_all_chains:
        for k in COMBINED_ALL_KEYS_EXTRA:
            out[k] = ec_summary.get(k, "")
    return out


def _filter_packing_for_export(flat: dict[str, Any]) -> dict[str, Any]:
    """Drop metric_descriptions_* keys and warning from per-stem packing section."""
    out: dict[str, Any] = {}
    for k, v in flat.items():
        if k == "warning" or str(k).startswith("metric_descriptions_"):
            continue
        out[k] = v
    return out


def _cell_nonempty(v: Any) -> bool:
    if v is None:
        return False
    if isinstance(v, (int, float)):
        return True
    return str(v).strip() != ""


def _write_per_stem_merged_full_csv(
    path: str,
    *,
    stem: str,
    flat_packing: dict[str, Any],
    summary_rows: list[dict[str, Any]],
    iface_rows: list[dict[str, Any]],
    report_txt: str,
) -> None:
    """One CSV: packing row, EC summary row(s), interface rows; prune all-empty columns."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    rows_out: list[dict[str, Any]] = []

    prow: dict[str, Any] = {"section": "packing", "structure_stem": stem}
    prow.update(_filter_packing_for_export(flat_packing))
    rows_out.append(prow)

    for s in summary_rows:
        summ = dict(s)
        summ["section"] = "ec_summary"
        summ["structure_stem"] = stem
        summ.setdefault("report_txt", report_txt)
        rows_out.append(summ)

    for r in iface_rows:
        ir = dict(r)
        ir["section"] = "ec_interface"
        ir["structure_stem"] = stem
        ir.setdefault("report_txt", report_txt)
        rows_out.append(ir)

    always_keep: frozenset[str] = frozenset({"section", "structure_stem"})
    cand: list[str] = ["section", "structure_stem"]
    idents = [
        "structure_basename",
        "report_txt",
        "set_label",
        "chain_pair",
        "interface_number",
        "chain_pair_canonical",
        "chain_pair_raw",
        "chain1_id",
        "chain2_id",
    ]
    u: set[str] = set()
    for r in rows_out:
        u.update(r.keys())
    for x in idents:
        if x in u and x not in cand:
            cand.append(x)
    pack_keys = sorted(k for k in prow if k not in cand and k != "section")
    cand.extend(pack_keys)
    for f in _MERGED_EC_FIELDS_ORDER:
        if f in u and f not in cand:
            cand.append(f)
    for k in sorted(u):
        if k not in cand:
            cand.append(k)

    pruned: list[str] = []
    for c in cand:
        if c in always_keep:
            pruned.append(c)
            continue
        if any(_cell_nonempty(r.get(c)) for r in rows_out):
            pruned.append(c)

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=pruned, extrasaction="ignore")
        w.writeheader()
        for row in rows_out:
            w.writerow({k: row.get(k, "") for k in pruned})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Merge lattice packing JSON with lattice EC text reports (stem-matched); "
            "combined CSV, metric matrices, heatmaps, per-structure {stem}_full.csv. "
            "Heatmaps default to structures on x, horizontal colour bar, short structure labels "
            "(see module docstring)."
        ),
    )
    ap.add_argument(
        "--packing",
        nargs="+",
        metavar="PATH",
        required=True,
        help="Lattice packing JSON paths, dirs, or globs (cf. lattice_packing_report_csv.py).",
    )
    ap.add_argument(
        "--ec-reports",
        nargs="+",
        metavar="PATH",
        required=True,
        help="interface_analyser_lattice_ec.py text outputs (paths or globs).",
    )
    ap.add_argument(
        "--reference-chain",
        required=True,
        metavar="ID",
        help="Reference chain ID (columns ref-partner, filter interfaces involving this chain).",
    )
    ap.add_argument(
        "-d",
        "--output-dir",
        metavar="DIR",
        required=True,
        help="Output directory for all artefacts.",
    )
    ap.add_argument(
        "--json-glob",
        action="append",
        dest="json_globs",
        metavar="GLOB",
        help="Extra glob(s) when scanning dirs for packing JSON (repeatable). "
        "Defaults include *_lattice_packing.json and packing_*.json.",
    )
    ap.add_argument(
        "--include-all-chains-summary",
        action="store_true",
        help="Add all_total_interfaces (etc.) columns from all-chains summary block.",
    )
    ap.add_argument(
        "--strict-match",
        action="store_true",
        help="Exit with error if any packing stem lacks EC summary or vice versa.",
    )
    ap.add_argument(
        "--structure-order",
        default=None,
        metavar="STEM,STEM,…",
        help=(
            "Row order for combined_lattice_vs_ec.csv and per-structure {stem}_full.csv "
            "(structure stems, not report filenames). Partial lists allowed; "
            "any unmatched stems follow in alphabetical order."
        ),
    )
    ap.add_argument(
        "--structure-order-file",
        default=None,
        metavar="PATH",
        help="One structure stem per line (after --structure-order, if both are set).",
    )
    ap.add_argument(
        "--skip-heatmaps",
        action="store_true",
        help="Only write CSV matrices, not figure files.",
    )
    ap.add_argument(
        "--no-default-heatmap-annotations",
        action="store_true",
        help="Do not force BSA→contact_count and EC→ec_density_per_1000_A2 cell labels.",
    )
    ap.add_argument(
        "--metrics",
        nargs="+",
        choices=list(METRIC_CHOICES),
        metavar="NAME",
        default=list(DEFAULT_METRICS),
        help="Matrices / heatmaps to write (default: bsa ec ec_density contacts).",
    )

    add_interface_matrix_heatmap_argument_group(
        ap,
        METRIC_CHOICES,
        default_structures_x_axis=True,
        default_colorbar_orientation="horizontal",
        default_short_labels=True,
    )
    args = ap.parse_args(argv)

    finalize_interface_matrix_heatmap_args(
        args,
        ap,
        metric_choices=METRIC_CHOICES,
        cell_text_allowed=frozenset((*PER_IFACE_COLUMNS, "same")),
    )

    if not args.no_default_heatmap_annotations:
        apply_default_lattice_compare_heatmap_annotations(args)

    od = os.path.abspath(args.output_dir)
    os.makedirs(od, exist_ok=True)
    prefix = os.path.join(od, "lattice_compare")

    json_patterns = list(args.json_globs or []) or ["*_lattice_packing.json", "packing_*.json"]
    expanded = _iter_paths(list(args.packing))
    json_files, _txt_unused = _discover_from_dirs(expanded, json_patterns, [])
    if not json_files:
        ap.error("No lattice packing JSON files found (--packing / --json-glob).")

    packing_rows = _ensure_percent_fields(_load_json_records(json_files))
    packing_by_stem: dict[str, dict[str, Any]] = {}
    for pr in packing_rows:
        stem = _packing_row_stem(pr)
        if stem:
            packing_by_stem[stem] = pr

    ec_paths = _expand_ec_reports(list(args.ec_reports))
    if not ec_paths:
        ap.error("No EC report files found.")

    iface_rows, summary_by_stem, stem_report = _collect_ec_rows_and_summaries(
        ec_paths,
        reference_chain=args.reference_chain,
    )

    ec_stems = set(summary_by_stem.keys())
    pack_stems = set(packing_by_stem.keys())
    matched = sorted(pack_stems & ec_stems)
    only_pack = sorted(pack_stems - ec_stems)
    only_ec = sorted(ec_stems - pack_stems)

    if args.strict_match and (only_pack or only_ec):
        ap.error(
            "Stem mismatch under --strict-match: "
            f"only packing {only_pack}, only EC {only_ec}"
        )
    if only_pack:
        print(f"Warning: packing stems without EC summary (skipped): {only_pack}", file=sys.stderr)
    if only_ec:
        print(f"Warning: EC stems without packing JSON (skipped): {only_ec}", file=sys.stderr)

    if not matched:
        ap.error("No stems matched between packing JSON and EC reports.")

    try:
        wanted_stems = _load_structure_order(
            getattr(args, "structure_order", None),
            getattr(args, "structure_order_file", None),
        )
    except FileNotFoundError as exc:
        ap.error(str(exc))
    if wanted_stems:
        try:
            matched = _apply_custom_order(matched, wanted_stems)
        except ValueError as exc:
            ap.error(str(exc))

    combined_rows: list[dict[str, Any]] = []
    for stem in matched:
        flat_full = dict(packing_by_stem[stem])

        ec_sum = summary_by_stem[stem]
        combined_rows.append(
            _build_combined_row(
                stem=stem,
                flat_packing=flat_full,
                ec_summary=ec_sum,
                include_all_chains=bool(args.include_all_chains_summary),
            )
        )

        stem_ifaces = [r for r in iface_rows if _stem_from_any(str(r.get("structure_basename") or "")) == stem]
        summ_by_sb: dict[str, dict[str, Any]] = {}
        for path in ec_paths:
            with open(path, encoding="utf-8", errors="replace") as f:
                text = f.read()
            _, sums = parse_ec_report_text(text)
            for s in sums:
                if _stem_from_any(str(s.get("structure_basename") or "")) != stem:
                    continue
                sb = str(s.get("structure_basename") or "").strip()
                if sb:
                    summ_by_sb[sb] = dict(s)
        summ_rows_for_file = list(summ_by_sb.values())
        merged_path = os.path.join(od, f"{stem}_full.csv")
        _write_per_stem_merged_full_csv(
            merged_path,
            stem=stem,
            flat_packing=flat_full,
            summary_rows=summ_rows_for_file,
            iface_rows=stem_ifaces,
            report_txt=stem_report.get(stem, ""),
        )

    combined_path = os.path.join(od, "combined_lattice_vs_ec.csv")
    comb_fieldnames = combined_csv_fieldnames(include_all_chains=bool(args.include_all_chains_summary))
    with open(combined_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=comb_fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in combined_rows:
            w.writerow({k: r.get(k, "") for k in comb_fieldnames})

    heat_rows = [r for r in iface_rows if _stem_from_any(str(r.get("structure_basename") or "")) in set(matched)]
    structures = sorted(
        {str(r["structure_basename"]) for r in heat_rows if str(r.get("structure_basename") or "").strip()}
    )
    ordered_pairs = _ordered_ref_pairs(heat_rows)

    transpose = bool(getattr(args, "heatmap_structures_x_axis", False))
    wanted_x = _parse_csv_list(getattr(args, "heatmap_x_order", None))
    try:
        if wanted_x:
            if transpose:
                if bool(getattr(args, "heatmap_short_labels", False)):
                    structures = _apply_custom_order_by_key(
                        structures,
                        wanted_x,
                        key_fn=lambda sb: _structure_axis_label(sb, short_fn=lambda x: x, use_short=True),
                    )
                else:
                    structures = _apply_custom_order(structures, wanted_x)
            else:
                ordered_pairs = _apply_custom_order(ordered_pairs, wanted_x)
    except ValueError as exc:
        ap.error(str(exc))

    metrics_sel = list(dict.fromkeys(args.metrics))
    value_maps = {k: _value_by_struct_pair(heat_rows, _METRIC_FIELD[k]) for k in METRIC_CHOICES}

    written = [combined_path]
    for key in metrics_sel:
        if key not in METRIC_CHOICES:
            continue
        stem_suffix = _METRIC_STEM_SUFFIX[key]
        p_mat = f"{prefix}_matrix_{stem_suffix}.csv"
        _write_matrix_csv(
            p_mat,
            structures,
            ordered_pairs,
            value_maps[key],
            transpose=transpose,
        )
        written.append(p_mat)

        if args.skip_heatmaps:
            continue
        p_hm_base = f"{prefix}_heatmap_{stem_suffix}"
        title_use = _METRIC_TITLE[key].replace("Pairwise interface ", "Lattice interface ")
        out_hm = _write_heatmap(
            p_hm_base,
            structures,
            ordered_pairs,
            value_maps[key],
            rows=heat_rows,
            metric_key=key,
            args=args,
            cbar_label=_METRIC_HEATMAP_LABEL[key],
            title=title_use,
            transpose=transpose,
        )
        if out_hm:
            written.append(out_hm)

    extra_merged = [os.path.join(od, f"{s}_full.csv") for s in matched]
    print(
        "Wrote:\n  "
        + "\n  ".join(written + extra_merged)
        + f"\nMatched stems: {len(matched)}; structures in matrices: {len(structures)}; "
        f"reference-first interfaces: {len(ordered_pairs)}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
