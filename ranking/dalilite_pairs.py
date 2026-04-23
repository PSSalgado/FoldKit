#!/usr/bin/env python3
"""
Run DaliLite pairwise or all-vs-all, optionally write superposed PDBs (target onto query frame),
and emit Dali-style scores, Z-matrix, ranking, Newick tree, and an optional matplotlib Z-score heatmap
(via ``foldkit_heatmap``; this script does not import ``rmsd_to_csv``).

Superposition uses the translation–rotation block from DaliLite output (--outfmt transrot):
for target coordinates X, the superposition onto the query frame is R @ X + t
(DaliLite manual; same as applymatrix.pl).

The CSV ``z_score`` column is DaliLite’s summary Z when that line is parsed; otherwise FoldKit
fills it with the empirical (Holm-style) Z from ``dali_score.compute_z_score``; progress output
then labels "empirical Z-score". ``raw_score`` is always the recomputed Dali sum on the
mapped core. Pairs and Z-matrix are natural-sorted by label; the ranking table is ordered by
Z only. Optional ``--equivalences-dir`` writes one TSV of aligned residue pairs per comparison;
the pairs CSV can include min/max core residue numbers on each chain.

Requires: BioPython, NumPy, DaliLite (DALILITE_HOME or --dalilite-path). ``import.pl`` needs
``mkdssp`` at the path in DaliLite’s ``bin/mpidali.pm`` (``$DSSP_EXE``); see ``foldkit_dali_like_scores.py`` notes.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
import warnings

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import foldkit_dali_utils as du
import foldkit_dali_like_scores as _ds
from structure_phylogeny import _natural_sort_key

from cli_log import add_log_args, setup_log_from_args

warnings.filterwarnings("ignore", category=UserWarning, module="Bio")

try:
    from Bio.PDB.PDBIO import PDBIO
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    warnings.simplefilter("ignore", PDBConstructionWarning)
except ImportError:
    PDBIO = None


def _sanitize_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-._" else "_" for c in s)


def _find_dalilite_bin(dalilite_path: str | None) -> str:
    p = du._find_dalilite_path(dalilite_path)
    return p or ""


def parse_dalilite_transrot(content: str) -> tuple | None:
    """
    Parse first U(1,.), U(2,.), U(3,.) block from DaliLite text output.
    Returns (R 3x3 ndarray, t length-3 ndarray) or None.
    """
    if _ds.np is None:
        return None
    np = _ds.np
    rows: dict[int, list[float]] = {}
    for ln in content.splitlines():
        m = re.search(
            r"U\(([123]),\.\)\s+"
            r"([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)\s+"
            r"([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)\s+"
            r"([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)\s+"
            r"([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",
            ln,
        )
        if not m:
            continue
        idx = int(m.group(1))
        rows[idx] = [float(m.group(k)) for k in range(2, 6)]
    if not all(i in rows for i in (1, 2, 3)):
        return None
    R = np.array([rows[1][:3], rows[2][:3], rows[3][:3]], dtype=float)
    t = np.array([rows[1][3], rows[2][3], rows[3][3]], dtype=float)
    return R, t


def _kabsch_rotate_translate(P, Q) -> tuple:
    """
    Kabsch: align rows of Q (n×3) onto P (n×3). Return (R 3×3, t length-3) with
    P ≈ Q @ R.T + t (row vectors).
    """
    np = _ds.np
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    cP = P.mean(axis=0)
    cQ = Q.mean(axis=0)
    Pc = P - cP
    Qc = Q - cQ
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = cP - (R @ cQ)
    return R, t


def _pair_result_from_biotite_fallback(
    pdb_a: str,
    pdb_b: str,
    chain_a: str,
    chain_b: str,
    work_cwd: str,
    prior_dali_content: str = "",
    dali_pipeline_error: str = "",
) -> tuple[dict | None, str]:
    """
    When DaliLite omits hits (e.g. Z<~2), use biotite structural alignment + Kabsch
    on matched Cα. Z in output is empirical (dali_score), not Dali Z.
    """
    if not _ds.BIOTITE_AVAILABLE:
        return None, "biotite is not installed (pip install biotite)"
    equiv_raw = _ds.alignment_from_biotite(pdb_a, pdb_b, chain_a, chain_b)
    if not equiv_raw:
        return None, "biotite superimpose_structural_homologs found no correspondence"
    parser = _ds.PDBParser(QUIET=True)
    sa = parser.get_structure("A", pdb_a)
    sb = parser.get_structure("B", pdb_b)
    coords_A = _ds.get_ca_coords_biopython(sa, chain_a)
    coords_B = _ds.get_ca_coords_biopython(sb, chain_b)
    equivs = _ds.normalize_equivalences(equiv_raw, coords_A, coords_B)
    if len(equivs) < 3:
        return None, f"biotite fallback: only {len(equivs)} matched residues after normalisation"
    np = _ds.np
    P = np.stack([coords_A[ka] for ka, _ in equivs])
    Q = np.stack([coords_B[kb] for _, kb in equivs])
    R, tvec = _kabsch_rotate_translate(P, Q)
    Qal = Q @ R.T + tvec
    rmsd = float(np.sqrt(((P - Qal) ** 2).sum() / len(P)))
    raw_score, n_core = _ds.compute_dali_score(coords_A, coords_B, equivs)
    L_A, L_B = len(coords_A), len(coords_B)
    note = (
        "# biotite_fallback: DaliLite reported no significant hit (often Z<~2) but "
        "biotite superimpose_structural_homologs found a correspondence (like Coot SSM, "
        "this can disagree with Dali’s statistical cutoff).\n"
        f"# core_rmsd_after_kabsch_angstrom: {rmsd:.4f}\n"
        "# z_score in CSV is empirical (Holm-style fit from foldkit_dali_like_scores), not Dali Z.\n"
    )
    if prior_dali_content:
        note += "# --- DaliLite output (no equivalences) ---\n" + prior_dali_content
    if dali_pipeline_error:
        note += "# --- DaliLite pipeline message ---\n" + dali_pipeline_error + "\n"
    os.makedirs(work_cwd, exist_ok=True)
    out_file = os.path.join(work_cwd, "biotite_fallback.txt")
    try:
        with open(out_file, "w", encoding="utf-8") as tf:
            tf.write(note)
    except OSError:
        out_file = ""
    return (
        {
            "z_score": None,
            "rmsd": rmsd,
            "lali": n_core,
            "nres": None,
            "pct_id": None,
            "hit_id": None,
            "description": "",
            "equivs": equivs,
            "raw_text": note,
            "rotation": R,
            "translation": _ds.np.asarray(tvec, dtype=float).reshape(3),
            "dali_txt_path": out_file,
            "resolved_chain_a": chain_a,
            "resolved_chain_b": chain_b,
            "alignment_source": "biotite_fallback",
        },
        "",
    )


def _resolve_chain_ids(
    pdb_a: str,
    pdb_b: str,
    chain_a: str | None,
    chain_b: str | None,
) -> tuple[str, str]:
    """Match foldkit_dali_like_scores.run(): default to first chain with CA in each structure."""
    parser = _ds.PDBParser(QUIET=True)
    sa = parser.get_structure("A", pdb_a)
    sb = parser.get_structure("B", pdb_b)
    ca = _ds.get_ca_coords_biopython(sa, chain_a)
    cb = _ds.get_ca_coords_biopython(sb, chain_b)
    cid_a = chain_a or _ds._get_first_chain_id(ca)
    cid_b = chain_b or _ds._get_first_chain_id(cb)
    return cid_a, cid_b


def run_dalilite_pair(
    pdb_a: str,
    pdb_b: str,
    work_cwd: str,
    chain_a: str,
    chain_b: str,
    dalilite_path: str | None,
    debug_dalilite: bool = False,
    fallback_biotite: bool = False,
) -> tuple[dict | None, str]:
    """
    Run DaliLite import+dali (short temp-dir staging inside dali_score); save full text
    to work_cwd/{mol1}{chain}.txt when possible. Returns (result_dict_or_None, err_msg).
    If ``fallback_biotite`` and DaliLite omits a hit (common for Z<~2), try biotite
    alignment + Kabsch (Z in scores is then empirical, not Dali Z).
    """
    pdb_a = os.path.abspath(pdb_a)
    pdb_b = os.path.abspath(pdb_b)
    if not os.path.isfile(pdb_a) or not os.path.isfile(pdb_b):
        return None, "PDB path not a file"

    ca, cb = chain_a, chain_b
    pdbid1 = "mol1"
    os.makedirs(work_cwd, exist_ok=True)

    if not _find_dalilite_bin(dalilite_path):
        if fallback_biotite and _ds.BIOTITE_AVAILABLE:
            fb, fe = _pair_result_from_biotite_fallback(
                pdb_a,
                pdb_b,
                ca,
                cb,
                work_cwd,
                dali_pipeline_error="DaliLite not found (set DALILITE_HOME or --dalilite-path)",
            )
            if fb:
                return fb, ""
            return None, fe or "biotite-only fallback failed"
        return None, "DaliLite not found: set DALILITE_HOME or pass --dalilite-path"

    content, ierr, ra, rb = du._dalilite_pair_via_dat(
        pdb_a,
        pdb_b,
        ca,
        cb,
        dalilite_path,
        work_cwd,
        "summary,equivalences,alignments,transrot",
        "dalilite_pairs",
        dali_timeout=600,
    )
    if content is None:
        if debug_dalilite:
            print(f"[debug] DaliLite pipeline: {ierr}", file=sys.stderr)
        if fallback_biotite and _ds.BIOTITE_AVAILABLE:
            fb, fe = _pair_result_from_biotite_fallback(
                pdb_a,
                pdb_b,
                ca,
                cb,
                work_cwd,
                dali_pipeline_error=ierr or "",
            )
            if fb:
                return fb, ""
            return None, (ierr or "DaliLite failed") + f" | biotite fallback: {fe}"
        return None, ierr or "DaliLite import/compare failed"

    out_name = f"{pdbid1}{ra}.txt"
    out_file = os.path.join(work_cwd, out_name)
    try:
        os.makedirs(work_cwd, exist_ok=True)
        with open(out_file, "w", encoding="utf-8") as tf:
            tf.write(content)
    except OSError:
        out_file = ""

    hit = du.parse_dalilite_summary_hit(content)
    z_score = rmsd = lali = None
    nres = pct_id = None
    hit_id = ""
    description = ""
    if hit:
        z_score = hit["z_score"]
        rmsd = hit["rmsd"]
        lali = hit["lali"]
        nres = hit["nres"]
        pct_id = hit["pct_id"]
        hit_id = hit["hit_id"]
        description = hit["description"]

    equivs = du._parse_dalilite_equivalences(content, ra, rb)
    if not equivs:
        if debug_dalilite:
            print(
                f"[debug] equivs parse failed; output head:\n{content[:2000]}",
                file=sys.stderr,
            )
        if not du._dalilite_text_has_summary_hit(content):
            hint = (
                "DaliLite found no significant structural match: the summary table has no "
                "hit row, so equivalences and transrot are empty. This is typical when Z "
                "would be below DaliLite’s reporting threshold (~2)—e.g. unrelated folds, "
                "different domains, or very divergent models. Try pairs you expect to be "
                "structurally homologous, or use foldkit_dali_like_scores.py with --no-dalilite for "
                "biotite/sequence-order alignment (different method, not Dali Z). "
            )
        else:
            hint = (
                "DaliLite reported a summary hit but equivalences could not be parsed "
                "(unexpected v5 text format). "
            )
        if fallback_biotite and _ds.BIOTITE_AVAILABLE:
            fb, fe = _pair_result_from_biotite_fallback(
                pdb_a,
                pdb_b,
                ca,
                cb,
                work_cwd,
                prior_dali_content=content,
                dali_pipeline_error=hint.strip()[:800],
            )
            if fb:
                return fb, ""
            return None, hint + f" | biotite fallback failed: {fe} | head: {content[:800]!r}"
        return None, hint + f"Output head: {content[:1200]!r}"

    rt = parse_dalilite_transrot(content)
    R, tvec = (None, None) if rt is None else rt

    return (
        {
            "z_score": z_score,
            "rmsd": rmsd,
            "lali": lali,
            "nres": nres,
            "pct_id": pct_id,
            "hit_id": hit_id,
            "description": description,
            "equivs": equivs,
            "raw_text": content,
            "rotation": R,
            "translation": tvec,
            "dali_txt_path": out_file,
            "resolved_chain_a": ra,
            "resolved_chain_b": rb,
            "alignment_source": "dalilite",
        },
        "",
    )


def write_target_superposed_pdb(
    pdb_b: str,
    rotation,
    translation,
    out_path: str,
) -> bool:
    """Apply R @ x + t to all atoms in pdb_b; write PDB. DaliLite target -> query frame."""
    if PDBIO is None or _ds.PDBParser is None or _ds.np is None:
        return False
    if rotation is None or translation is None:
        return False
    R = _ds.np.asarray(rotation, dtype=float)
    t = _ds.np.asarray(translation, dtype=float).reshape(3)
    parser = _ds.PDBParser(QUIET=True)
    try:
        struct = parser.get_structure("b", pdb_b)
    except Exception:
        return False
    for atom in struct.get_atoms():
        c = _ds.np.array(atom.coord, dtype=float)
        atom.coord = R @ c + t
    try:
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        io = PDBIO()
        io.set_structure(struct)
        io.save(out_path)
    except OSError:
        return False
    return True


def write_fixed_copy(pdb_a: str, out_path: str) -> bool:
    try:
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        shutil.copy2(pdb_a, out_path)
        return True
    except OSError:
        return False


def score_pair_from_dalilite(
    pdb_a: str,
    pdb_b: str,
    chain_a: str | None,
    chain_b: str | None,
    dali: dict,
) -> dict:
    """Raw Dali-like score via dali_score; Z from DaliLite when present else empirical."""
    parser = _ds.PDBParser(QUIET=True)
    struct_a = parser.get_structure("A", pdb_a)
    struct_b = parser.get_structure("B", pdb_b)
    ca_use = dali.get("resolved_chain_a", chain_a)
    cb_use = dali.get("resolved_chain_b", chain_b)
    coords_A = _ds.get_ca_coords_biopython(struct_a, ca_use)
    coords_B = _ds.get_ca_coords_biopython(struct_b, cb_use)
    ca = ca_use or _ds._get_first_chain_id(coords_A)
    cb = cb_use or _ds._get_first_chain_id(coords_B)
    equivs = _ds.normalize_equivalences(dali["equivs"], coords_A, coords_B)
    raw_score, n_core = _ds.compute_dali_score(coords_A, coords_B, equivs)
    span = _ds.core_residue_spans_from_equivs(equivs)
    L_A, L_B = len(coords_A), len(coords_B)
    z_dali = dali.get("z_score")
    z_is_empirical = z_dali is None
    if z_dali is not None:
        z_out = float(z_dali)
    else:
        z_out = float(_ds.compute_z_score(raw_score, L_A, L_B))
    src = dali.get("alignment_source", "dalilite")
    return {
        "raw_score": raw_score,
        "z_score": z_out,
        "z_is_empirical": z_is_empirical,
        "n_core": n_core,
        "L_A": L_A,
        "L_B": L_B,
        "dalilite_rmsd": dali.get("rmsd") if src == "dalilite" else None,
        "core_rmsd": dali.get("rmsd"),
        "alignment_source": src,
        "lali": dali.get("lali"),
        "nres": dali.get("nres"),
        "pct_id": dali.get("pct_id"),
        "dalilite_hit_id": dali.get("hit_id"),
        "dali_description": dali.get("description"),
        "equivs": equivs,
        **span,
    }


def run_all_pairs(
    files: list[str],
    out_dir: str,
    chain_a: str | None,
    chain_b: str | None,
    dalilite_path: str | None,
    keep_work: bool,
    write_superpose: bool,
    verbose: bool,
    debug_dalilite: bool = False,
    fallback_biotite: bool = False,
    equivalences_dir: str | None = None,
) -> dict:
    # Match rmsd_to_csv / structure_phylogeny: axis order is natural-sorted structure labels.
    files = sorted(
        files,
        key=lambda p: (_natural_sort_key(du._label_from_path(p)), _natural_sort_key(p)),
    )
    labels = [du._label_from_path(f) for f in files]
    zscores: dict = {}
    raw_scores: dict = {}
    n_core: dict = {}
    results_rows: list[dict] = []
    tmp_root = os.path.join(out_dir, "_dalilite_work")
    super_dir = os.path.join(out_dir, "superposed")

    n_pairs = len(files) * (len(files) - 1) // 2
    done = 0

    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            pa, pb = files[i], files[j]
            la, lb = labels[i], labels[j]
            done += 1
            pair_tag = f"{_sanitize_label(la)}__{_sanitize_label(lb)}"
            work = os.path.join(tmp_root, pair_tag)
            if not keep_work and os.path.isdir(work):
                shutil.rmtree(work, ignore_errors=True)

            cid_a, cid_b = _resolve_chain_ids(pa, pb, chain_a, chain_b)
            dali, derr = run_dalilite_pair(
                pa,
                pb,
                work,
                cid_a,
                cid_b,
                dalilite_path,
                debug_dalilite=debug_dalilite,
                fallback_biotite=fallback_biotite,
            )
            if dali is None:
                if verbose:
                    print(
                        f"[{done}/{n_pairs}] {la} vs {lb}: DaliLite failed — {derr}",
                        file=sys.stderr,
                    )
                continue

            meta = score_pair_from_dalilite(pa, pb, chain_a, chain_b, dali)
            z_emp_flag = bool(meta.get("z_is_empirical"))
            equivs = meta.pop("equivs", [])
            meta.pop("z_is_empirical", None)
            if equivalences_dir and equivs:
                os.makedirs(equivalences_dir, exist_ok=True)
                du.write_equivalences_tsv(
                    os.path.join(equivalences_dir, f"{pair_tag}.tsv"), equivs
                )
            z = meta["z_score"]
            zscores[(la, lb)] = zscores[(lb, la)] = z
            raw_scores[(la, lb)] = raw_scores[(lb, la)] = meta["raw_score"]
            n_core[(la, lb)] = n_core[(lb, la)] = meta["n_core"]

            if verbose:
                tag = ""
                if meta.get("alignment_source") == "biotite_fallback":
                    tag = " [biotite fallback; core_rmsd ≈ Kabsch on matched Cα]"
                if z_emp_flag:
                    z_part = f"empirical Z-score={z:.2f} (FoldKit dali-like)"
                else:
                    z_part = f"Z={z:.2f} (DaliLite)"
                print(
                    f"[{done}/{n_pairs}] {la} vs {lb}: {z_part} "
                    f"raw={meta['raw_score']:.2f} n_core={meta['n_core']}{tag}",
                    file=sys.stderr,
                )

            results_rows.append(
                {
                    "pdb_a": pa,
                    "pdb_b": pb,
                    "label_a": la,
                    "label_b": lb,
                    **{k: meta[k] for k in ("raw_score", "z_score", "n_core")},
                    "core_resnum_min_a": meta.get("core_resnum_min_a", ""),
                    "core_resnum_max_a": meta.get("core_resnum_max_a", ""),
                    "core_resnum_min_b": meta.get("core_resnum_min_b", ""),
                    "core_resnum_max_b": meta.get("core_resnum_max_b", ""),
                    "dalilite_rmsd": meta.get("dalilite_rmsd"),
                    "core_rmsd": meta.get("core_rmsd"),
                    "alignment_source": meta.get("alignment_source", "dalilite"),
                    "lali": meta.get("lali"),
                    "nres": meta.get("nres"),
                    "pct_id": meta.get("pct_id"),
                    "dalilite_hit_id": meta.get("dalilite_hit_id"),
                    "dali_description": meta.get("dali_description"),
                }
            )

            if write_superpose and dali.get("rotation") is not None:
                os.makedirs(super_dir, exist_ok=True)
                fix_p = os.path.join(super_dir, f"{pair_tag}_query.pdb")
                mob_p = os.path.join(super_dir, f"{pair_tag}_target_on_query.pdb")
                if write_fixed_copy(pa, fix_p) and write_target_superposed_pdb(
                    pb, dali["rotation"], dali["translation"], mob_p
                ):
                    if verbose:
                        print(f"  wrote {fix_p} and {mob_p}", file=sys.stderr)
                elif verbose:
                    print("  superposed PDB write failed", file=sys.stderr)

            if not keep_work:
                shutil.rmtree(work, ignore_errors=True)

    return {
        "zscores": zscores,
        "raw_scores": raw_scores,
        "n_core": n_core,
        "labels": labels,
        "results_rows": results_rows,
    }


def _matrix_from_zscores(labels: list[str], zscores: dict) -> list[list[float | None]]:
    """Square matrix of Z-scores; diagonal None; missing pairs None."""
    mat: list[list[float | None]] = []
    for i, la in enumerate(labels):
        row: list[float | None] = []
        for j, lb in enumerate(labels):
            if i == j:
                row.append(None)
                continue
            z = zscores.get((la, lb))
            if z is None:
                z = zscores.get((lb, la))
            row.append(float(z) if z is not None else None)
        mat.append(row)
    return mat


def _matrix_from_ncore(labels: list[str], n_core: dict) -> list[list[int | float | None]]:
    """Square matrix of n_core; diagonal None; missing pairs None."""
    mat: list[list[int | float | None]] = []
    for i, la in enumerate(labels):
        row: list[int | float | None] = []
        for j, lb in enumerate(labels):
            if i == j:
                row.append(None)
                continue
            v = n_core.get((la, lb), n_core.get((lb, la)))
            row.append(int(v) if v is not None else None)
        mat.append(row)
    return mat


def _write_dalilite_z_heatmap(
    labels: list[str],
    zscores: dict,
    heatmap_path: str,
    *,
    title: str,
    cmap: str,
    vmin: float | None,
    vmax: float | None,
    short_labels: bool,
    diverging_center: str | None,
    colorbar_orientation: str,
    y_axis_right: bool,
    n_core: dict | None = None,
    n_core_heatmap: bool = False,
    n_core_n_equal_bins: int = 4,
    n_core_bin_edges: list[float] | None = None,
) -> None:
    """Square Dali Z heatmap via ``foldkit_heatmap.plot_heatmap`` (hatch legend: aligned Cα, Dali n_core)."""
    try:
        from foldkit_heatmap import plot_heatmap as _plot_heatmap
    except ImportError as e:
        print("Error: --heatmap requires foldkit_heatmap in the same package.", file=sys.stderr)
        raise SystemExit(1) from e
    mat = _matrix_from_zscores(labels, zscores)
    n_core_mat: list | None = None
    if n_core_heatmap and n_core and len(labels) >= 2:
        n_core_mat = _matrix_from_ncore(labels, n_core)
    kwargs: dict = {
        "cbar_label": "Dali Z",
        "autoscale_positive_offdiag_only": False,
    }
    if n_core_mat is not None:
        kwargs["hatch_value_matrix"] = n_core_mat
        kwargs["hatch_n_equal_bins"] = max(1, int(n_core_n_equal_bins))
        kwargs["hatch_channel_name"] = "Dali n_core (aligned Cα)"
        if n_core_bin_edges is not None:
            kwargs["hatch_bin_edges"] = n_core_bin_edges
    _plot_heatmap(
        labels,
        mat,
        os.path.abspath(heatmap_path),
        title=title,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        short_labels=short_labels,
        diverging_center=diverging_center or None,
        vcenter=None,
        colorbar_orientation=colorbar_orientation,
        y_axis_right=y_axis_right,
        **kwargs,
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Run DaliLite, write superposed PDBs (target on query frame via transrot), "
            "pairwise scores / Z-matrix / tree (via dali_score), and an optional matplotlib Z heatmap."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Pairwise:
    %(prog)s /path/to/project/model_01.pdb /path/to/project/model_02.pdb \\
      -d /path/to/project/out --dalilite-path /path/to/DaliLite

  All-vs-all (at least two paths; quote glob filters; directory roots only unless --recursive):
    %(prog)s --all-vs-all /path/to/project/models/ /path/to/project/models/ \\
      -d /path/to/project/run_01/ --filter '*.pdb' --dalilite-path /path/to/DaliLite

  All-vs-all including subdirectories:
    %(prog)s --all-vs-all /path/to/project/models/ /path/to/project/models/ \\
      -d /path/to/project/run_01/ --filter '*.pdb' --recursive --dalilite-path /path/to/DaliLite

  All-vs-all with Z-score heatmap (SVG) and median-centred colour scale:
    %(prog)s --all-vs-all /path/to/project/models/ /path/to/project/models/ \\
      -d /path/to/project/run_01/ --filter '*.pdb' --dalilite-path /path/to/DaliLite \\
      --heatmap /path/to/project/run_01/z_matrix.svg --heatmap-diverging-center median

  When DaliLite omits a hit (Z<~2) but models still superpose in Coot (e.g. SSM):
    %(prog)s ... --fallback-biotite   # needs: pip install biotite

Environment: DALILITE_HOME or --dalilite-path. DaliLite runs from a short path under
the system temp directory (very long paths can hit Fortran 80-character limits inside import.pl).

mkdssp: DaliLite's import.pl calls mkdssp via bin/mpidali.pm ($DSSP_EXE). That binary
must exist (symlink a compatible mkdssp to the expected path, or edit mpidali.pm).
Optional: MKDSSP=/path/to/mkdssp helps dali_score diagnostics match your install.
""",
    )
    ap.add_argument(
        "paths",
        nargs="*",
        help=(
            "Pairwise: two structure files. All-vs-all: directories and/or files; "
            "per-directory files default to the directory root only (see --recursive)."
        ),
    )
    ap.add_argument(
        "-d",
        "--output-dir",
        metavar="DIR",
        required=True,
        help="Output directory (scores, tree, superposed/).",
    )
    ap.add_argument(
        "--all-vs-all",
        action="store_true",
        help="All unordered pairs from collected PDB/CIF paths.",
    )
    ap.add_argument(
        "--recursive",
        action="store_true",
        help=(
            "With --all-vs-all: search each given directory recursively for "
            "*.pdb / *.cif / *.ent (default: only the root of each directory, not subfolders)."
        ),
    )
    ap.add_argument(
        "--filter",
        metavar="PATTERN",
        help="Filename filter (see foldkit_dali_like_scores.py).",
    )
    ap.add_argument("--chain-a", metavar="ID", help="Chain ID for structure A (first/query).")
    ap.add_argument("--chain-b", metavar="ID", help="Chain ID for structure B (second/target).")
    ap.add_argument(
        "--dalilite-path",
        metavar="DIR",
        help="DaliLite installation directory (else DALILITE_HOME); mkdssp must satisfy bin/mpidali.pm.",
    )
    ap.add_argument(
        "--no-superpose-pdb",
        action="store_true",
        help="Skip writing superposed PDB files.",
    )
    ap.add_argument(
        "--keep-dalilite-work",
        action="store_true",
        help="Keep per-pair DaliLite work dirs under _dalilite_work/.",
    )
    ap.add_argument(
        "--equivalences-dir",
        metavar="DIR",
        default=None,
        help="Optional: write one TSV per pair listing aligned Cα pairs (structural core used for n_core and raw score).",
    )
    ap.add_argument(
        "--fallback-biotite",
        action="store_true",
        help=(
            "If DaliLite reports no hit (often Z<~2), align with biotite + Kabsch; "
            "Z in CSV is empirical (not Dali Z). Requires biotite."
        ),
    )
    ap.add_argument("-o", "--output", metavar="FILE", help="Pairwise results CSV basename.")
    ap.add_argument(
        "--output-tree",
        metavar="FILE",
        default="dali_tree.nwk",
        help="Newick tree (default: dali_tree.nwk in output directory).",
    )
    ap.add_argument("--output-plot", metavar="FILE", help="Dendrogram image path.")
    ap.add_argument("--output-matrix", metavar="FILE", help="Z-score matrix CSV path.")
    ap.add_argument(
        "--output-ranking",
        metavar="FILE",
        default="dali_ranking.csv",
        help="Ranking CSV (default: dali_ranking.csv in output directory).",
    )
    ap.add_argument(
        "--transform",
        choices=["inv", "maxminus", "exp"],
        default="inv",
        help="Z to distance for NJ tree (default: inv).",
    )
    ap.add_argument(
        "--exp-scale",
        type=float,
        default=10.0,
        metavar="SCALE",
        help="Scale for exp transform: d=exp(-Z/scale).",
    )
    ap.add_argument("--root", metavar="LABEL", help="Outgroup label for rooting.")
    ap.add_argument(
        "--no-midpoint-root",
        action="store_true",
        help="Disable midpoint rooting when no --root.",
    )
    ap.add_argument(
        "--output-tree-nomid",
        metavar="FILE",
        default=None,
        help=(
            "All-vs-all: optional second Newick (NJ on Z) built without midpoint-rooting, "
            "i.e. --no-midpoint-root (ignored when the primary tree is already un-midpoint-rooted). "
            "Use with --output-plot-nomid to also render a second dendrogram."
        ),
    )
    ap.add_argument(
        "--output-plot-nomid",
        metavar="FILE",
        default=None,
        help="All-vs-all: optional second dendrogram image for --output-tree-nomid (same renderer as --output-plot).",
    )
    ap.add_argument(
        "--heatmap",
        metavar="PATH",
        default=None,
        help=(
            "Write a pairwise Dali Z-score heatmap (PNG, PDF, or SVG from the extension). "
            "All-vs-all: full matrix; pairwise: 2×2. Heatmap via foldkit_heatmap (same options as rmsd_to_csv); "
            "requires matplotlib."
        ),
    )
    ap.add_argument(
        "--heatmap-title",
        metavar="STR",
        default="DaliLite pairwise Z-scores",
        help="Figure title for --heatmap (default: DaliLite pairwise Z-scores).",
    )
    ap.add_argument(
        "--cmap",
        default="viridis_r",
        help="Matplotlib colour map for --heatmap (default: viridis_r).",
    )
    ap.add_argument(
        "--vmin",
        type=float,
        default=None,
        metavar="VAL",
        help="Lower colour scale limit for --heatmap (Dali Z).",
    )
    ap.add_argument(
        "--vmax",
        type=float,
        default=None,
        metavar="VAL",
        help="Upper colour scale limit for --heatmap (Dali Z).",
    )
    ap.add_argument(
        "--short-heatmap-labels",
        action="store_true",
        help="Shorten --heatmap axis labels (same convention as ranking/rmsd_to_csv.py).",
    )
    ap.add_argument(
        "--heatmap-diverging-center",
        choices=("none", "median"),
        default="none",
        help="For --heatmap: linear scale (none) or TwoSlopeNorm at the median Z (median).",
    )
    ap.add_argument(
        "--heatmap-colorbar-orientation",
        choices=("vertical", "horizontal"),
        default="vertical",
        help="Colour bar placement for --heatmap (default: vertical).",
    )
    ap.add_argument(
        "--heatmap-y-axis-right",
        action="store_true",
        help="Draw row labels on the right for --heatmap (vertical bar moves to the left).",
    )
    ap.add_argument(
        "--heatmap-n-core-patterns",
        action="store_true",
        help=(
            "On the Z heatmap, encode n_core in hatch patterns; colour still encodes Dali Z. "
            "Uses equal-width n_core ranges (--heatmap-n-core-bins) or explicit "
            "boundaries (--heatmap-n-core-edges)."
        ),
    )
    ap.add_argument(
        "--heatmap-n-core-bins",
        type=int,
        default=4,
        metavar="N",
        help="With --heatmap-n-core-patterns: number of equal-width n_core ranges (default: 4).",
    )
    ap.add_argument(
        "--heatmap-n-core-edges",
        metavar="E0,E1,…",
        default=None,
        help=(
            "With --heatmap-n-core-patterns: optional comma-separated bin boundaries, low to high; "
            "they are extended to the observed n_core min and max on the heatmap if needed."
        ),
    )
    ap.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Reduce console output.",
    )
    ap.add_argument(
        "--debug-dalilite",
        action="store_true",
        help="Print dali.pl stderr when a pair fails (nonzero exit, missing .txt, or no equivalences).",
    )
    add_log_args(ap)
    args = ap.parse_args()
    summary_log = setup_log_from_args(
        args, script_path=__file__, inputs=[getattr(args, "input", "")], pattern=None
    )
    if summary_log is not None:
        summary_log.task("DaliLite superpose + scores (dalilite_pairs.py)")

    if not _ds.BIOPYTHON_AVAILABLE or _ds.PDBParser is None:
        print("Error: BioPython is required.", file=sys.stderr)
        sys.exit(1)
    if _ds.np is None:
        print("Error: NumPy is required.", file=sys.stderr)
        sys.exit(1)
    if not args.no_superpose_pdb and PDBIO is None:
        print(
            "Error: BioPython PDBIO is required to write superposed PDBs. "
            "Use --no-superpose-pdb for scores/tree only, or install BioPython.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.fallback_biotite and not _ds.BIOTITE_AVAILABLE:
        print(
            "Error: --fallback-biotite requires biotite (pip install biotite).",
            file=sys.stderr,
        )
        sys.exit(1)

    if not _find_dalilite_bin(args.dalilite_path) and not args.fallback_biotite:
        print(
            "Error: DaliLite not found. Set DALILITE_HOME to the install root "
            "(containing bin/dali.pl) or pass --dalilite-path. "
            "Or use --fallback-biotite to run without DaliLite (biotite only).",
            file=sys.stderr,
        )
        sys.exit(1)

    out_dir = os.path.abspath(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    _dalilite_run_body(ap, args, out_dir, summary_log=summary_log)


def _dalilite_run_body(ap, args, out_dir: str, summary_log=None) -> None:
    def out_path(name: str | None) -> str | None:
        if not name:
            return None
        return name if os.path.isabs(name) else os.path.join(out_dir, name)

    def n_core_heatmap_kwargs() -> dict:
        """Argparse --heatmap-n-core-* into ``_write_dalilite_z_heatmap``."""
        if not args.heatmap_n_core_patterns:
            return {}
        if int(args.heatmap_n_core_bins) < 1:
            ap.error("--heatmap-n-core-bins must be at least 1")
        edges: list[float] | None = None
        if args.heatmap_n_core_edges:
            parts = [p.strip() for p in str(args.heatmap_n_core_edges).split(",") if p.strip()]
            if len(parts) < 2:
                ap.error("--heatmap-n-core-edges needs at least two comma-separated values")
            edges = [float(p) for p in parts]
        return {
            "n_core_heatmap": True,
            "n_core_n_equal_bins": max(1, int(args.heatmap_n_core_bins)),
            "n_core_bin_edges": edges,
        }

    verbose = not args.quiet
    write_superpose = not args.no_superpose_pdb

    if args.all_vs_all:
        if len(args.paths) < 2:
            ap.error("--all-vs-all requires at least two paths (dirs or structure files)")
        files = du._collect_pdb_files(args.paths, args.filter, recursive=args.recursive)
        if len(files) < 2:
            ap.error(f"need at least 2 structure files; found {len(files)}")
        if verbose:
            print(
                f"All-vs-all: {len(files)} structures, {len(files)*(len(files)-1)//2} pairs",
                file=sys.stderr,
            )
        if summary_log is not None:
            summary_log.task(f"All-vs-all: {len(files)} structures")
        data = run_all_pairs(
            files,
            out_dir,
            args.chain_a,
            args.chain_b,
            args.dalilite_path,
            args.keep_dalilite_work,
            write_superpose,
            verbose,
            debug_dalilite=args.debug_dalilite,
            fallback_biotite=args.fallback_biotite,
            equivalences_dir=out_path(args.equivalences_dir)
            if args.equivalences_dir
            else None,
        )
        zscores = data["zscores"]
        labels = data["labels"]
        if not zscores:
            print("Error: no pairwise DaliLite results.", file=sys.stderr)
            sys.exit(1)

        pairs_csv = out_path(args.output or "pairs.csv")
        assert pairs_csv
        with open(pairs_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "pdb_a",
                    "pdb_b",
                    "label_a",
                    "label_b",
                    "raw_score",
                    "z_score",
                    "n_core",
                    "core_resnum_min_a",
                    "core_resnum_max_a",
                    "core_resnum_min_b",
                    "core_resnum_max_b",
                    "alignment_source",
                    "dalilite_rmsd",
                    "core_rmsd",
                    "lali",
                    "nres",
                    "pct_id",
                    "dalilite_hit_id",
                    "dali_description",
                ]
            )
            for row in sorted(
                data["results_rows"],
                key=lambda r: (
                    _natural_sort_key(r["label_a"]),
                    _natural_sort_key(r["label_b"]),
                    _natural_sort_key(r["pdb_a"]),
                    _natural_sort_key(r["pdb_b"]),
                ),
            ):
                cr = row.get("core_rmsd")
                w.writerow(
                    [
                        row["pdb_a"],
                        row["pdb_b"],
                        row["label_a"],
                        row["label_b"],
                        f"{row['raw_score']:.6f}",
                        f"{row['z_score']:.6f}",
                        row["n_core"],
                        row.get("core_resnum_min_a", ""),
                        row.get("core_resnum_max_a", ""),
                        row.get("core_resnum_min_b", ""),
                        row.get("core_resnum_max_b", ""),
                        row.get("alignment_source", "dalilite"),
                        row["dalilite_rmsd"] if row["dalilite_rmsd"] is not None else "",
                        f"{cr:.4f}" if cr is not None else "",
                        row["lali"] if row["lali"] is not None else "",
                        row["nres"] if row.get("nres") is not None else "",
                        row["pct_id"] if row.get("pct_id") is not None else "",
                        row.get("dalilite_hit_id") or "",
                        row.get("dali_description") or "",
                    ]
                )
        print(f"Wrote {pairs_csv}")

        # Ranking CSV: rows strictly by Z (avg then max), not label or pair-table order.
        ranking = sorted(
            du._rank_structures(zscores),
            key=lambda r: (-float(r["avg_z"]), -float(r["max_z"]), str(r["label"])),
        )
        rank_p = out_path(args.output_ranking)
        assert rank_p
        with open(rank_p, "w", newline="") as f:
            w = csv.DictWriter(
                f, fieldnames=["label", "n_pairs", "avg_z", "max_z"]
            )
            w.writeheader()
            w.writerows(ranking)
        print(f"Wrote {rank_p}")

        mat_p = out_path(args.output_matrix)
        if mat_p:
            with open(mat_p, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["Model"] + labels)
                for la in labels:
                    row = [la]
                    for lb in labels:
                        z = zscores.get((la, lb), zscores.get((lb, la)))
                        row.append(f"{z:.4f}" if z is not None else "-")
                    w.writerow(row)
            print(f"Wrote {mat_p}")

        newick = du._build_tree_from_zscores(
            zscores,
            args.transform,
            args.exp_scale,
            args.root or "",
            not args.no_midpoint_root,
        )
        tree_p = out_path(args.output_tree)
        assert tree_p
        with open(tree_p, "w") as f:
            f.write(newick)
        print(f"Wrote {tree_p}")

        nomid_tree_p = out_path(args.output_tree_nomid) if args.output_tree_nomid else None
        nomid_plot_p = out_path(args.output_plot_nomid) if args.output_plot_nomid else None
        if nomid_tree_p and not args.no_midpoint_root:
            newick_nomid = du._build_tree_from_zscores(
                zscores,
                args.transform,
                args.exp_scale,
                args.root or "",
                False,
            )
            with open(nomid_tree_p, "w") as f:
                f.write(newick_nomid)
            print(f"Wrote {nomid_tree_p}")

            if nomid_plot_p:
                if du._plot_tree(
                    newick_nomid, nomid_plot_p, "DaliLite structural similarity tree (no midpoint root)"
                ):
                    print(f"Wrote {nomid_plot_p}")
                else:
                    print(
                        "Plot skipped (--output-plot-nomid): install ete3 or biopython+matplotlib",
                        file=sys.stderr,
                    )
        elif nomid_tree_p and args.no_midpoint_root:
            print(
                "Warning: --output-tree-nomid ignored (primary tree already uses --no-midpoint-root).",
                file=sys.stderr,
            )

        if args.heatmap:
            hp = out_path(args.heatmap)
            assert hp
            _write_dalilite_z_heatmap(
                labels,
                zscores,
                hp,
                title=args.heatmap_title,
                cmap=args.cmap,
                vmin=args.vmin,
                vmax=args.vmax,
                short_labels=args.short_heatmap_labels,
                diverging_center=args.heatmap_diverging_center,
                colorbar_orientation=args.heatmap_colorbar_orientation,
                y_axis_right=args.heatmap_y_axis_right,
                n_core=data["n_core"] if args.heatmap_n_core_patterns else None,
                **n_core_heatmap_kwargs(),
            )

        plot_p = out_path(args.output_plot)
        if plot_p:
            if du._plot_tree(
                newick, plot_p, "DaliLite structural similarity tree"
            ):
                print(f"Wrote {plot_p}")
            else:
                print(
                    "Plot skipped: install ete3 or biopython+matplotlib",
                    file=sys.stderr,
                )
        sys.exit(0)

    # Pairwise
    if len(args.paths) != 2:
        ap.error("pairwise mode: provide exactly two PDB/CIF paths")
    pdb_a, pdb_b = args.paths[0], args.paths[1]
    pair_tag = (
        f"{_sanitize_label(du._label_from_path(pdb_a))}__"
        f"{_sanitize_label(du._label_from_path(pdb_b))}"
    )
    work = os.path.join(out_dir, "_dalilite_work", pair_tag)
    cid_a, cid_b = _resolve_chain_ids(pdb_a, pdb_b, args.chain_a, args.chain_b)
    dali, derr = run_dalilite_pair(
        pdb_a,
        pdb_b,
        work,
        cid_a,
        cid_b,
        args.dalilite_path,
        debug_dalilite=args.debug_dalilite,
        fallback_biotite=args.fallback_biotite,
    )
    if dali is None:
        print(f"Error: DaliLite run failed — {derr}", file=sys.stderr)
        sys.exit(1)

    meta = score_pair_from_dalilite(pdb_a, pdb_b, args.chain_a, args.chain_b, dali)
    z_emp_pair = bool(meta.pop("z_is_empirical", False))
    equivs_pw = meta.pop("equivs", [])
    if args.equivalences_dir and equivs_pw:
        ed = out_path(args.equivalences_dir)
        os.makedirs(ed, exist_ok=True)
        du.write_equivalences_tsv(os.path.join(ed, f"{pair_tag}.tsv"), equivs_pw)
        print(f"Wrote equivalences: {ed}/{pair_tag}.tsv", file=sys.stderr)
    print(f"alignment: {meta.get('alignment_source', 'dalilite')}")
    print(f"raw_score: {meta['raw_score']:.4f}")
    if z_emp_pair:
        print(f"empirical Z-score: {meta['z_score']:.4f}  (FoldKit dali-like; not the DaliLite summary Z)")
    else:
        print(f"Z-score: {meta['z_score']:.4f}  (DaliLite reported)")
    print(f"n_core:    {meta['n_core']}")
    if meta.get("core_resnum_min_a") != "":
        print(
            f"core PDB residue span: A {meta.get('core_resnum_min_a')}-{meta.get('core_resnum_max_a')}, "
            f"B {meta.get('core_resnum_min_b')}-{meta.get('core_resnum_max_b')}"
        )
    if meta.get("dalilite_rmsd") is not None:
        print(f"dalilite_rmsd: {meta['dalilite_rmsd']:.2f} Å")
    if meta.get("core_rmsd") is not None:
        print(f"core_rmsd: {meta['core_rmsd']:.2f} Å (Dali or Kabsch on matched Cα)")

    if write_superpose and dali.get("rotation") is not None:
        super_dir = os.path.join(out_dir, "superposed")
        os.makedirs(super_dir, exist_ok=True)
        fix_p = os.path.join(super_dir, f"{pair_tag}_query.pdb")
        mob_p = os.path.join(super_dir, f"{pair_tag}_target_on_query.pdb")
        if write_fixed_copy(pdb_a, fix_p) and write_target_superposed_pdb(
            pdb_b, dali["rotation"], dali["translation"], mob_p
        ):
            print(f"Wrote {fix_p}")
            print(f"Wrote {mob_p}")
        else:
            print("Warning: could not write superposed PDBs.", file=sys.stderr)
    elif write_superpose:
        print(
            "Warning: no transrot block; superposed PDBs not written.",
            file=sys.stderr,
        )

    txt_copy = os.path.join(out_dir, f"{pair_tag}_dalilite.txt")
    with open(txt_copy, "w", encoding="utf-8") as f:
        f.write(dali["raw_text"])
    print(f"Wrote {txt_copy}")

    pair_csv = out_path(args.output or "pair.csv")
    if pair_csv:
        with open(pair_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "pdb_a",
                    "pdb_b",
                    "raw_score",
                    "z_score",
                    "n_core",
                    "core_resnum_min_a",
                    "core_resnum_max_a",
                    "core_resnum_min_b",
                    "core_resnum_max_b",
                    "alignment_source",
                    "dalilite_rmsd",
                    "core_rmsd",
                    "lali",
                    "nres",
                    "pct_id",
                    "dalilite_hit_id",
                    "dali_description",
                ]
            )
            cr = meta.get("core_rmsd")
            w.writerow(
                [
                    pdb_a,
                    pdb_b,
                    f"{meta['raw_score']:.6f}",
                    f"{meta['z_score']:.6f}",
                    meta["n_core"],
                    meta.get("core_resnum_min_a", ""),
                    meta.get("core_resnum_max_a", ""),
                    meta.get("core_resnum_min_b", ""),
                    meta.get("core_resnum_max_b", ""),
                    meta.get("alignment_source", "dalilite"),
                    meta["dalilite_rmsd"] if meta["dalilite_rmsd"] is not None else "",
                    f"{cr:.4f}" if cr is not None else "",
                    meta["lali"] if meta["lali"] is not None else "",
                    meta["nres"] if meta.get("nres") is not None else "",
                    meta["pct_id"] if meta.get("pct_id") is not None else "",
                    meta.get("dalilite_hit_id") or "",
                    meta.get("dali_description") or "",
                ]
            )
        print(f"Wrote {pair_csv}")

    if args.heatmap:
        la = du._label_from_path(pdb_a)
        lb = du._label_from_path(pdb_b)
        z = float(meta["z_score"])
        zm = {(la, lb): z, (lb, la): z}
        n_core_pair: dict | None = None
        if args.heatmap_n_core_patterns:
            nc = int(meta["n_core"])
            n_core_pair = {(la, lb): nc, (lb, la): nc}
        labels_pair = sorted([la, lb], key=_natural_sort_key)
        hp = out_path(args.heatmap)
        assert hp
        _write_dalilite_z_heatmap(
            labels_pair,
            zm,
            hp,
            title=args.heatmap_title,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            short_labels=args.short_heatmap_labels,
            diverging_center=args.heatmap_diverging_center,
            colorbar_orientation=args.heatmap_colorbar_orientation,
            y_axis_right=args.heatmap_y_axis_right,
            n_core=n_core_pair,
            **n_core_heatmap_kwargs(),
        )

    if not args.keep_dalilite_work:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
