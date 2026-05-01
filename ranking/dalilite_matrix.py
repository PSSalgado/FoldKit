#!/usr/bin/env python3
"""
DaliLite **matrix** mode driver (``dali.pl --matrix``).

This script **only** runs DaliLite (``import.pl`` + ``dali.pl --matrix``). Scores and trees
come from Dali’s own output files—there is no separate FoldKit/Dali re-scoring pass.

**Always written (Dali standard layout, copied or mirrored under the directory given with ``-d``):**
  native ``ordered`` matrix, per-query ``xxxxX.txt``, ``newick``, ``newick_unrooted``, etc.
  (see DaliLite manual for names).

  DaliLite’s Fortran tools enforce path lengths (often **≤ 80** characters) for
  ``import.pl``; the run therefore uses a **short directory under the system temp** folder
  (e.g. ``/tmp``), and copies results to ``-d``. With ``--keep-dalilite-work``, a full copy
  of that work tree is saved under ``<output>/dali_matrix_work/`` before the temporary directory is removed.

**Optional re-exports (all derived only from Dali’s files, for convenience):**
  pairwise details table, square ordered-matrix CSV, ranking from the matrix, and/or a
  matplotlib heatmap. Enable with ``--export-pair-table``, ``--export-z-matrix``,
  ``--export-ranking``, and/or ``--heatmap`` (see ``--help``).

Examples:
  # Dali outputs only
  python ranking/dalilite_matrix.py --all-vs-all /path/to/models/ /path/to/models/ \\
    -d /path/to/out --dalilite-path /path/to/DaliLite

  # Two heatmaps (different explicit paths): Z colours + n_core in cells; %%ID colours + %%ID in cells
  python ranking/dalilite_matrix.py --all-vs-all /path/to/models/ /path/to/models/ \\
    -d /path/to/out --dalilite-path /path/to/DaliLite \\
    --export-all-tables --heatmap dali_z.svg --heatmap-pct-id dali_pct_id.svg --heatmap-n-core-bins 4

  # Plot Dali’s native unrooted Newick (ete3 or biopython+matplotlib)
  python ranking/dalilite_matrix.py --all-vs-all /path/to/models/ /path/to/models/ \\
    -d /path/to/out --dalilite-path /path/to/DaliLite \\
    --output-dali-newick-plot dali_newick_unrooted.png
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

_REPO_ROOT = os.path.dirname(_SCRIPT_DIR)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# dali_score: shared helpers only; matrix values still come only from Dali --matrix (ordered + .txt).
import foldkit_dali_utils as _ds
from structure_phylogeny import _natural_sort_key

from utils.cli_log import add_log_args, setup_log_from_args


@dataclass(frozen=True)
class _Struct:
    path: str
    label: str
    dali_pdbid4: str
    dali_id5: str  # e.g. abcdA


def _ensure_unique_pdbid4(used: set[str], want: str) -> str:
    base = re.sub(r"[^A-Za-z0-9]", "", (want or "")[:4])
    if len(base) < 4:
        base = (base + "AAAA")[:4]
    cand = base[:4]
    n = 0
    while cand.lower() in used:
        n += 1
        # deterministic fallback: aaaa, aaab, ...
        suff = f"{n:04d}"[-3:]
        cand = (base[0] + base[1] + suff[:2])[:4]
    used.add(cand.lower())
    return cand


def _default_chain(pdb_path: str, explicit: str | None) -> str:
    if explicit:
        return explicit
    try:
        from Bio.PDB import PDBParser, MMCIFParser
    except Exception:
        return "A"
    p = os.path.abspath(pdb_path)
    try:
        if p.lower().endswith((".cif", ".mmcif")):
            s = MMCIFParser(QUIET=True).get_structure("m", p)
        else:
            s = PDBParser(QUIET=True).get_structure("m", p)
        chains = [c.id for c in s[0].get_chains()]
        if not chains:
            return "A"
        if "A" in chains:
            return "A"
        return chains[0]
    except Exception:
        return "A"


def _import_many_and_build_query_list(
    files: list[str],
    out_run_root: str,
    dat_dir: str,
    import_pl: str,
    *,
    default_chain: str | None,
    chain_from_args_a: str | None,
    verbose: bool = True,
) -> tuple[list[_Struct], str, str]:
    """
    Import all structures to ``dat_dir`` and write ``query.list`` into ``out_run_root``.

    Dali uses a 4-letter + chain identifier. Stable 4-letter ids are built from the filename stem,
    with collision handling, then append the chain letter as the 5th character, matching
    ``import.pl`` expectations (see DaliLite manual: ``--pdbid 1ppt`` + chain -> ``1pptA``).
    """
    used: set[str] = set()
    structs: list[_Struct] = []

    for fp in files:
        lab = _ds._label_from_path(fp)
        stem = os.path.splitext(os.path.basename(fp))[0]
        pid4 = _ensure_unique_pdbid4(used, stem)
        if chain_from_args_a is not None:
            ch = chain_from_args_a
        else:
            ch = _default_chain(fp, default_chain)
        d5 = f"{pid4}{ch}"
        structs.append(_Struct(path=os.path.abspath(fp), label=lab, dali_pdbid4=pid4, dali_id5=d5))

    # import each
    n_st = len(structs)
    for i, s in enumerate(structs, start=1):
        if verbose:
            print(
                f"[{i}/{n_st}] import.pl → {s.dali_id5}  ({s.label} ← {os.path.basename(s.path)})",
                file=sys.stderr,
            )
        staged = _ds._stage_pdb_for_dalilite_import(
            s.path, os.path.join(out_run_root, f"in_{s.dali_pdbid4}")
        )
        cmd = [
            import_pl,
            "--pdbfile",
            staged,
            "--pdbid",
            s.dali_pdbid4,
            "--dat",
            dat_dir,
            # Do NOT use --clean here: it would remove previously imported chains on each call.
        ]
        try:
            proc = subprocess_run(cmd, cwd=out_run_root, timeout=300)
        except Exception as e:
            raise SystemExit(f"import.pl failed to start: {e}") from e
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "")[-4000:]
            raise SystemExit(
                f"import.pl failed for {s.path} (pdbid {s.dali_pdbid4}): exit {proc.returncode}: {tail!r}"
            )
        if not any(
            f.startswith(s.dali_pdbid4) and f.endswith(".dat") for f in os.listdir(dat_dir)
        ):
            raise SystemExit(
                f"import.pl produced no {s.dali_pdbid4}*.dat in {dat_dir} for {s.path}. "
                f"stderr tail: {(proc.stderr or '')[-2000:]!r}"
            )

    # query list: first field is 5-char id; may include a free-text label after tabs/spaces
    qpath = os.path.join(out_run_root, "query.list")
    with open(qpath, "w", encoding="utf-8") as f:
        for s in structs:
            f.write(f"{s.dali_id5}\t{s.label}\n")
    return structs, qpath, dat_dir


def subprocess_run(cmd: list[str], cwd: str, timeout: int) -> "subprocess.CompletedProcess":
    import subprocess

    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def _detach_dalilite_workdir_temp(work: str, out_dir: str, keep_copy: bool) -> None:
    """
    DaliLite’s Fortran import step requires paths to stay short; ``work`` is under
    system temp. Optionally copy the full run tree to ``out_dir/dali_matrix_work``,
    then always remove ``work``.
    """
    if keep_copy and os.path.isdir(work):
        dest = os.path.join(out_dir, "dali_matrix_work")
        try:
            if os.path.isdir(dest):
                shutil.rmtree(dest, ignore_errors=True)
            shutil.copytree(work, dest, dirs_exist_ok=True)
        except OSError as e:
            print(f"Warning: could not copy Dali work dir to {dest!r}: {e}", file=sys.stderr)
    shutil.rmtree(work, ignore_errors=True)


def _copy_if_exists(src: str, dst_dir: str) -> str | None:
    if not src or not os.path.isfile(src):
        return None
    os.makedirs(dst_dir, exist_ok=True)
    base = os.path.basename(src)
    out = os.path.join(dst_dir, base)
    shutil.copy2(src, out)
    return out


def _read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def _parse_dali_ordered(path: str) -> tuple[list[str], "object"]:
    """
    Parse Dali ``ordered`` file (DaliLite v5 layout):
      n
      id1  v11 v12 ...
      ...
    Returns (row_ids_in_order, matrix ndarray float64) where M[i][j] matches row order.
    """
    try:
        import numpy as np
    except ImportError as e:
        raise SystemExit("NumPy is required to parse the Dali 'ordered' matrix.") from e

    lines = [ln for ln in _read_text(path).splitlines() if ln.strip()]
    if not lines:
        raise SystemExit(f"empty ordered file: {path}")
    n = int(lines[0].strip().split()[0])
    if n <= 0:
        raise SystemExit(f"invalid n in ordered file: {path}")
    if len(lines) < n + 1:
        raise SystemExit(f"ordered file too short: expected {n+1} lines, got {len(lines)}")

    row_ids: list[str] = []
    rows: list[list[float]] = []
    for ln in lines[1 : n + 1]:
        parts = ln.split()
        if not parts:
            continue
        rid = parts[0]
        row_ids.append(rid)
        rows.append([float(x) for x in parts[1:]])
    m = np.asarray(rows, dtype=float)
    if m.shape != (n, n):
        raise SystemExit(
            f"ordered matrix shape mismatch: got {m.shape} for n={n} from {path}"
        )
    return row_ids, m


def _align_matrix_to_labels(
    dali_ids: list[str], mat, labels: list[str], struct_by_dali5: dict[str, _Struct]
) -> tuple[list[str], "object", dict[str, str]]:
    """
    Reorder the parsed matrix to match ``labels`` (FoldKit), using 5- or 4-letter keys.
    """
    import numpy as np

    dmap: dict[str, int] = {}
    for i, did in enumerate(dali_ids):
        dmap[did] = i
        if len(did) >= 4:
            dmap.setdefault(did[:4], i)

    pos: list[int] = []
    used = set()
    dali5_for_label: dict[str, str] = {}
    for lab in labels:
        s = next((x for x in struct_by_dali5.values() if x.label == lab), None)
        if s is None:
            raise SystemExit(f"internal error: missing struct for label {lab!r}")
        for key in (s.dali_id5, s.dali_pdbid4 + s.dali_id5[-1:]):
            if key in dmap and dmap[key] not in used:
                pos.append(dmap[key])
                used.add(dmap[key])
                dali5_for_label[lab] = s.dali_id5
                break
        else:
            raise SystemExit(
                f"could not find matrix row for label {lab!r} (expected Dali id {s.dali_id5})"
            )
    m2 = mat[np.ix_(pos, pos)]
    return labels, m2, dali5_for_label


def _write_square_csv(path: str, labels: list[str], matrix) -> None:
    import numpy as np

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Model"] + labels)
        for i, la in enumerate(labels):
            row = [la]
            for j in range(len(labels)):
                v = float(matrix[i, j]) if (i, j) != (j, i) else float(matrix[i, j])
                if i == j:
                    row.append("-")
                else:
                    row.append(f"{v:.4f}" if np.isfinite(v) else "-")
            w.writerow(row)


def _build_pair_rows_from_matrix(
    files_by_label: dict[str, str],
    labels: list[str],
    mat,
    pair_extras: dict[tuple[str, str], dict] | None = None,
) -> list[dict]:
    import numpy as np

    out: list[dict] = []
    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = labels[i], labels[j]
            v = float(mat[i, j])
            if not np.isfinite(v):
                continue
            ex = (pair_extras or {}).get((a, b)) or (pair_extras or {}).get((b, a), {}) or {}
            out.append(
                {
                    "pdb_a": files_by_label.get(a, ""),
                    "pdb_b": files_by_label.get(b, ""),
                    "label_a": a,
                    "label_b": b,
                    "raw_score": v,
                    "z_score": v,
                    "n_core": ex.get("n_core", ""),
                    "core_resnum_min_a": ex.get("core_resnum_min_a", ""),
                    "core_resnum_max_a": ex.get("core_resnum_max_a", ""),
                    "core_resnum_min_b": ex.get("core_resnum_min_b", ""),
                    "core_resnum_max_b": ex.get("core_resnum_max_b", ""),
                    "alignment_source": ex.get("alignment_source", "dali_matrix"),
                    "dalilite_rmsd": ex.get("dalilite_rmsd"),
                    "core_rmsd": ex.get("core_rmsd"),
                    "lali": ex.get("lali"),
                    "nres": ex.get("nres"),
                    "pct_id": ex.get("pct_id"),
                    "dalilite_hit_id": ex.get("dalilite_hit_id", ""),
                    "dali_description": ex.get("dali_description", "from Dali --matrix 'ordered' similarity matrix"),
                }
            )
    return out


def _parse_float_list_csv(s: str | None) -> list[float] | None:
    if not s:
        return None
    parts = [p.strip() for p in str(s).split(",") if p.strip()]
    if not parts:
        return None
    return [float(p) for p in parts]


def _parse_dalilite_summary_hits_from_text(content: str) -> list[dict]:
    """
    Return DaliLite summary hit rows (``n:  hit_id  Z  rmsd  lali  ...``) from a Dali .txt.
    """
    out: list[dict] = []
    for ln in content.splitlines():
        s = ln.strip()
        if not re.match(r"^\d+:\s+\S", s):
            continue
        h = _ds.parse_dalilite_summary_hit(s + "\n")
        if h:
            out.append(h)
    return out


def _dali5_key_for_hit(hit_id: str) -> str:
    """
    Map a Dali summary ``hit_id`` token onto the same 5-char key space as struct.dali_id5
    (4-letter PDB + 1-letter chain), best-effort.
    """
    tok = (hit_id or "").strip()
    if not tok:
        return ""
    alnum = re.sub(r"[^A-Za-z0-9]", "", tok)
    if len(alnum) >= 5:
        return alnum[:4].lower() + alnum[4:5]
    if len(alnum) == 4:
        return alnum.lower() + "a"
    return (alnum + "aaaaa")[:5].lower()


def _best_hit_for_target(hits: list[dict], want_d5: str) -> dict | None:
    want = (want_d5 or "").strip()
    if not want or not hits:
        return None
    want_l = want.lower()
    want4 = want_l[:4]
    want_ch = want_l[-1:] if len(want_l) == 5 else ""
    cands: list[dict] = []
    for h in hits:
        hid = (h.get("hit_id") or "").strip()
        if not hid:
            continue
        k = _dali5_key_for_hit(hid)
        if k == want_l:
            cands.append(h)
            continue
        alnum = re.sub(r"[^A-Za-z0-9]", "", hid)
        if len(alnum) >= 5 and alnum[:4].lower() == want4 and alnum[4:5].lower() == (want_ch or alnum[4:5]).lower():
            cands.append(h)
    if not cands:
        return None
    cands.sort(key=lambda x: float(x.get("z_score", 0.0)), reverse=True)
    return cands[0]


def _read_pair_dalilite_hits(
    work_dir: str, structs: list[_Struct]
) -> tuple[dict[tuple[str, str], dict], list[str]]:
    """
    Parse each ``{dali_id5}.txt`` produced under work_dir, extract per-target hit metadata.

    Returns (extras_for_pairs, warnings).
    """
    extras: dict[tuple[str, str], dict] = {}
    warns: list[str] = []

    def _canon(a: str, b: str) -> tuple[str, str]:
        return (a, b) if a <= b else (b, a)

    for s in structs:
        p = os.path.join(work_dir, f"{s.dali_id5}.txt")
        if not os.path.isfile(p):
            warns.append(
                f"missing Dali per-query file {os.path.basename(p)} under {work_dir!r} "
                f"(skipping n_core / pair metadata parsed from Dali .txt files)"
            )
            continue
        txt = _read_text(p)
        hits = _parse_dalilite_summary_hits_from_text(txt)
        if not hits:
            continue
        for s2 in structs:
            if s2.dali_id5 == s.dali_id5:
                continue
            h = _best_hit_for_target(hits, s2.dali_id5)
            if h is None:
                continue
            key = _canon(s.label, s2.label)
            cz = float(h.get("z_score", 0.0) or 0.0)
            n_core_v = int(h.get("lali", 0) or 0)
            cand = {
                "n_core": n_core_v,
                "lali": h.get("lali"),
                "nres": h.get("nres"),
                "pct_id": h.get("pct_id"),
                "dalilite_rmsd": h.get("rmsd"),
                "core_rmsd": h.get("rmsd"),
                "dalilite_hit_id": h.get("hit_id", ""),
                "dali_description": (h.get("description") or "").strip()
                or "DaliLite summary line from matrix-mode per-query .txt",
                "alignment_source": "dalilite_matrix_txt",
            }
            prev = extras.get(key)
            if prev is None or float(prev.get("z_for_pick", 0.0) or 0.0) < cz:
                cand2 = dict(cand)
                cand2["z_for_pick"] = cz
                extras[key] = cand2
    for v in list(extras.values()):
        v.pop("z_for_pick", None)
    for (a, b), v in list(extras.items()):
        extras[(b, a)] = v
    return extras, warns


def _n_core_dict_from_extras(
    labels: list[str], extras: dict[tuple[str, str], dict]
) -> dict[tuple[str, str], int]:
    out: dict[tuple[str, str], int] = {}
    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = labels[i], labels[j]
            ex = extras.get((a, b), {})
            v = ex.get("n_core", None)
            if v is None or v == "":
                continue
            try:
                out[(a, b)] = int(v)
            except (TypeError, ValueError):
                continue
    return out


def _matrix_from_ncore(labels: list[str], n_core: dict[tuple[str, str], int]) -> list[list[int | float | None]]:
    m: list[list[int | float | None]] = []
    for i, la in enumerate(labels):
        row: list[int | float | None] = []
        for j, lb in enumerate(labels):
            if i == j:
                row.append(None)
                continue
            v = n_core.get((la, lb), n_core.get((lb, la)))
            row.append(int(v) if v is not None else None)
        m.append(row)
    return m


def _reorder_matrix_by_label_order(M, mat_labels: list[str], order: str):
    """
    Return (M_perm, labels_perm).
    - ``natural``: numeric-aware sort (S1, S2, …, S10; same key as structure_phylogeny / rmsd_to_csv).
    - ``alphanumeric``: plain ``sorted()`` (lexicographic: S1, S10, S2, …).
    - ``dali``: keep row/column order from Dali's ``ordered`` output after relabelling to model names.
    """
    o = (order or "natural").strip().lower()
    if o == "dali" or len(mat_labels) < 2:
        return M, mat_labels
    import numpy as np

    if o == "natural":
        new_labels = sorted(mat_labels, key=_natural_sort_key)
    elif o == "alphanumeric":
        new_labels = sorted(mat_labels)
    else:
        return M, mat_labels
    try:
        idx = [mat_labels.index(x) for x in new_labels]
    except ValueError as e:
        raise SystemExit(f"--label-order {order!r}: duplicate or missing label: {e}") from e
    m2 = M[np.ix_(idx, idx)]
    return m2, new_labels


def _pct_id_dict_from_extras(
    labels: list[str], extras: dict[tuple[str, str], dict]
) -> dict[tuple[str, str], float]:
    out: dict[tuple[str, str], float] = {}
    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = labels[i], labels[j]
            ex = extras.get((a, b), {})
            v = ex.get("pct_id", None)
            if v is None or v == "":
                continue
            try:
                out[(a, b)] = float(v)
            except (TypeError, ValueError):
                continue
    return out


def _matrix_from_pct_id(labels: list[str], pct: dict[tuple[str, str], float]) -> list[list[float | None]]:
    m: list[list[float | None]] = []
    for i, la in enumerate(labels):
        row: list[float | None] = []
        for j, lb in enumerate(labels):
            if i == j:
                row.append(None)
                continue
            v = pct.get((la, lb), pct.get((lb, la)))
            row.append(float(v) if v is not None else None)
        m.append(row)
    return m


def _n_core_hatch_bin_options_from_args(args) -> dict:
    """
    ``hatch_n_equal_bins`` and ``hatch_bin_edges`` for binned n_core hatches
    (Z heatmap and %ID heatmap second channel; values from Dali .txt per-query hits).
    """
    if int(args.heatmap_n_core_bins) < 1:
        raise SystemExit("--heatmap-n-core-bins must be at least 1")
    edges: list[float] | None = None
    if args.heatmap_n_core_edges:
        parts = [p.strip() for p in str(args.heatmap_n_core_edges).split(",") if p.strip()]
        if len(parts) < 2:
            raise SystemExit("--heatmap-n-core-edges needs at least two comma-separated values")
        edges = [float(p) for p in parts]
    return {
        "hatch_n_equal_bins": max(1, int(args.heatmap_n_core_bins)),
        "hatch_bin_edges": edges,
    }


def _n_core_heatmap_kwargs_from_args(args) -> dict:
    if not args.heatmap_n_core_patterns:
        return {}
    return _n_core_hatch_bin_options_from_args(args)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Run DaliLite ``dali.pl --matrix`` only. Copies Dali’s native outputs; "
            "optionally exports tables/heatmaps reformatting Dali’s ``ordered`` matrix and ``*.txt`` hits."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "paths",
        nargs="*",
        help="PDB/CIF files and/or directories (same collection rules as foldkit_dali_utils).",
    )
    ap.add_argument("-d", "--output-dir", required=True, metavar="DIR", help="Output directory.")
    ap.add_argument(
        "--all-vs-all",
        action="store_true",
        help="Use all structure files discovered from the given paths (unordered pairs implied by Dali --matrix).",
    )
    ap.add_argument(
        "--recursive",
        action="store_true",
        help="With --all-vs-all, search each directory argument recursively (default: root only).",
    )
    ap.add_argument(
        "--filter",
        metavar="PATTERN",
        help="Optional basename filter (see foldkit_dali_utils._collect_pdb_files).",
    )
    ap.add_argument("--chain", metavar="ID", help="Default chain for every structure (else auto-detect).")
    ap.add_argument("--chain-a", metavar="ID", help="(compat) if set, forces this chain for all models.")
    ap.add_argument("--chain-b", metavar="ID", help="Ignored in matrix mode (present for flag compatibility).")
    ap.add_argument("--dalilite-path", metavar="DIR", help="DaliLite install directory (or set DALILITE_HOME).")
    ap.add_argument(
        "--keep-dalilite-work",
        action="store_true",
        help="Copy the full Dali work tree from the short-path temp run into <output>/dali_matrix_work/ (then remove temp).",
    )
    ap.add_argument(
        "--dali-title",
        default="Dali matrix mode",
        help="Dali --title string.",
    )
    ap.add_argument(
        "--dali-np",
        type=int,
        default=1,
        help="Dali --np (MPI workers). Default 1 (serial).",
    )

    ap.add_argument(
        "-o",
        "--output",
        default="dali_matrix_pairs.csv",
        metavar="FILE",
        help="With --export-pair-table: pairwise CSV (Dali 'ordered' matrix + per-pair Dali .txt details when available).",
    )
    ap.add_argument(
        "--output-matrix-ordered",
        default="dali_matrix_ordered",
        help="Path prefix for the copied native 'ordered' file; '.txt' is appended if no extension is given.",
    )
    ap.add_argument(
        "--output-dali-newick",
        default="dali_newick.nwk",
        help="Copy/rename the native Dali 'newick' file to this name (if present).",
    )
    ap.add_argument(
        "--output-dali-newick-unrooted",
        default="dali_newick_unrooted.nwk",
        help="Copy/rename the native Dali 'newick_unrooted' file to this name (if present).",
    )
    ap.add_argument(
        "--output-dali-newick-plot",
        default=None,
        metavar="FILE",
        help="Optional PNG/SVG/PDF of the native Dali 'newick_unrooted' if ete3 or biopython+matplotlib is installed.",
    )
    ap.add_argument(
        "--output-dali-newick-plot-rooted",
        default=None,
        metavar="FILE",
        help="Optional PNG/SVG/PDF of the native Dali 'newick' (rooted) tree, if that file is present.",
    )
    ap.add_argument(
        "--output-ranking",
        default="dali_matrix_ranking.csv",
        metavar="FILE",
        help="With --export-ranking: ranking CSV (avg/max off-diagonals of Dali 'ordered' matrix).",
    )
    ap.add_argument(
        "--z-matrix",
        default="dali_matrix_ordered_square.csv",
        metavar="FILE",
        help="With --export-z-matrix: square CSV (Model + labels; values from Dali 'ordered' file).",
    )
    ap.add_argument(
        "--export-pair-table",
        action="store_true",
        help="After Dali: write a pairwise table CSV (see -o) from Dali's ordered matrix + per-query .txt details.",
    )
    ap.add_argument(
        "--export-z-matrix",
        action="store_true",
        help="After Dali: write a square matrix CSV of Dali's 'ordered' matrix (relabelled to supplied model names).",
    )
    ap.add_argument(
        "--export-ranking",
        action="store_true",
        help="After Dali: write a ranking CSV from avg/max of off-diagonal values in the ordered matrix.",
    )
    ap.add_argument(
        "--export-all-tables",
        action="store_true",
        help="Shorthand for --export-pair-table --export-z-matrix --export-ranking.",
    )
    ap.add_argument(
        "--label-order",
        choices=("natural", "alphanumeric", "dali"),
        default="natural",
        help=(
            "Row/column order for exported square matrix, pairwise table, ranking, and heatmaps. "
            "natural (default): S1, S2, S3, … before S10, S11 (numeric-aware; same as rmsd_to_csv). "
            "alphanumeric: strict string order — S10 and S11 come before S2 (often undesirable for numbered labels). "
            "dali: keep Dali's 'ordered' matrix row order after relabelling to model names."
        ),
    )

    # Heatmap
    ap.add_argument(
        "--heatmap",
        metavar="PATH",
        default=None,
        help="Optional: matplotlib image of Dali's 'ordered' matrix (PNG/SVG/PDF from file extension).",
    )
    ap.add_argument(
        "--heatmap-pct-id",
        default=None,
        metavar="PATH",
        help=(
            "Output path for the second heatmap: %%ID as colour + %%ID in cells (Dali per-query .txt). "
            "Must be a different path than --heatmap. Omit to write only the Z heatmap."
        ),
    )
    ap.add_argument(
        "--pct-id-heatmap",
        action="store_true",
        help=(
            "Ignored except for backward compatibility checks. Write the %%ID heatmap by passing "
            "--heatmap-pct-id PATH (required; must differ from --heatmap)."
        ),
    )
    ap.add_argument(
        "--heatmap-pct-id-title",
        default="Dali %ID (per-query .txt)",
        help="Title for --heatmap-pct-id (figure title; colour bar is fixed to match).",
    )
    ap.add_argument(
        "--heatmap-pct-id-fmt",
        default="{:.1f}",
        metavar="FMT",
        help="Format for %%ID values drawn in cells on the %%ID heatmap (colour still encodes %%ID).",
    )
    ap.add_argument(
        "--heatmap-pct-id-vmin",
        type=float,
        default=None,
        metavar="VAL",
        help="Lower colour scale limit for the %%ID heatmap (e.g. 50).",
    )
    ap.add_argument(
        "--heatmap-pct-id-vmax",
        type=float,
        default=None,
        metavar="VAL",
        help="Upper colour scale limit for the %%ID heatmap (e.g. 85).",
    )
    ap.add_argument(
        "--heatmap-pct-id-hatch-name",
        default=None,
        metavar="STR",
        help=(
            "Unused (reserved). The %%ID heatmap has no n_core hatch; use --heatmap-hatch-name / --heatmap-n-core-patterns "
            "only on the Z heatmap."
        ),
    )
    ap.add_argument(
        "--heatmap-title",
        default="DaliLite pairwise Z-scores",
        help="Figure title for --heatmap (same default as dalilite_pairs).",
    )
    ap.add_argument("--cmap", default="viridis_r", help="Heatmap colormap.")
    ap.add_argument("--vmin", type=float, default=None)
    ap.add_argument("--vmax", type=float, default=None)
    ap.add_argument("--short-heatmap-labels", action="store_true")
    ap.add_argument(
        "--heatmap-diverging-center",
        "--heatmap-diverging-centre",
        choices=("none", "median"),
        default="none",
    )
    ap.add_argument(
        "--heatmap-colorbar-orientation",
        "--heatmap-colourbar-orientation",
        choices=("vertical", "horizontal"),
        default="vertical",
    )
    ap.add_argument("--heatmap-y-axis-right", action="store_true")
    ap.add_argument(
        "--heatmap-cbar-label",
        default="Dali Z",
        help="Colour bar label for the Z / ordered-matrix heatmap (same as dalilite_pairs).",
    )
    ap.add_argument(
        "--heatmap-cbar-ticks",
        default=None,
        metavar="T0,T1,…",
        help="Optional comma-separated colour bar tick positions (foldkit_heatmap: --cbar-ticks).",
    )
    ap.add_argument(
        "--heatmap-cbar-tick-step",
        type=float,
        default=None,
        metavar="STEP",
        help="Optional colour bar tick step (foldkit_heatmap: --cbar-tick-step).",
    )
    ap.add_argument(
        "--heatmap-triangle",
        choices=("full", "lower"),
        default="full",
        help="Plot the full matrix, or the lower triangle including the diagonal (foldkit_heatmap: --triangle).",
    )
    ap.add_argument(
        "--heatmap-cell-border-linewidth-pt",
        type=float,
        default=0.3,
        metavar="PT",
        help="Per-cell border line width in points (foldkit_heatmap default: 0.3).",
    )
    ap.add_argument(
        "--heatmap-cell-border-color",
        "--heatmap-cell-border-colour",
        default="black",
        help="Per-cell border colour.",
    )
    ap.add_argument(
        "--heatmap-autoscale-positive-offdiag-only",
        action="store_true",
        help="If set, autoscale the colour range from positive off-diagonal values only. "
        "By default, all finite off-diagonal values are used (recommended for the Dali 'ordered' matrix).",
    )
    ap.add_argument(
        "--heatmap-annotate",
        choices=("none", "main", "hatch", "main+hatch", "auto"),
        default="auto",
        help=(
            "Cell labels on the Z heatmap: none, main (Z only), hatch (n_core numbers only; Z is colour only), "
            "main+hatch (Z + n_core on two lines), or auto (hatch when n_core is available from Dali .txt, else none). "
            "Use --heatmap-n-core-patterns only to add binned n_core hatches on the Z figure."
        ),
    )
    ap.add_argument(
        "--heatmap-annotate-fmt",
        default="{:.1f}",
        help="Format for the main (Z) channel when annotating (foldkit_heatmap: --annotate-fmt).",
    )
    ap.add_argument(
        "--heatmap-annotate-hatch-fmt",
        default="{:.0f}",
        metavar="FMT",
        help="Format for n_core in cells: second line of main+hatch, or the only line when annotate is hatch/auto.",
    )
    ap.add_argument(
        "--heatmap-annotate-fontsize",
        type=float,
        default=6.0,
        help="Font size (points) for --heatmap-annotate.",
    )
    ap.add_argument(
        "--heatmap-annotate-color-light",
        "--heatmap-annotate-colour-light",
        default="black",
        help="Annotation colour on light cells (foldkit_heatmap: --annotate-colour-light).",
    )
    ap.add_argument(
        "--heatmap-annotate-color-dark",
        "--heatmap-annotate-colour-dark",
        default="white",
        help="Annotation colour on dark cells (foldkit_heatmap: --annotate-colour-dark).",
    )
    ap.add_argument(
        "--heatmap-annotate-color-threshold",
        "--heatmap-annotate-colour-threshold",
        type=float,
        default=0.45,
        metavar="LUMA",
        help="Luminance threshold for auto text colour (foldkit_heatmap: --annotate-colour-threshold).",
    )
    ap.add_argument(
        "--heatmap-n-core-patterns",
        action="store_true",
        help=(
            "Optional second channel: bin n_core (Dali lali, aligned Cα) into matplotlib hatch patterns. "
            "Values are read from the per-query Dali .txt files when available (best-effort)."
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
        default=None,
        metavar="E0,E1,…",
        help="With --heatmap-n-core-patterns: optional comma-separated bin boundaries (low to high; extended to data min/max as needed).",
    )
    ap.add_argument(
        "--heatmap-hatch-name",
        default="N. aligned Cα",
        metavar="STR",
        help="Legend title for n_core hatch patterns when --heatmap-n-core-patterns is set (foldkit_heatmap: --hatch-name).",
    )
    ap.add_argument(
        "--heatmap-hatch-repeat-cells",
        type=int,
        default=4,
        metavar="N",
        help="Hatch line density in cells (foldkit_heatmap: --hatch-repeat-cells; smaller = more spaced lines).",
    )
    ap.add_argument(
        "--heatmap-hatch-repeat-legend",
        type=int,
        default=5,
        metavar="N",
        help="Hatch line density in legend patches (foldkit_heatmap: --hatch-repeat-legend).",
    )
    ap.add_argument(
        "--heatmap-hatch-linewidth-pt",
        type=float,
        default=0.15,
        metavar="PT",
        help="Hatch line width in points (foldkit_heatmap: --hatch-linewidth-pt).",
    )
    ap.add_argument(
        "--heatmap-hatch-color-light",
        "--heatmap-hatch-colour-light",
        default="black",
        help="Hatch colour on light cells (foldkit_heatmap: --hatch-colour-light).",
    )
    ap.add_argument(
        "--heatmap-hatch-color-dark",
        "--heatmap-hatch-colour-dark",
        default="white",
        help="Hatch colour on dark cells (foldkit_heatmap: --hatch-colour-dark).",
    )
    ap.add_argument(
        "--heatmap-hatch-color-threshold",
        "--heatmap-hatch-colour-threshold",
        type=float,
        default=0.45,
        metavar="LUMA",
        help="Luminance threshold for hatch line colour (foldkit_heatmap: --hatch-colour-threshold).",
    )
    ap.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Reduce console output (no per-model import or dali.pl progress on stderr).",
    )

    add_log_args(ap)
    args = ap.parse_args()

    if args.export_all_tables:
        args.export_pair_table = True
        args.export_z_matrix = True
        args.export_ranking = True

    if not _ds.BIOPYTHON_AVAILABLE or _ds.PDBParser is None:
        print("Error: BioPython is required.", file=sys.stderr)
        raise SystemExit(1)

    out_dir = os.path.abspath(
        os.path.expanduser(os.path.expandvars(str(args.output_dir or "").strip() or "."))
    )
    os.makedirs(out_dir, exist_ok=True)

    if bool(args.pct_id_heatmap) and not args.heatmap_pct_id:
        print(
            "Error: --pct-id-heatmap no longer sets an output name. Pass an explicit --heatmap-pct-id PATH "
            "(different from --heatmap) for the %%ID heatmap.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    if args.heatmap and args.heatmap_pct_id:
        zp = os.path.abspath(os.path.join(out_dir, str(args.heatmap)))
        ip = os.path.abspath(os.path.join(out_dir, str(args.heatmap_pct_id)))
        if zp == ip:
            print(
                "Error: --heatmap and --heatmap-pct-id must be two different output files (got the same path).",
                file=sys.stderr,
            )
            raise SystemExit(1)
    summary_log = setup_log_from_args(
        args, script_path=__file__, inputs=[x for x in [getattr(args, "input", None)] if x], pattern=None
    )
    if summary_log is not None:
        summary_log.task("DaliLite matrix mode (dali.pl --matrix)")

    if not args.all_vs_all or len(args.paths) < 1:
        ap.error("Provide at least one path with --all-vs-all (directories and/or files).")

    if str(getattr(args, "label_order", "") or "").strip().lower() == "alphanumeric" and not bool(
        getattr(args, "quiet", False)
    ):
        print(
            "Note: --label-order alphanumeric sorts names as plain strings (e.g. S10 before S2). "
            "For S1, S2, …, S10 use the default --label-order natural (or omit the flag).",
            file=sys.stderr,
        )

    files = _ds._collect_pdb_files(
        [os.path.abspath(p) for p in args.paths],
        args.filter,
        recursive=bool(args.recursive),
    )
    files = sorted(
        files,
        key=lambda p: (_natural_sort_key(_ds._label_from_path(p)), _natural_sort_key(p)),
    )
    if len(files) < 2:
        ap.error("need at least 2 structure files for a matrix run")

    verbose = not bool(args.quiet)
    if verbose:
        npr = len(files) * (len(files) - 1) // 2
        print(
            f"Dali matrix mode: {len(files)} structures ({npr} off-diagonals in ordered matrix; one Dali run).",
            file=sys.stderr,
        )

    bin_dir = _ds._find_dalilite_path(args.dalilite_path)
    if not bin_dir:
        ap.error("DaliLite not found: pass --dalilite-path or set DALILITE_HOME")
    import_pl = os.path.join(bin_dir, "import.pl")
    dali_pl = os.path.join(bin_dir, "dali.pl")
    if not os.path.isfile(import_pl) or not os.path.isfile(dali_pl):
        ap.error("DaliLite bin missing import.pl or dali.pl; check --dalilite-path")

    pre = _ds._dalilite_preflight_mkdssp(os.path.join(bin_dir))
    if pre:
        print("Error: " + pre, file=sys.stderr)
        raise SystemExit(1)

    work = tempfile.mkdtemp(prefix="dlm_", dir=tempfile.gettempdir())
    dat_dir = os.path.join(work, "DAT")
    os.makedirs(dat_dir, exist_ok=True)
    try:
        chain_default = args.chain
        chain_a = args.chain_a
        # chain-b is ignored, but if user set only chain-b in error, use it
        if chain_a is None and args.chain_b and not args.chain:
            chain_a = args.chain_b
    
        structs, qlist, dat_dir2 = _import_many_and_build_query_list(
            files,
            work,
            dat_dir,
            import_pl,
            default_chain=chain_default,
            chain_from_args_a=chain_a,
            verbose=verbose,
        )
        assert dat_dir2 == dat_dir
    
        # Run dali.pl --matrix
        lock = os.path.join(work, "dali.lock")
        if os.path.exists(lock):
            try:
                os.remove(lock)
            except OSError:
                pass
        np_par = int(args.dali_np) if int(args.dali_np) > 0 else 1
        dcmd = [dali_pl, "--matrix", "--query", qlist, "--dat1", dat_dir, "--title", str(args.dali_title), "--clean"]
        if np_par != 1:
            dcmd = [dali_pl, "--np", str(np_par), *dcmd[1:]]
        if verbose:
            print(
                f"Running dali.pl --matrix (np={np_par}, {len(structs)} models, work={work!r}) …",
                file=sys.stderr,
            )
        try:
            dproc = subprocess_run(dcmd, cwd=work, timeout=3600 * 6)
        except Exception as e:
            print(f"Error: dali.pl failed to start: {e}", file=sys.stderr)
            raise SystemExit(1) from e
        if dproc.returncode != 0:
            tail = (dproc.stderr or dproc.stdout or "")[-4000:]
            print(f"Error: dali.pl exit {dproc.returncode}: {tail!r}", file=sys.stderr)
            raise SystemExit(1)
        if verbose:
            print("dali.pl --matrix finished.", file=sys.stderr)
    
        # Native artifacts: ordered/newick* + per-cd .txt
        def _p(name: str) -> str:
            return os.path.join(work, name)
    
        ordered_path = _p("ordered")
        newick_path = _p("newick")
        newick_u_path = _p("newick_unrooted")
        for must in (ordered_path,):
            if not os.path.isfile(must):
                # Some installs may use different file names: try case-insensitive search in cwd
                tail = (dproc.stderr or dproc.stdout or "")[-2000:]
                raise SystemExit(
                    f"expected {must} after dali --matrix, but it is missing. stdout/stderr tail: {tail!r}"
                )
    
        # ------------------------------------------------------------------
        # Dali standard outputs: copy Dali’s native files (no re-scoring)
        # ------------------------------------------------------------------
        o_prefix = args.output_matrix_ordered
        if o_prefix and not os.path.splitext(o_prefix)[1]:
            o_prefix = o_prefix + ".txt"
        native_dir = os.path.join(out_dir, "dali_native")
        c1 = _copy_if_exists(ordered_path, native_dir)
        if c1 and o_prefix and not os.path.isabs(args.output_matrix_ordered):
            tdest = os.path.join(out_dir, os.path.basename(o_prefix)) if o_prefix else None
            if tdest and (not os.path.exists(tdest)):
                shutil.copy2(ordered_path, tdest)
                print(f"Wrote {tdest}")
        if c1:
            print(f"Wrote {c1}")
        c2 = _copy_if_exists(newick_path, native_dir)
        c3 = _copy_if_exists(newick_u_path, native_dir)
        if c2:
            dst2 = os.path.join(out_dir, args.output_dali_newick)
            shutil.copy2(newick_path, dst2)
            print(f"Wrote {dst2}")
        if c3:
            dst3 = os.path.join(out_dir, args.output_dali_newick_unrooted)
            shutil.copy2(newick_u_path, dst3)
            print(f"Wrote {dst3}")
    
        if args.output_dali_newick_plot and c3:
            pout = os.path.join(out_dir, args.output_dali_newick_plot)
            if not _ds._plot_tree(
                _read_text(newick_u_path), pout, "Dali matrix mode (native newick_unrooted)"
            ):
                print("Warning: could not render native newick (install ete3 or biopython+matplotlib).", file=sys.stderr)
            else:
                print(f"Wrote {pout}")
    
        if args.output_dali_newick_plot_rooted and c2:
            pr = os.path.join(out_dir, args.output_dali_newick_plot_rooted)
            if not _ds._plot_tree(
                _read_text(newick_path), pr, "Dali matrix mode (native newick, rooted)"
            ):
                print("Warning: could not render native rooted newick (install ete3 or biopython+matplotlib).", file=sys.stderr)
            else:
                print(f"Wrote {pr}")
    
        # Optional: reformat Dali "ordered" / .txt (requires NumPy; heatmap also needs matplotlib)
        need_matrix = (
            bool(args.heatmap)
            or bool(args.heatmap_pct_id)
            or bool(args.export_z_matrix)
            or bool(args.export_pair_table)
            or bool(args.export_ranking)
        )
        need_pair_txt = bool(args.export_pair_table) or bool(args.heatmap) or bool(args.heatmap_pct_id)
    
        mat_labels: list[str] | None = None
        M = None
        files_by_label: dict[str, str] = {}
        pair_extras: dict = {}
    
        if need_matrix:
            if _ds.np is None:
                print(
                    "Error: NumPy is required for --export-*/--heatmap (reading Dali 'ordered' matrix).",
                    file=sys.stderr,
                )
                raise SystemExit(1)
            by5 = {s.dali_id5: s for s in structs}
            labels = [s.label for s in structs]
            dali_ids, mat0 = _parse_dali_ordered(ordered_path)
            mat_labels, M, _dali5_map = _align_matrix_to_labels(dali_ids, mat0, labels, by5)
            files_by_label = {s.label: s.path for s in structs}
            if need_pair_txt:
                pair_extras, pair_extra_warns = _read_pair_dalilite_hits(work, structs)
                for w in pair_extra_warns[:5]:
                    print(f"Warning: {w}", file=sys.stderr)
                if len(pair_extra_warns) > 5:
                    print(
                        f"Warning: {len(pair_extra_warns) - 5} additional Dali .txt warnings suppressed.",
                        file=sys.stderr,
                    )
            M, mat_labels = _reorder_matrix_by_label_order(M, mat_labels, str(args.label_order))
    
        if need_matrix and mat_labels is not None and M is not None:
            if args.export_z_matrix:
                zmat = os.path.join(out_dir, args.z_matrix)
                _write_square_csv(zmat, mat_labels, M)
                print(f"Wrote {zmat}")
    
            if args.export_pair_table:
                results_rows = _build_pair_rows_from_matrix(files_by_label, mat_labels, M, pair_extras)
                pairs_path = os.path.join(out_dir, args.output)
                with open(pairs_path, "w", newline="") as f:
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
                    for row in results_rows:
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
                                row.get("alignment_source", "dali_matrix"),
                                row["dalilite_rmsd"] if row["dalilite_rmsd"] is not None else "",
                                f"{cr:.4f}" if cr is not None else "",
                                row["lali"] if row["lali"] is not None else "",
                                row["nres"] if row.get("nres") is not None else "",
                                row["pct_id"] if row.get("pct_id") is not None else "",
                                row.get("dalilite_hit_id") or "",
                                row.get("dali_description") or "",
                            ]
                        )
                print(f"Wrote {pairs_path}")
    
            if args.export_ranking:
                zscores: dict = {}
                n = len(mat_labels)
                for i in range(n):
                    for j in range(i + 1, n):
                        a, b = mat_labels[i], mat_labels[j]
                        v = float(M[i, j])
                        zscores[(a, b)] = v
                ranking = sorted(
                    _ds._rank_structures(zscores),
                    key=lambda r: (
                        -float(r["avg_z"]),
                        -float(r["max_z"]),
                        _natural_sort_key(str(r["label"])),
                    ),
                )
                rank_p = os.path.join(out_dir, args.output_ranking)
                with open(rank_p, "w", newline="") as f:
                    w2 = csv.DictWriter(f, fieldnames=["label", "n_pairs", "avg_z", "max_z"])
                    w2.writeheader()
                    w2.writerows(ranking)
                print(f"Wrote {rank_p}")
    
        # Heatmap(s): Dali ordered matrix and optional %ID (both use --label-order for axes)
        if args.heatmap or args.heatmap_pct_id:
            if M is None or mat_labels is None:
                print(
                    "Error: heatmaps require loading Dali's ordered matrix (NumPy).",
                    file=sys.stderr,
                )
                raise SystemExit(1)
            from utils.foldkit_heatmap import plot_heatmap

        if args.heatmap:
            hp = os.path.join(out_dir, args.heatmap)
            mlist = M.tolist()
            cbar_ticks = _parse_float_list_csv(getattr(args, "heatmap_cbar_ticks", None))

            n_core_d = _n_core_dict_from_extras(mat_labels, pair_extras)
            hatch_mat: list[list[int | float | None]] | None = (
                _matrix_from_ncore(mat_labels, n_core_d) if n_core_d else None
            )
            if args.heatmap_n_core_patterns and not hatch_mat:
                print(
                    "Warning: --heatmap-n-core-patterns was set, but no per-pair n_core values "
                    "could be parsed from Dali's per-query .txt files; drawing Z heatmap without n_core hatches.",
                    file=sys.stderr,
                )

            z_mode = str(args.heatmap_annotate or "auto").strip().lower()
            if z_mode == "auto":
                # Z as colour; one numeric line = n_core only (not Z + n_core, and never %%ID on this figure).
                z_mode = "hatch" if hatch_mat is not None else "none"
            elif z_mode == "main+hatch" and hatch_mat is None:
                print(
                    "Warning: --heatmap-annotate main+hatch requires n_core; using main (Z only).",
                    file=sys.stderr,
                )
                z_mode = "main"
            elif z_mode == "hatch" and hatch_mat is None:
                print(
                    "Warning: --heatmap-annotate hatch requires n_core from Dali .txt; using none.",
                    file=sys.stderr,
                )
                z_mode = "none"
            z_annotate = None if z_mode in ("", "none") else z_mode

            z_cell_fmt = str(args.heatmap_annotate_fmt)
            if z_annotate == "hatch":
                z_cell_fmt = str(getattr(args, "heatmap_annotate_hatch_fmt", "{:.0f}"))

            hkwargs: dict = {}
            if hatch_mat is not None and (
                bool(args.heatmap_n_core_patterns)
                or z_annotate in ("main+hatch", "hatch")
            ):
                hkwargs["hatch_value_matrix"] = hatch_mat
                hkwargs["hatch_draw_patterns"] = bool(args.heatmap_n_core_patterns)
                hkwargs["hatch_channel_name"] = str(args.heatmap_hatch_name)
                if args.heatmap_n_core_patterns:
                    hkwargs.update(_n_core_heatmap_kwargs_from_args(args))
                    hkwargs["hatch_repeat_cells"] = max(1, int(args.heatmap_hatch_repeat_cells))
                    hkwargs["hatch_repeat_legend"] = max(1, int(args.heatmap_hatch_repeat_legend))
                    hkwargs["hatch_linewidth_pt"] = float(args.heatmap_hatch_linewidth_pt)
                    hkwargs["hatch_color_light"] = str(args.heatmap_hatch_color_light)
                    hkwargs["hatch_color_dark"] = str(args.heatmap_hatch_color_dark)
                    hkwargs["hatch_color_luma_threshold"] = float(args.heatmap_hatch_color_threshold)

            plot_heatmap(
                mat_labels,
                mlist,
                hp,
                title=args.heatmap_title,
                cmap=str(args.cmap),
                vmin=args.vmin,
                vmax=args.vmax,
                short_labels=bool(args.short_heatmap_labels),
                diverging_center=str(args.heatmap_diverging_center or "none"),
                vcenter=None,
                colorbar_orientation=str(args.heatmap_colorbar_orientation),
                y_axis_right=bool(args.heatmap_y_axis_right),
                cbar_label=str(args.heatmap_cbar_label),
                cbar_ticks=cbar_ticks,
                cbar_tick_step=getattr(args, "heatmap_cbar_tick_step", None),
                autoscale_positive_offdiag_only=bool(args.heatmap_autoscale_positive_offdiag_only),
                triangle=str(args.heatmap_triangle),
                cell_border_linewidth_pt=float(args.heatmap_cell_border_linewidth_pt),
                cell_border_color=str(args.heatmap_cell_border_color),
                annotate=z_annotate,
                annotate_fmt=z_cell_fmt,
                annotate_hatch_fmt=str(getattr(args, "heatmap_annotate_hatch_fmt", "{:.0f}")),
                annotate_fontsize=float(args.heatmap_annotate_fontsize),
                annotate_color_light=str(args.heatmap_annotate_color_light),
                annotate_color_dark=str(args.heatmap_annotate_color_dark),
                annotate_color_luma_threshold=float(args.heatmap_annotate_color_threshold),
                **hkwargs,
            )
            print(f"Wrote {hp}", flush=True)

        if args.heatmap_pct_id:
            pct_d = _pct_id_dict_from_extras(mat_labels, pair_extras)
            if not pct_d:
                print(
                    "Warning: no %ID in Dali .txt hit lines; skipped --heatmap-pct-id",
                    file=sys.stderr,
                )
            else:
                pm = _matrix_from_pct_id(mat_labels, pct_d)
                pp = os.path.join(out_dir, args.heatmap_pct_id)
                cbar_t = _parse_float_list_csv(getattr(args, "heatmap_cbar_ticks", None))
                plot_heatmap(
                    mat_labels,
                    pm,
                    pp,
                    title=str(args.heatmap_pct_id_title),
                    cmap=str(args.cmap),
                    vmin=getattr(args, "heatmap_pct_id_vmin", None),
                    vmax=getattr(args, "heatmap_pct_id_vmax", None),
                    short_labels=bool(args.short_heatmap_labels),
                    diverging_center="none",
                    vcenter=None,
                    colorbar_orientation=str(args.heatmap_colorbar_orientation),
                    y_axis_right=bool(args.heatmap_y_axis_right),
                    cbar_label="Dali %ID (per-query .txt)",
                    cbar_ticks=cbar_t,
                    cbar_tick_step=getattr(args, "heatmap_cbar_tick_step", None),
                    autoscale_positive_offdiag_only=False,
                    triangle=str(args.heatmap_triangle),
                    cell_border_linewidth_pt=float(args.heatmap_cell_border_linewidth_pt),
                    cell_border_color=str(args.heatmap_cell_border_color),
                    annotate="main",
                    annotate_fmt=str(getattr(args, "heatmap_pct_id_fmt", "{:.1f}")),
                    annotate_fontsize=float(args.heatmap_annotate_fontsize),
                    annotate_color_light=str(args.heatmap_annotate_color_light),
                    annotate_color_dark=str(args.heatmap_annotate_color_dark),
                    annotate_color_luma_threshold=float(args.heatmap_annotate_color_threshold),
                )
                print(f"Wrote {pp}", flush=True)

    finally:
        _detach_dalilite_workdir_temp(work, out_dir, args.keep_dalilite_work)


if __name__ == "__main__":
    main()
