#!/usr/bin/env python3
"""
Utilities for running and parsing **true DaliLite** outputs.

This module is intentionally **not** the "FoldKit Dali-like empirical scoring" implementation.
It contains only helpers needed by:
  - `ranking/dalilite_matrix.py`
  - `ranking/dalilite_pairs.py`

Most functions are extracted from the previous `ranking/foldkit_dali_like_scores.py` to avoid coupling
true-DaliLite pipelines to the empirical Z-like scoring script.
"""

from __future__ import annotations

import csv
import fnmatch
import glob
import io
import os
import re
import shutil
import subprocess
import tempfile
import warnings

from structure_phylogeny import _natural_sort_key

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

try:
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.PDBIO import PDBIO, Select
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning, PDBIOException

    warnings.simplefilter("ignore", PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:  # pragma: no cover
    MMCIFParser = None  # type: ignore[misc, assignment]
    PDBIO = None  # type: ignore[misc, assignment]
    Select = None  # type: ignore[misc, assignment]
    PDBParser = None
    PDBIOException = None  # type: ignore[misc, assignment]
    BIOPYTHON_AVAILABLE = False


def _label_row_sort_key(label: str) -> tuple:
    """Stable label order matching rmsd_to_csv / structure_phylogeny."""
    return (_natural_sort_key(label), label)


def _canonical_label_pair(la: str, lb: str) -> tuple[str, str]:
    """Undirected edge key: natural order, then str for ties."""
    ka, kb = _label_row_sort_key(la), _label_row_sort_key(lb)
    return (la, lb) if ka <= kb else (lb, la)


def write_equivalences_tsv(path: str, equivs) -> None:
    """
    One row per structurally equivalent pair (query / target) after normalisation.
    Columns: chain_A, resnum_A, icode_A, chain_B, resnum_B, icode_B (BioPython keys).
    """
    out_dir = os.path.dirname(path) or "."
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", newline="") as f:
        f.write(
            "# Structural core: aligned Cα pairs used for raw Dali score and n_core.\n"
            "# Z_score / lali / nres / pct_id come from the DaliLite summary line when present.\n"
        )
        f.write("chain_A\tresnum_A\ticode_A\tchain_B\tresnum_B\ticode_B\n")
        for ka, kb in equivs:
            ca, ra, ia = ka
            cb, rb, ib = kb
            ia = ia if isinstance(ia, str) else (ia or " ")
            ib = ib if isinstance(ib, str) else (ib or " ")
            f.write(f"{ca!s}\t{int(ra)}\t{ia!s}\t{cb!s}\t{int(rb)}\t{ib!s}\n")


def normalize_equivalences(equivs, coords_A, coords_B):
    """
    Convert equivalences to keys that exist in coords_A and coords_B.
    For single-chain: (chain, resnum, icode) may use ' ' for chain; match by resnum.
    """
    by_resnum_A = {}
    for k in coords_A.keys():
        chain, resnum, _icode = k
        by_resnum_A[(chain, resnum)] = k
        by_resnum_A[resnum] = k  # fallback for alignment with resnum only
    by_resnum_B = {}
    for k in coords_B.keys():
        chain, resnum, _icode = k
        by_resnum_B[(chain, resnum)] = k
        by_resnum_B[resnum] = k

    out = []
    for eA, eB in equivs:
        cA, rA, _iA = eA
        cB, rB, _iB = eB
        keyA = (cA, rA) if cA and cA != " " else rA
        keyB = (cB, rB) if cB and cB != " " else rB
        kA = by_resnum_A.get(keyA)
        kB = by_resnum_B.get(keyB)
        if kA is not None and kB is not None:
            out.append((kA, kB))
    return out


def _find_dalilite_path(dalilite_path: str | None = None) -> str:
    """Return path to DaliLite bin directory, or empty string if not found."""
    candidates = [p for p in (dalilite_path, os.environ.get("DALILITE_HOME", "")) if p]
    for base in candidates:
        for check in (base, os.path.join(base, "bin")):
            dali_pl = os.path.join(check, "dali.pl")
            if os.path.isfile(dali_pl):
                return os.path.dirname(dali_pl)
    return ""


def _parse_dalilite_equivalences(txt: str, chain_a: str = "A", chain_b: str = "A"):
    """
    Parse DaliLite structural equivalences block.
    Returns list of ((chain_a, resnum_a, ' '), (chain_b, resnum_b, ' ')).
    """
    equivs = []
    pattern = re.compile(
        r"^\s*\d+:\s+\S+\s+\S+\s+(\d+)\s*-\s*(\d+)(?:\s*<=\>\s*|\s+)(\d+)\s*-\s*(\d+)",
        re.MULTILINE,
    )
    for m in pattern.finditer(txt):
        q_start, q_end = int(m.group(1)), int(m.group(2))
        t_start, t_end = int(m.group(3)), int(m.group(4))
        n = min(q_end - q_start + 1, t_end - t_start + 1)
        for i in range(n):
            equivs.append(((chain_a, q_start + i, " "), (chain_b, t_start + i, " ")))
    return equivs


def parse_dalilite_summary_hit(content: str) -> dict | None:
    """
    Parse the first DaliLite summary hit line:
      ``  1:  xxxxX  Z  rmsd  lali  nres  %id  Description...``
    """
    m = re.search(
        r"^\s*\d+:\s+(\S+)\s+"
        r"([-\d.Ee+]+)\s+"
        r"([\d.Ee+]+)\s+"
        r"(\d+)\s+"
        r"(\d+)\s+"
        r"(\d+)\s*"
        r"(.*)$",
        content,
        re.MULTILINE,
    )
    if m:
        desc = (m.group(7) or "").strip()
        return {
            "hit_id": m.group(1),
            "z_score": float(m.group(2)),
            "rmsd": float(m.group(3)),
            "lali": int(m.group(4)),
            "nres": int(m.group(5)),
            "pct_id": int(m.group(6)),
            "description": desc,
        }
    m2 = re.search(
        r"^\s*\d+:\s+(\S+)\s+([-\d.Ee+]+)\s+([\d.Ee+]+)\s+(\d+)",
        content,
        re.MULTILINE,
    )
    if not m2:
        return None
    return {
        "hit_id": m2.group(1),
        "z_score": float(m2.group(2)),
        "rmsd": float(m2.group(3)),
        "lali": int(m2.group(4)),
        "nres": None,
        "pct_id": None,
        "description": "",
    }


def _dalilite_text_has_summary_hit(txt: str) -> bool:
    """True if DaliLite summary contains at least one hit line."""
    return parse_dalilite_summary_hit(txt) is not None


def _get_first_chain_id(coords_dict) -> str:
    """Return the first chain ID seen in coords, or 'A'."""
    if not coords_dict:
        return "A"
    for k in coords_dict:
        if k[0] and k[0] != " ":
            return str(k[0])
    return "A"


def _read_dalilite_dssp_exe_from_mpidali(bin_dir: str) -> str:
    pm = os.path.join(bin_dir, "mpidali.pm")
    try:
        with open(pm, "r", encoding="utf-8", errors="replace") as f:
            head = f.read(20000)
    except OSError:
        return ""
    m = re.search(r'^\s*my\s+\$DSSP_EXE\s*=\s*"([^"]+)"\s*;', head, re.MULTILINE)
    return m.group(1) if m else ""


def _dalilite_preflight_mkdssp(bin_dir: str) -> str | None:
    exe = _read_dalilite_dssp_exe_from_mpidali(bin_dir)
    if not exe:
        return None
    if os.path.isfile(exe) and os.access(exe, os.X_OK):
        return None
    return (
        "DaliLite import.pl requires mkdssp (see bin/mpidali.pm $DSSP_EXE). "
        f"Expected {exe!r} but it is missing or not executable."
    )


def _find_mkdssp_executable(dalilite_bin_dir: str | None = None) -> str:
    env_mk = (os.environ.get("MKDSSP") or "").strip()
    if env_mk and os.path.isfile(env_mk) and os.access(env_mk, os.X_OK):
        return env_mk
    if dalilite_bin_dir:
        exe = _read_dalilite_dssp_exe_from_mpidali(dalilite_bin_dir)
        if exe and os.path.isfile(exe) and os.access(exe, os.X_OK):
            return exe
    mk = shutil.which("mkdssp")
    return mk or ""


def _mkdssp_diagnostic_suffix(pdb_path: str, workdir: str, dalilite_bin_dir: str | None = None) -> str:
    mk = _find_mkdssp_executable(dalilite_bin_dir)
    if not mk:
        return " | mkdssp probe: no usable mkdssp (set MKDSSP or fix mpidali.pm $DSSP_EXE)."
    try:
        fd, out = tempfile.mkstemp(suffix=".dssp", dir=workdir)
        os.close(fd)
    except (OSError, FileNotFoundError) as e:
        return f" | mkdssp probe: could not create temp dssp in {workdir!r}: {e}"
    try:
        proc = subprocess.run([mk, pdb_path, out], capture_output=True, text=True, timeout=120)
        sz = os.path.getsize(out) if os.path.isfile(out) else 0
        tail_err = (proc.stderr or "")[-800:]
        if proc.returncode == 0 and sz > 0:
            return f" | mkdssp probe OK ({mk} wrote {sz} bytes) but import.pl still produced no .dat"
        return f" | mkdssp probe ({mk}): exit={proc.returncode} dssp_bytes={sz} stderr_tail={tail_err!r}"
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        return f" | mkdssp probe failed to run {mk!r}: {e}"
    finally:
        try:
            os.remove(out)
        except OSError:
            pass


def _stage_pdb_for_dalilite_import(src_path: str, dest_base: str) -> str:
    """
    Copy or rewrite structure into dest_base + extension for DaliLite import.
    (Same behaviour as the original helper in `foldkit_dali_like_scores.py`.)
    """
    ext = os.path.splitext(src_path)[1].lower()
    if ext not in (".cif", ".mmcif", ".mcif", ".pdb", ".ent", ""):
        ext = ".pdb"
    if not ext:
        ext = ".pdb"
    dest = dest_base + ext

    if ext in (".cif", ".mmcif", ".mcif"):
        dest_pdb = dest_base + ".pdb"
        if BIOPYTHON_AVAILABLE and MMCIFParser is not None and PDBIO is not None:
            try:
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure("stg", src_path)
                models = list(structure.get_models())
                if not models:
                    raise ValueError("no models in mmCIF")
                first_model_id = models[0].id

                class _FirstModelOnly(Select):
                    def accept_model(self, model):
                        return 1 if model.id == first_model_id else 0

                buf = io.StringIO()
                pdb_out = PDBIO()
                pdb_out.set_structure(structure)
                pdb_out.save(buf, select=_FirstModelOnly())
                header = (
                    "HEADER    STRUCTURE FOR DALILITE/DSSP                      "
                    "01-JAN-00   UNKN\n"
                )
                with open(dest_pdb, "w", encoding="ascii", errors="replace") as fout:
                    fout.write(header)
                    fout.write(buf.getvalue())
                return dest_pdb
            except (OSError, ValueError, KeyError, IndexError, PDBIOException):
                pass
        try:
            shutil.copy2(src_path, dest)
        except OSError:
            return src_path
        return dest

    try:
        with open(src_path, "rb") as f:
            first = f.readline()
    except OSError:
        return src_path
    if not first:
        try:
            shutil.copy2(src_path, dest)
        except OSError:
            return src_path
        return dest
    stripped = first.lstrip()
    if stripped.startswith(b"HEADER") or stripped.startswith(b"data_"):
        try:
            shutil.copy2(src_path, dest)
        except OSError:
            return src_path
        return dest
    header = b"HEADER    STRUCTURE FOR DALILITE/DSSP                      01-JAN-00   UNKN\n"
    try:
        with open(src_path, "rb") as fin, open(dest, "wb") as fout:
            fout.write(header)
            shutil.copyfileobj(fin, fout)
    except OSError:
        return src_path
    return dest


def _dalilite_pair_via_dat(
    pdb_a: str,
    pdb_b: str,
    chain_a: str,
    chain_b: str,
    dalilite_path: str,
    workdir: str,
    outfmt: str,
    title: str,
    dali_timeout: int = 300,
) -> tuple:
    """
    DaliLite: import.pl + dali.pl --cd1/--cd2. Returns (text|None, err, chain_a_res, chain_b_res).
    """
    bin_dir = _find_dalilite_path(dalilite_path)
    if not bin_dir:
        return None, "DaliLite path not set", None, None
    import_pl = os.path.join(bin_dir, "import.pl")
    dali_pl = os.path.join(bin_dir, "dali.pl")
    if not os.path.isfile(import_pl) or not os.path.isfile(dali_pl):
        return None, "import.pl or dali.pl missing in DaliLite bin", None, None

    pre = _dalilite_preflight_mkdssp(bin_dir)
    if pre:
        return None, pre, None, None

    pdb_a = os.path.abspath(pdb_a)
    pdb_b = os.path.abspath(pdb_b)
    if not os.path.isfile(pdb_a) or not os.path.isfile(pdb_b):
        return None, "PDB file missing", None, None

    try:
        run_root = tempfile.mkdtemp(prefix="dli", dir=tempfile.gettempdir())
    except (OSError, FileNotFoundError) as e:
        return None, f"cannot create temp directory for DaliLite: {e}", None, None
    try:
        dat_dir = os.path.join(run_root, "DAT")
        os.makedirs(dat_dir, exist_ok=True)
        qa = _stage_pdb_for_dalilite_import(pdb_a, os.path.join(run_root, "a"))
        qb = _stage_pdb_for_dalilite_import(pdb_b, os.path.join(run_root, "b"))
        pdbid1, pdbid2 = "mol1", "mol2"

        def _pick_chain_after_import(pid: str, expect_chain: str, proc: subprocess.CompletedProcess, pdb_path: str):
            tail = (proc.stderr or proc.stdout or "")[-3000:]
            try:
                names = sorted(os.listdir(dat_dir))
            except OSError:
                names = []
            candidates = []
            for f in names:
                if f.endswith(".dat") and f[:-4].startswith(pid):
                    stem = f[:-4]
                    chain_id = stem[len(pid) :]
                    if chain_id:
                        candidates.append((f, chain_id))
            preferred = f"{pid}{expect_chain}.dat"
            if preferred in [c[0] for c in candidates]:
                return None, expect_chain
            if len(candidates) == 1:
                return None, candidates[0][1]
            if not candidates:
                base = (
                    f"import.pl finished but no {pid}*.dat in {dat_dir}. "
                    f"DAT listing={names!r} import stderr tail={tail!r}"
                )
                base += _mkdssp_diagnostic_suffix(pdb_path, run_root, bin_dir)
                return base, None
            return f"import.pl wrote multiple {pid}*.dat: {[c[0] for c in candidates]!r}; stderr tail={tail!r}", None

        def _import_one(pdb_path: str, pid: str, expect_chain: str):
            cmd = [import_pl, "--pdbfile", pdb_path, "--pdbid", pid, "--dat", dat_dir, "--clean"]
            try:
                proc = subprocess.run(cmd, cwd=run_root, capture_output=True, text=True, timeout=180)
            except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
                return f"import.pl failed: {e}", None
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or "")[-3000:]
                return f"import.pl exit {proc.returncode}: {tail!r}", None
            err, resolved = _pick_chain_after_import(pid, expect_chain, proc, pdb_path)
            return (err or ""), resolved

        err, res_a = _import_one(qa, pdbid1, chain_a)
        if err:
            return None, err, None, None
        err, res_b = _import_one(qb, pdbid2, chain_b)
        if err:
            return None, err, None, None

        chain_a, chain_b = res_a, res_b
        cd1 = f"{pdbid1}{chain_a}"
        cd2 = f"{pdbid2}{chain_b}"
        out_file = os.path.join(run_root, f"{cd1}.txt")
        cmd = [
            dali_pl,
            "--cd1",
            cd1,
            "--cd2",
            cd2,
            "--dat1",
            dat_dir,
            "--dat2",
            dat_dir,
            "--outfmt",
            outfmt,
            "--title",
            title,
            "--clean",
        ]
        try:
            proc = subprocess.run(cmd, cwd=run_root, capture_output=True, text=True, timeout=dali_timeout)
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
            return None, f"dali.pl failed: {e}", None, None
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "")[-4000:]
            return None, f"dali.pl exit {proc.returncode}: {tail!r}", None, None
        if not os.path.isfile(out_file):
            tail = (proc.stderr or proc.stdout or "")[-4000:]
            return None, f"missing output {out_file}. stderr tail: {tail!r}", None, None
        with open(out_file, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
        return content, "", chain_a, chain_b
    finally:
        shutil.rmtree(run_root, ignore_errors=True)


def run_dalilite(pdb_a: str, pdb_b: str, dalilite_path: str | None = None, chain_a: str | None = None, chain_b: str | None = None):
    """
    Run DaliLite pairwise comparison and parse output. Returns dict with summary fields + equivs.
    """
    bin_dir = _find_dalilite_path(dalilite_path)
    if not bin_dir:
        return None
    chain_a = chain_a or "A"
    chain_b = chain_b or "A"
    with tempfile.TemporaryDirectory() as tmpdir:
        content, err, ra, rb = _dalilite_pair_via_dat(
            pdb_a, pdb_b, chain_a, chain_b, dalilite_path or "", tmpdir, "summary,equivalences,alignments", "foldkit", dali_timeout=300
        )
        if content is None:
            return None
    hit = parse_dalilite_summary_hit(content)
    equivs = _parse_dalilite_equivalences(content, ra, rb)
    if not equivs:
        return None
    if hit:
        return {**hit, "equivs": equivs}
    return {"z_score": None, "rmsd": None, "lali": None, "nres": None, "pct_id": None, "hit_id": None, "description": "", "equivs": equivs}


def alignment_from_dalilite(pdb_a: str, pdb_b: str, dalilite_path: str | None = None, chain_a: str | None = None, chain_b: str | None = None):
    """Get residue equivalences from DaliLite. Returns (equivs, dalilite_result) or (None, None)."""
    res = run_dalilite(pdb_a, pdb_b, dalilite_path, chain_a, chain_b)
    if res is None or not res.get("equivs"):
        return None, None
    return res["equivs"], res


def _collect_pdb_files(paths: list, filter_pattern: str | None = None, *, recursive: bool = False) -> list:
    """
    Collect PDB/CIF files from directories and/or files. Optional filter supports glob patterns.
    """
    files: list[str] = []
    seen: set[str] = set()

    def _add(f: str) -> None:
        abspath = os.path.abspath(f)
        if abspath not in seen:
            seen.add(abspath)
            files.append(abspath)

    for p in paths:
        p = os.path.abspath(p)
        if os.path.isfile(p):
            if p.lower().endswith((".pdb", ".cif", ".ent")):
                _add(p)
        elif os.path.isdir(p):
            for ext in ("*.pdb", "*.cif", "*.ent"):
                for f in glob.glob(os.path.join(p, ext)):
                    _add(f)
                if recursive:
                    for f in glob.glob(os.path.join(p, "**", ext), recursive=True):
                        _add(f)

    if filter_pattern:
        def _matches(basename: str, pattern: str) -> bool:
            stem, _ = os.path.splitext(basename)
            if any(ch in pattern for ch in ["*", "?", "["]):
                return fnmatch.fnmatch(basename, pattern) or fnmatch.fnmatch(stem, pattern)
            return pattern in basename or pattern in stem
        files = [f for f in files if _matches(os.path.basename(f), str(filter_pattern))]

    return sorted(set(files))


def _label_from_path(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def _build_tree_from_zscores(zscores: dict, transform: str, exp_scale: float, root: str, midpoint_root: bool) -> str:
    """Build Newick tree from a Z-score dict (DaliLite outputs)."""
    labels, matrix = _zscores_to_distance_matrix(zscores, transform, exp_scale)
    if len(labels) < 2:
        return "();"
    return _build_tree_newick(labels, matrix, root, midpoint_root)


def _zscores_to_distance_matrix(zscores: dict, transform: str, exp_scale: float):
    """Build (labels, matrix) from Z-scores. Z higher = more similar -> distance lower."""
    import math

    labels = sorted(set(a for ab in zscores for a in ab), key=_label_row_sort_key)
    n = len(labels)
    if n == 0:
        return [], []
    seen_obs: set[tuple[str, str]] = set()
    observed: list[float] = []
    for (a, b), v in zscores.items():
        if a == b:
            continue
        ca, cb = _canonical_label_pair(a, b)
        if (ca, cb) in seen_obs:
            continue
        seen_obs.add((ca, cb))
        observed.append(float(v))
    zmax = max(observed) if observed else 0.0
    zmin = min(observed) if observed else 0.0

    def to_dist(z):
        z = max(0.0, float(z))
        if transform == "inv":
            return 1.0 / (1.0 + z)
        if transform == "maxminus":
            return max(0.0, (zmax - z) / max(1e-9, zmax))
        if transform == "exp":
            return math.exp(-z / max(1e-9, exp_scale))
        return 1.0 / (1.0 + z)

    matrix = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            if i != j:
                z = zscores.get((a, b), zmin)
                matrix[i][j] = matrix[j][i] = to_dist(z)
    return labels, matrix


def _build_tree_newick(labels: list, matrix: list, root: str, midpoint_root: bool) -> str:
    """Build Newick tree. Prefer scikit-bio, then Bio.Phylo; fallback to UPGMA."""
    n = len(labels)
    if n < 2:
        return "();"
    try:
        import skbio

        dm = skbio.DistanceMatrix(matrix, labels)
        tree = skbio.tree.nj(dm)
        if root and root in labels:
            tree = tree.root_at({root})
        elif midpoint_root:
            tree = tree.root_at_midpoint()
        return str(tree)
    except ImportError:
        pass
    try:
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
        from io import StringIO
        from Bio import Phylo

        lt = [[matrix[i][j] for j in range(i + 1)] for i in range(len(labels))]
        dm = DistanceMatrix(labels, lt)
        tree = DistanceTreeConstructor().nj(dm)
        if root and root in labels:
            for c in tree.get_terminals():
                if c.name == root:
                    tree.root_with_outgroup(c)
                    break
        s = StringIO()
        Phylo.write(tree, s, "newick")
        return s.getvalue()
    except ImportError:
        pass
    return _upgma_newick_fallback(labels, matrix)


def _upgma_newick_fallback(labels: list, matrix: list) -> str:
    """Pure-Python UPGMA; returns Newick."""
    n = len(labels)
    if n == 0:
        return "();"
    if n == 1:
        return f"({labels[0]});"
    if n == 2:
        d = float(matrix[0][1])
        return f"({labels[0]}:{d/2:.6f},{labels[1]}:{d/2:.6f});"
    cluster_size = [1] * n + [0] * (n - 1)
    height = [0.0] * (2 * n - 1)
    parent = [None] * (2 * n - 1)
    dist = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = float(matrix[i][j])
    next_node = n
    active = set(range(n))
    for _ in range(n - 1):
        best, best_pair = float("inf"), None
        for i in active:
            for j in active:
                if i >= j:
                    continue
                d = dist.get((min(i, j), max(i, j)), float("inf"))
                if d < best:
                    best, best_pair = d, (i, j)
        if best_pair is None:
            break
        i, j = best_pair
        ni, nj = cluster_size[i], cluster_size[j]
        h = best / 2.0
        height[next_node] = h
        cluster_size[next_node] = ni + nj
        parent[i] = parent[j] = next_node
        active.discard(i)
        active.discard(j)
        active.add(next_node)
        for k in active:
            if k == next_node:
                continue
            ki, kj = (min(i, k), max(i, k)), (min(j, k), max(j, k))
            d_ik = dist.get(ki, 0.0)
            d_jk = dist.get(kj, 0.0)
            new_d = (d_ik * ni + d_jk * nj) / (ni + nj)
            dist[(min(k, next_node), max(k, next_node))] = new_d
        next_node += 1

    def newick_node(node):
        if node < n:
            return labels[node]
        children = [c for c in range(2 * n - 1) if parent[c] == node]
        if len(children) != 2:
            return "();"
        L, R = children[0], children[1]
        h = height[node]
        bl_L = h - (height[L] if L >= n else 0.0)
        bl_R = h - (height[R] if R >= n else 0.0)
        return f"({newick_node(L)}:{bl_L:.6f},{newick_node(R)}:{bl_R:.6f})"

    return newick_node(next_node - 1) + ";"


def _zscore_for_undirected_pair(zscores: dict, a: str, b: str) -> float | None:
    t = zscores.get((a, b))
    if t is not None:
        return float(t)
    t = zscores.get((b, a))
    if t is not None:
        return float(t)
    return None


def _rank_structures(zscores: dict) -> list:
    labels = sorted(set(a for ab in zscores for a in ab))
    out = []
    for a in labels:
        vals: list[float] = []
        for b in labels:
            if b == a:
                continue
            z = _zscore_for_undirected_pair(zscores, a, b)
            if z is not None:
                vals.append(z)
        out.append({"label": a, "n_pairs": len(vals), "avg_z": sum(vals) / len(vals) if vals else 0.0, "max_z": max(vals) if vals else 0.0})
    out.sort(key=lambda r: (-r["avg_z"], -r["max_z"], r["label"]))
    return out


def _plot_tree(newick: str, out_path: str, title: str) -> bool:
    try:
        from ete3 import Tree, TreeStyle, TextFace

        t = Tree(newick)
        ts = TreeStyle()
        ts.title.add_face(TextFace(title, fsize=12), column=0)
        t.render(out_path, tree_style=ts)
        return True
    except ImportError:
        pass
    try:
        from Bio import Phylo
        from io import StringIO
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        tree = Phylo.read(StringIO(newick), "newick")
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title(title)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        return True
    except ImportError:
        return False

