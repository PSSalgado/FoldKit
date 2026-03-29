#!/usr/bin/env python3
from __future__ import annotations

"""
Dali-like structural similarity score
====================================

Compute a Dali-equivalent score from two superposed/embedded structures using
residue equivalences. The score is based on intramolecular Cα–Cα distance
similarities (Holm & Sander, J. Mol. Biol. 1993) and does not require
superposition—only a set of equivalent residue pairs.

Input: two PDB files + residue equivalences (from alignment file, DaliLite,
biotite, or sequence-order fallback).

When DaliLite is available (path specified by user via --dalilite-path or
DALILITE_HOME), it is used as the primary alignment source—it provides the
canonical Dali Z-score and residue equivalences.

DaliLite's import.pl calls mkdssp via the path in ``bin/mpidali.pm`` (``$DSSP_EXE``,
often ``/usr/local/bin/mkdssp`` in upstream installs). That binary must exist and
be executable or imports produce empty ``DAT/*.dat``. Optionally set environment
variable ``MKDSSP`` to the full path of mkdssp for diagnostics in this script.

Formula (Dali raw score):
  S(A,B) = Σ_{i,j in core} (θ - |d_ij^A - d_ij^B|/d*_ij) * exp(-(d*_ij/20)²)
  where θ=0.2, d*_ij = (d_ij^A + d_ij^B)/2, and d_ij are intramolecular distances.

Z-score: Z = (S - m(L)) / σ(L), with L = √(L_A * L_B).
m(L) and σ(L) from empirical fits (Holm & Sander).
"""

import argparse
import csv
import fnmatch
import glob
import os
import shutil
import re
import subprocess
import sys
import tempfile
import warnings

try:
    import numpy as np
except ImportError:
    np = None

try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.simplefilter('ignore', PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    PDBParser = None
    PDBConstructionWarning = None
    BIOPYTHON_AVAILABLE = False

try:
    from biotite.structure.io.pdb import PDBFile
    from biotite.structure import filter_backbone
    from biotite.structure.superimpose import superimpose_structural_homologs
    BIOTITE_AVAILABLE = True
except ImportError:
    PDBFile = None
    filter_backbone = None
    superimpose_structural_homologs = None
    BIOTITE_AVAILABLE = False

# Dali constants
THETA = 0.2
R0 = 20.0


def _mean_score_fit(L: float) -> float:
    """
    Empirical mean score m(L) for random pairs (Holm & Sander).
    Equation 5: m(L) ~ 7.95 + 0.71*L - 2.59e-4*L² - 1.92e-6*L³ if L ≤ 400.
    m(L) = m(400) + (L - 400) if L > 400.
    """
    if L <= 400:
        return 7.95 + 0.71 * L - 2.59e-4 * (L ** 2) - 1.92e-6 * (L ** 3)
    return _mean_score_fit(400) + (L - 400)


def _std_score_fit(L: float) -> float:
    """Empirical standard deviation: σ(L) = 0.5 * m(L)."""
    return 0.5 * _mean_score_fit(L)


def _residue_pair_score(d_A: float, d_B: float) -> float:
    """
    Dali residue-pair score φ(i,j).
    φ(i,j) = (θ - diff(i,j)) * exp(-(d*_ij/20)²)
    diff = |d_ij^A - d_ij^B| / d*_ij, d*_ij = (d_A + d_B) / 2
    """
    if d_A <= 0 and d_B <= 0:
        return 0.0
    d_star = (d_A + d_B) / 2.0
    if d_star <= 0:
        return 0.0
    diff = abs(d_A - d_B) / d_star
    env = np.exp(-(d_star / R0) ** 2)
    raw = THETA - diff
    return max(0.0, raw * env)


def get_ca_coords_biopython(structure, chain_id: str = None):
    """
    Extract (residue_id, CA_coord) from a BioPython structure.
    residue_id = (chain_id, resseq, icode).
    Returns dict: residue_id -> np.array([x,y,z])
    """
    if not BIOPYTHON_AVAILABLE or PDBParser is None:
        return {}
    coords = {}
    for model in structure:
        for chain in model:
            cid = chain.id
            if chain_id is not None and cid != chain_id:
                continue
            for res in chain:
                if res.id[0] != ' ':  # skip hetero
                    continue
                try:
                    ca = res['CA']
                except KeyError:
                    continue
                key = (cid, res.id[1], res.id[2] if len(res.id) > 2 else ' ')
                coords[key] = np.array(ca.coord)
    return coords


def parse_alignment_file(path: str):
    """
    Parse alignment file. Supported formats:
    - TSV/CSV: resnum_A, resnum_B (for single chain per structure)
    - TSV/CSV: chain_A, resnum_A, chain_B, resnum_B
    - Lines starting with # are comments.
    - Header row is auto-detected if first line looks like a header.

    Returns list of ((chain_A, resnum_A, icode_A), (chain_B, resnum_B, icode_B)).
    """
    equivs = []
    with open(path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or (len(row) > 0 and str(row[0]).strip().startswith('#')):
                continue
            # Skip header-like row
            if len(equivs) == 0 and len(row) >= 2:
                first = str(row[0]).strip().upper()
                if first in ('CHAIN_A', 'CHAIN', 'RESNUM_A', 'RESNUM', 'RESA', 'RESB'):
                    continue
            try:
                if len(row) >= 4:
                    cA, rA, cB, rB = row[0], row[1], row[2], row[3]
                    try:
                        resnumA = int(rA)
                        resnumB = int(rB)
                    except ValueError:
                        continue
                    equivs.append(((str(cA).strip(), resnumA, ' '), (str(cB).strip(), resnumB, ' ')))
                elif len(row) >= 2:
                    rA, rB = row[0], row[1]
                    try:
                        resnumA = int(rA)
                        resnumB = int(rB)
                    except ValueError:
                        continue
                    equivs.append(((' ', resnumA, ' '), (' ', resnumB, ' ')))
            except (ValueError, IndexError):
                continue
    return equivs


def normalize_equivalences(equivs, coords_A, coords_B):
    """
    Convert equivalences to keys that exist in coords_A and coords_B.
    For single-chain: (chain, resnum, icode) may use ' ' for chain; match by resnum.
    """
    by_resnum_A = {}
    for k, v in coords_A.items():
        chain, resnum, icode = k
        by_resnum_A[(chain, resnum)] = k
        by_resnum_A[resnum] = k  # fallback for alignment with resnum only
    by_resnum_B = {}
    for k, v in coords_B.items():
        chain, resnum, icode = k
        by_resnum_B[(chain, resnum)] = k
        by_resnum_B[resnum] = k

    out = []
    for eA, eB in equivs:
        cA, rA, iA = eA
        cB, rB, iB = eB
        keyA = (cA, rA) if cA and cA != ' ' else rA
        keyB = (cB, rB) if cB and cB != ' ' else rB
        kA = by_resnum_A.get(keyA)
        kB = by_resnum_B.get(keyB)
        if kA is not None and kB is not None:
            out.append((kA, kB))
    return out


def _find_dalilite_path(dalilite_path: str = None) -> str:
    """Return path to DaliLite bin directory, or empty string if not found.
    Path must be specified by user via --dalilite-path or DALILITE_HOME env.
    """
    candidates = [p for p in (dalilite_path, os.environ.get('DALILITE_HOME', '')) if p]
    for base in candidates:
        for check in (base, os.path.join(base, 'bin')):
            dali_pl = os.path.join(check, 'dali.pl')
            if os.path.isfile(dali_pl):
                return os.path.dirname(dali_pl)
    return ''


def _parse_dalilite_equivalences(txt: str, chain_a: str = 'A', chain_b: str = 'A'):
    """
    Parse DaliLite structural equivalences block.
    Formats:
      "   1: query-A target-A     1 -  33     1 -  33   (...)"
      v5: "   1: mol1-A mol2-A     6 -  29 <=>   13 -  36   (...)"
    Returns list of ((chain_a, resnum_a, ' '), (chain_b, resnum_b, ' ')).
    """
    equivs = []
    # Two range pairs separated by whitespace or by " <=> " (DaliLite v5).
    pattern = re.compile(
        r'^\s*\d+:\s+\S+\s+\S+\s+(\d+)\s*-\s*(\d+)(?:\s*<=\>\s*|\s+)(\d+)\s*-\s*(\d+)',
        re.MULTILINE,
    )
    for m in pattern.finditer(txt):
        q_start, q_end = int(m.group(1)), int(m.group(2))
        t_start, t_end = int(m.group(3)), int(m.group(4))
        length = q_end - q_start + 1
        t_length = t_end - t_start + 1
        n = min(length, t_length)
        for i in range(n):
            equivs.append(
                ((chain_a, q_start + i, ' '), (chain_b, t_start + i, ' '))
            )
    return equivs


def _dalilite_text_has_summary_hit(txt: str) -> bool:
    """
    True if DaliLite summary block contains at least one hit line (Z, rmsd, lali, ...).
    When False, equivalences/transrot are usually empty—DaliLite omits alignments for
    non-significant pairs (commonly Z < ~2).
    """
    return bool(
        re.search(
            r'^\s*\d+:\s+\S+\s+([-\d.]+)\s+([\d.]+)\s+(\d+)\s+\d+\s+\d+',
            txt,
            re.MULTILINE,
        )
    )


def _get_first_chain_id(coords_dict) -> str:
    """Return the first chain ID seen in coords, or 'A'."""
    if not coords_dict:
        return 'A'
    for k in coords_dict:
        if k[0] and k[0] != ' ':
            return str(k[0])
    return 'A'


def _read_dalilite_dssp_exe_from_mpidali(bin_dir: str) -> str:
    """Return ``$DSSP_EXE`` path from DaliLite's mpidali.pm, or ''."""
    pm = os.path.join(bin_dir, 'mpidali.pm')
    try:
        with open(pm, 'r', encoding='utf-8', errors='replace') as f:
            head = f.read(20000)
    except OSError:
        return ''
    m = re.search(r'^\s*my\s+\$DSSP_EXE\s*=\s*"([^"]+)"\s*;', head, re.MULTILINE)
    return m.group(1) if m else ''


def _dalilite_preflight_mkdssp(bin_dir: str) -> str | None:
    """
    If mpidali.pm names an mkdssp path that is missing, return an error message.
    Otherwise return None.
    """
    exe = _read_dalilite_dssp_exe_from_mpidali(bin_dir)
    if not exe:
        return None
    if os.path.isfile(exe) and os.access(exe, os.X_OK):
        return None
    return (
        'DaliLite import.pl requires mkdssp (see bin/mpidali.pm $DSSP_EXE). '
        f'Expected {exe!r} but it is missing or not executable. Install a '
        'compatible mkdssp, symlink it to that path, or edit $DSSP_EXE in '
        'mpidali.pm to your mkdssp binary.'
    )


def _find_mkdssp_executable(dalilite_bin_dir: str | None = None) -> str:
    """Path to mkdssp for diagnostics (prefer MKDSSP env, then DaliLite mpidali.pm)."""
    env_mk = (os.environ.get('MKDSSP') or '').strip()
    if env_mk and os.path.isfile(env_mk) and os.access(env_mk, os.X_OK):
        return env_mk
    if dalilite_bin_dir:
        exe = _read_dalilite_dssp_exe_from_mpidali(dalilite_bin_dir)
        if exe and os.path.isfile(exe) and os.access(exe, os.X_OK):
            return exe
    mk = shutil.which('mkdssp')
    if mk:
        return mk
    ccp4 = (os.environ.get('CCP4') or '').strip()
    if ccp4:
        c = os.path.join(ccp4, 'bin', 'mkdssp')
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    for cand in ('/usr/local/bin/mkdssp', '/usr/bin/mkdssp'):
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return ''


def _mkdssp_diagnostic_suffix(
    pdb_path: str, workdir: str, dalilite_bin_dir: str | None = None,
) -> str:
    """
    When import.pl leaves DAT/*.dat empty, run mkdssp the same way DaliLite often does
    (two-arg: pdb out.dssp) and append the real error to the user message.
    """
    mk = _find_mkdssp_executable(dalilite_bin_dir)
    if not mk:
        hint = _read_dalilite_dssp_exe_from_mpidali(dalilite_bin_dir) if dalilite_bin_dir else ''
        extra = f' mpidali.pm expects {hint!r}.' if hint else ''
        return (
            ' | mkdssp probe: no usable mkdssp (set MKDSSP to the binary path, '
            'or install mkdssp and ensure bin/mpidali.pm $DSSP_EXE points to it).'
            + extra
        )
    try:
        fd, out = tempfile.mkstemp(suffix='.dssp', dir=workdir)
        os.close(fd)
    except (OSError, FileNotFoundError) as e:
        return f' | mkdssp probe: could not create temp dssp in {workdir!r}: {e}'
    try:
        proc = subprocess.run(
            [mk, pdb_path, out],
            capture_output=True,
            text=True,
            timeout=120,
        )
        sz = os.path.getsize(out) if os.path.isfile(out) else 0
        tail_err = (proc.stderr or '')[-2500:]
        tail_out = (proc.stdout or '')[-800:]
        if proc.returncode == 0 and sz > 0:
            return (
                f' | mkdssp probe OK ({mk} wrote {sz} bytes) but import.pl still produced '
                f'no .dat—check DaliLite install, permissions, and import.pl logs.'
            )
        return (
            f' | mkdssp probe ({mk}): exit={proc.returncode} dssp_bytes={sz} '
            f'stderr_tail={tail_err!r} stdout_tail={tail_out!r}'
        )
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        return f' | mkdssp probe failed to run {mk!r}: {e}'
    finally:
        try:
            os.remove(out)
        except OSError:
            pass


def _stage_pdb_for_dalilite_import(src_path: str, dest_base: str) -> str:
    """
    Copy or rewrite structure into dest_base + extension for DaliLite import.

    mkdssp (inside import.pl) treats files without a PDB HEADER as mmCIF; prepend
    HEADER when needed. mmCIF is copied as-is.
    """
    ext = os.path.splitext(src_path)[1].lower()
    if ext not in ('.cif', '.mmcif', '.mcif', '.pdb', '.ent', ''):
        ext = '.pdb'
    if not ext:
        ext = '.pdb'
    dest = dest_base + ext

    if ext in ('.cif', '.mmcif', '.mcif'):
        try:
            shutil.copy2(src_path, dest)
        except OSError:
            return src_path
        return dest

    try:
        with open(src_path, 'rb') as f:
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
    if stripped.startswith(b'HEADER') or stripped.startswith(b'data_'):
        try:
            shutil.copy2(src_path, dest)
        except OSError:
            return src_path
        return dest
    header = (
        b'HEADER    STRUCTURE FOR DALILITE/DSSP                      01-JAN-00   UNKN\n'
    )
    try:
        with open(src_path, 'rb') as fin, open(dest, 'wb') as fout:
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
    DaliLite: import.pl + dali.pl --cd1/--cd2.

    Runs under a short-lived directory under the system temp folder (prefix ``dli``):
    DaliLite Fortran enforces paths <= 80 characters ("STOP: file names incl. path
    must be max 80 long"), so very long input paths must not be passed to import.pl.
    ``workdir`` is kept for API compatibility; DAT and outputs use the temp run root.

    Returns (output_text_or_None, error_message, chain_a_resolved, chain_b_resolved).
    On failure: (None, msg, None, None). On success: (text, '', ca, cb) for Dali cd ids.
    """
    bin_dir = _find_dalilite_path(dalilite_path)
    if not bin_dir:
        return None, 'DaliLite path not set', None, None
    import_pl = os.path.join(bin_dir, 'import.pl')
    dali_pl = os.path.join(bin_dir, 'dali.pl')
    if not os.path.isfile(import_pl) or not os.path.isfile(dali_pl):
        return None, 'import.pl or dali.pl missing in DaliLite bin', None, None

    pre = _dalilite_preflight_mkdssp(bin_dir)
    if pre:
        return None, pre, None, None

    pdb_a = os.path.abspath(pdb_a)
    pdb_b = os.path.abspath(pdb_b)
    if not os.path.isfile(pdb_a) or not os.path.isfile(pdb_b):
        return None, 'PDB file missing', None, None

    try:
        run_root = tempfile.mkdtemp(prefix='dli', dir=tempfile.gettempdir())
    except (OSError, FileNotFoundError) as e:
        return None, f'cannot create temp directory for DaliLite: {e}', None, None

    try:
        dat_dir = os.path.join(run_root, 'DAT')
        os.makedirs(dat_dir, exist_ok=True)

        qa = _stage_pdb_for_dalilite_import(pdb_a, os.path.join(run_root, 'a'))
        qb = _stage_pdb_for_dalilite_import(pdb_b, os.path.join(run_root, 'b'))

        pdbid1, pdbid2 = 'mol1', 'mol2'

        def _pick_chain_after_import(
            pid: str,
            expect_chain: str,
            proc: subprocess.CompletedProcess,
            pdb_path: str = '',
        ):
            """
            import.pl may write {pid}{chain}.dat where chain != expect_chain (e.g. B).
            If expect_chain file missing, use sole {pid}*.dat or match expect_chain.
            Returns (error_message_or_None, resolved_chain_or_None).
            """
            tail = (proc.stderr or proc.stdout or '')[-3000:]
            try:
                names = sorted(os.listdir(dat_dir))
            except OSError:
                names = []
            candidates = []
            for f in names:
                if not f.endswith('.dat'):
                    continue
                stem = f[:-4]
                if not stem.startswith(pid):
                    continue
                chain_id = stem[len(pid):]
                if not chain_id:
                    continue
                candidates.append((f, chain_id))
            preferred = f'{pid}{expect_chain}.dat'
            if preferred in [c[0] for c in candidates]:
                return None, expect_chain
            if len(candidates) == 1:
                return None, candidates[0][1]
            for _f, cid in candidates:
                if cid == expect_chain:
                    return None, expect_chain
            if not candidates:
                base = (
                    'import.pl finished but no {0}*.dat in {1}. DaliLite needs DSSP; '
                    'mkdssp usually fails on CA-only PDBs (needs N, CA, C, O per residue). '
                    'Use full-backbone PDB/mmCIF (e.g. mmCIF with all backbone atoms), or '
                    'rebuild backbone before import. DAT listing={2!r} import stderr tail={3!r}'
                ).format(pid, dat_dir, names, tail)
                if pdb_path:
                    base += _mkdssp_diagnostic_suffix(pdb_path, run_root, bin_dir)
                return base, None
            return (
                'import.pl wrote multiple {0}*.dat: {1!r}; set --chain-a/--chain-b to match '
                'PDB chain id. stderr tail={2!r}'
            ).format(pid, [c[0] for c in candidates], tail), None

        def _import_one(pdb_path: str, pid: str, expect_chain: str):
            cmd = [
                import_pl, '--pdbfile', pdb_path, '--pdbid', pid,
                '--dat', dat_dir, '--clean',
            ]
            try:
                proc = subprocess.run(
                    cmd, cwd=run_root, capture_output=True, text=True, timeout=180,
                )
            except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
                return f'import.pl failed: {e}', None
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or '')[-3000:]
                return f'import.pl exit {proc.returncode}: {tail!r}', None
            err, resolved = _pick_chain_after_import(pid, expect_chain, proc, pdb_path)
            if err:
                return err, None
            return '', resolved

        err, res_a = _import_one(qa, pdbid1, chain_a)
        if err:
            return None, err, None, None
        err, res_b = _import_one(qb, pdbid2, chain_b)
        if err:
            return None, err, None, None

        chain_a, chain_b = res_a, res_b
        cd1 = f'{pdbid1}{chain_a}'
        cd2 = f'{pdbid2}{chain_b}'
        out_file = os.path.join(run_root, f'{cd1}.txt')
        cmd = [
            dali_pl,
            '--cd1', cd1,
            '--cd2', cd2,
            '--dat1', dat_dir,
            '--dat2', dat_dir,
            '--outfmt', outfmt,
            '--title', title,
            '--clean',
        ]
        try:
            proc = subprocess.run(
                cmd, cwd=run_root, capture_output=True, text=True, timeout=dali_timeout,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
            return None, f'dali.pl failed: {e}', None, None

        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or '')[-4000:]
            return None, f'dali.pl exit {proc.returncode}: {tail!r}', None, None

        if not os.path.isfile(out_file):
            tail = (proc.stderr or proc.stdout or '')[-4000:]
            return None, f'missing output {out_file}. stderr tail: {tail!r}', None, None

        with open(out_file, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
        return content, '', chain_a, chain_b
    finally:
        shutil.rmtree(run_root, ignore_errors=True)


def run_dalilite(pdb_a: str, pdb_b: str, dalilite_path: str = None,
                 chain_a: str = None, chain_b: str = None):
    """
    Run DaliLite pairwise comparison and parse output.
    Returns dict with z_score, rmsd, lali, equivs, or None on failure.

    Uses import.pl to build DAT/*.dat then dali.pl --cd1/--cd2 (required for
    typical DaliLite installs; bare --pdbfile1/2 often yields empty DAT).
    """
    bin_dir = _find_dalilite_path(dalilite_path)
    if not bin_dir:
        return None

    pdb_a = os.path.abspath(pdb_a)
    pdb_b = os.path.abspath(pdb_b)
    if not os.path.isfile(pdb_a) or not os.path.isfile(pdb_b):
        return None

    chain_a = chain_a or 'A'
    chain_b = chain_b or 'A'

    with tempfile.TemporaryDirectory() as tmpdir:
        content, err, ra, rb = _dalilite_pair_via_dat(
            pdb_a, pdb_b, chain_a, chain_b, dalilite_path, tmpdir,
            'summary,equivalences,alignments', 'dali_score', dali_timeout=300,
        )
        if content is None:
            return None

    # Parse summary: "   1:  xxxxX  Z  rmsd  lali  nres  ..."
    # Fields: No, ids, Z, rmsd, lali, nres, %id, Description
    z_score = None
    rmsd = None
    lali = None
    summary_pattern = re.compile(
        r'^\s*\d+:\s+\S+\s+([-\d.]+)\s+([\d.]+)\s+(\d+)\s+\d+\s+\d+',
        re.MULTILINE
    )
    m = summary_pattern.search(content)
    if m:
        z_score = float(m.group(1))
        rmsd = float(m.group(2))
        lali = int(m.group(3))

    equivs = _parse_dalilite_equivalences(content, ra, rb)
    if not equivs:
        return None

    return {
        'z_score': z_score,
        'rmsd': rmsd,
        'lali': lali,
        'equivs': equivs,
    }


def alignment_from_dalilite(pdb_a: str, pdb_b: str, dalilite_path: str = None,
                            chain_a: str = None, chain_b: str = None):
    """
    Get residue equivalences from DaliLite. Returns (equivs, dalilite_result)
    or (None, None) on failure.
    """
    res = run_dalilite(pdb_a, pdb_b, dalilite_path, chain_a, chain_b)
    if res is None or not res.get('equivs'):
        return None, None
    return res['equivs'], res


def alignment_from_biotite(pdb_a: str, pdb_b: str, chain_a: str = None, chain_b: str = None):
    """
    Compute residue equivalences using biotite's superimpose_structural_homologs.
    Returns list of (keyA, keyB) where keys are (chain_id, resseq, icode).
    """
    if not BIOTITE_AVAILABLE:
        return None
    try:
        struct_a = PDBFile.read(pdb_a)
        struct_b = PDBFile.read(pdb_b)
        a = struct_a.get_structure()[0]
        b = struct_b.get_structure()[0]
        ca_a = a[a.atom_name == 'CA']
        ca_b = b[b.atom_name == 'CA']
        if chain_a is not None:
            ca_a = ca_a[ca_a.chain_id == chain_a]
        if chain_b is not None:
            ca_b = ca_b[ca_b.chain_id == chain_b]
        if len(ca_a) < 3 or len(ca_b) < 3:
            return None
        _, _, fix_idx, mob_idx = superimpose_structural_homologs(ca_a, ca_b)
        equivs = []
        fix_idx = np.atleast_1d(np.asarray(fix_idx))
        mob_idx = np.atleast_1d(np.asarray(mob_idx))
        for idx in range(len(fix_idx)):
            ia, ib = int(fix_idx[idx]), int(mob_idx[idx])
            chain_a_val = str(ca_a.chain_id[ia])
            chain_b_val = str(ca_b.chain_id[ib])
            res_a = ca_a.res_id[ia]
            res_b = ca_b.res_id[ib]
            resseq_a = int(res_a) if np.isscalar(res_a) else int(res_a[0])
            resseq_b = int(res_b) if np.isscalar(res_b) else int(res_b[0])
            icode_a = ' '
            icode_b = ' '
            if 'ins_code' in ca_a.get_annotation_categories():
                icode_a = str(ca_a.ins_code[ia]) if ca_a.ins_code[ia] else ' '
            if 'ins_code' in ca_b.get_annotation_categories():
                icode_b = str(ca_b.ins_code[ib]) if ca_b.ins_code[ib] else ' '
            key_a = (chain_a_val, resseq_a, icode_a)
            key_b = (chain_b_val, resseq_b, icode_b)
            equivs.append((key_a, key_b))
        return equivs
    except Exception:
        return None


def alignment_from_sequence_order(coords_A, coords_B, chain_a: str = None, chain_b: str = None):
    """
    Fallback: assume 1:1 correspondence by residue number when both structures
    have the same set of residue numbers. Only for single-chain or when chain
    specified and residues share numbering.
    """
    # Filter by chain if specified
    keys_a = sorted([k for k in coords_A if (chain_a is None or k[0] == chain_a)])
    keys_b = sorted([k for k in coords_B if (chain_b is None or k[0] == chain_b)])
    # Match by (chain, resnum)
    set_a = set((k[0], k[1]) for k in keys_a)
    set_b = set((k[0], k[1]) for k in keys_b)
    common = set_a & set_b
    if len(common) < 3:
        return None
    key_map_a = {(k[0], k[1]): k for k in keys_a}
    key_map_b = {(k[0], k[1]): k for k in keys_b}
    equivs = []
    for c, r in sorted(common):
        equivs.append((key_map_a[(c, r)], key_map_b[(c, r)]))
    return equivs


def compute_dali_score(coords_A, coords_B, equivs):
    """
    Compute raw Dali score S(A,B) from intramolecular Cα–Cα distances.
    equivs: list of (keyA, keyB) where keys index into coords_A and coords_B.

    Returns (raw_score: float, n_core: int). Requires NumPy; callers must check
    ``np is not None`` before calling (``run()`` does this).
    """
    if np is None:
        raise RuntimeError('NumPy is required for compute_dali_score')
    core = equivs
    n = len(core)
    if n < 2:
        return 0.0, n

    # Build index list for core residues (order matters for pair enumeration)
    idx_to_key_a = [c[0] for c in core]
    idx_to_key_b = [c[1] for c in core]

    # Precompute all pairwise distances within each structure
    def dist_matrix(keys, coords):
        m = len(keys)
        D = np.zeros((m, m))
        for i in range(m):
            for j in range(i + 1, m):
                d = np.linalg.norm(coords[keys[i]] - coords[keys[j]])
                D[i, j] = D[j, i] = d
        return D

    D_A = dist_matrix(idx_to_key_a, coords_A)
    D_B = dist_matrix(idx_to_key_b, coords_B)

    S = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d_A = D_A[i, j]
            d_B = D_B[i, j]
            phi = _residue_pair_score(d_A, d_B)
            S += phi
    # Dali sums over unordered pairs; each pair (i,j) with i<j is counted once,
    # but the formula often writes Σ Σ. Check original: it's a double sum over
    # i,j in core, so each pair (i,j) and (j,i) both contribute. Typically
    # φ(i,j)=φ(j,i), so S = 2 * sum_{i<j} φ(i,j). Our loop counts both,
    # so we're good.
    return float(S), n


def _collect_pdb_files(paths: list, filter_pattern: str = None) -> list:
    """Collect PDB/CIF files from directories and/or files. Optional --filter."""
    files = []
    seen = set()

    def _add(f):
        abspath = os.path.abspath(f)
        if abspath not in seen:
            seen.add(abspath)
            files.append(abspath)

    for p in paths:
        p = os.path.abspath(p)
        if os.path.isfile(p):
            if p.lower().endswith(('.pdb', '.cif', '.ent')):
                _add(p)
        elif os.path.isdir(p):
            for ext in ('*.pdb', '*.cif', '*.ent'):
                for f in glob.glob(os.path.join(p, ext)):
                    _add(f)
                for f in glob.glob(os.path.join(p, '**', ext), recursive=True):
                    _add(f)

    if filter_pattern:
        def _matches(basename, pattern):
            stem, _ = os.path.splitext(basename)
            if any(ch in pattern for ch in ['*', '?', '[']):
                return fnmatch.fnmatch(basename, pattern) or fnmatch.fnmatch(stem, pattern)
            return pattern in basename or pattern in stem

        files = [f for f in files if _matches(os.path.basename(f), filter_pattern)]

    return sorted(set(files))


def _label_from_path(path: str) -> str:
    """Use stem (no extension) as structure label."""
    return os.path.splitext(os.path.basename(path))[0]


def compute_z_score(raw_score: float, L_A: int, L_B: int) -> float:
    """Z-score = (S - m(L)) / σ(L), L = sqrt(L_A * L_B)."""
    if raw_score is None:
        raise TypeError('raw_score must be numeric, not None')
    L = (L_A * L_B) ** 0.5
    m = _mean_score_fit(L)
    sigma = _std_score_fit(L)
    if sigma <= 0:
        return 0.0
    return (raw_score - m) / sigma


def run(pdb_a: str, pdb_b: str, alignment_file: str = None,
        chain_a: str = None, chain_b: str = None,
        use_dalilite: bool = True, use_biotite: bool = True,
        use_sequence_order: bool = True, dalilite_path: str = None):
    """
    Compute Dali-like score and Z-score for structure pair.

    Returns:
        dict with keys: raw_score, z_score, n_core, L_A, L_B, alignment_source,
        plus dalilite_z_score, dalilite_rmsd when DaliLite is used.
    """
    if not BIOPYTHON_AVAILABLE or PDBParser is None:
        return {'error': 'BioPython is required'}
    if np is None:
        return {
            'error': 'NumPy is required for Dali score computation',
            'hint': 'pip install numpy',
        }

    parser = PDBParser(QUIET=True)
    try:
        struct_a = parser.get_structure('A', pdb_a)
        struct_b = parser.get_structure('B', pdb_b)
    except Exception as e:
        return {'error': f'Failed to parse PDB files: {e}'}

    coords_A = get_ca_coords_biopython(struct_a, chain_a)
    coords_B = get_ca_coords_biopython(struct_b, chain_b)

    if not coords_A or not coords_B:
        return {'error': 'No CA atoms found in one or both structures'}

    L_A = len(coords_A)
    L_B = len(coords_B)
    equivs = None
    source = None
    dalilite_result = None

    # 1. DaliLite (canonical Dali implementation)
    if use_dalilite:
        cid_a = chain_a or _get_first_chain_id(coords_A)
        cid_b = chain_b or _get_first_chain_id(coords_B)
        equivs_dali, dalilite_result = alignment_from_dalilite(
            pdb_a, pdb_b, dalilite_path, cid_a, cid_b
        )
        if equivs_dali:
            equivs = normalize_equivalences(
                equivs_dali, coords_A, coords_B
            )
            if equivs:
                source = 'dalilite'

    # 2. From alignment file
    if (equivs is None or len(equivs) < 3) and alignment_file and os.path.isfile(alignment_file):
        raw_equivs = parse_alignment_file(alignment_file)
        equivs = normalize_equivalences(raw_equivs, coords_A, coords_B)
        if equivs:
            source = 'alignment_file'

    # 2. Biotite structural alignment
    if (equivs is None or len(equivs) < 3) and use_biotite and BIOTITE_AVAILABLE:
        equivs_bio = alignment_from_biotite(pdb_a, pdb_b, chain_a, chain_b)
        if equivs_bio:
            equivs = [(kA, kB) for kA, kB in equivs_bio
                      if kA in coords_A and kB in coords_B]
            if equivs:
                source = 'biotite'

    # 3. Sequence-order fallback (same residue numbering)
    if (equivs is None or len(equivs) < 3) and use_sequence_order:
        equivs = alignment_from_sequence_order(coords_A, coords_B, chain_a, chain_b)
        if equivs:
            source = 'sequence_order'

    if equivs is None or len(equivs) < 2:
        return {
            'error': 'Could not determine residue equivalences',
            'hint': 'Use DaliLite (--dalilite-path), provide --alignment, or install biotite',
        }

    raw_score, n_core = compute_dali_score(coords_A, coords_B, equivs)

    # Use DaliLite Z-score when available (canonical); else our empirical fit
    if dalilite_result and dalilite_result.get('z_score') is not None:
        z_score = dalilite_result['z_score']
    else:
        z_score = compute_z_score(raw_score, L_A, L_B)

    out = {
        'raw_score': raw_score,
        'z_score': z_score,
        'n_core': n_core,
        'L_A': L_A,
        'L_B': L_B,
        'alignment_source': source or 'unknown',
    }
    if dalilite_result:
        out['dalilite_z_score'] = dalilite_result.get('z_score')
        out['dalilite_rmsd'] = dalilite_result.get('rmsd')
    return out


def run_all_vs_all(files: list, chain_a: str = None, chain_b: str = None,
                   use_dalilite: bool = True, use_biotite: bool = True,
                   use_sequence_order: bool = True, dalilite_path: str = None,
                   verbose: bool = True) -> dict:
    """
    Run pairwise Dali comparison for all pairs. Returns dict:
      zscores: (label_a, label_b) -> z_score
      raw_scores: (label_a, label_b) -> raw_score
      n_core: (label_a, label_b) -> n_core
      results: full run() results per pair
    """
    labels = [_label_from_path(f) for f in files]
    zscores = {}
    raw_scores = {}
    n_core = {}
    n_pairs = len(files) * (len(files) - 1) // 2
    done = 0

    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            pa, pb = files[i], files[j]
            la, lb = labels[i], labels[j]
            result = run(pa, pb, chain_a=chain_a, chain_b=chain_b,
                        use_dalilite=use_dalilite, use_biotite=use_biotite,
                        use_sequence_order=use_sequence_order,
                        dalilite_path=dalilite_path)
            done += 1
            if verbose and n_pairs > 1:
                print(f"[{done}/{n_pairs}] {la} vs {lb}: Z={result.get('z_score', 'N/A')}", file=sys.stderr)
            if 'error' in result:
                continue
            z = result['z_score']
            zscores[(la, lb)] = zscores[(lb, la)] = z
            raw_scores[(la, lb)] = raw_scores[(lb, la)] = result['raw_score']
            n_core[(la, lb)] = n_core[(lb, la)] = result['n_core']

    return {'zscores': zscores, 'raw_scores': raw_scores, 'n_core': n_core, 'labels': labels}


def _zscores_to_distance_matrix(zscores: dict, transform: str, exp_scale: float):
    """Build (labels, matrix) from Z-scores. Z higher = more similar -> distance lower."""
    import math

    labels = sorted(set(a for ab in zscores for a in ab))
    n = len(labels)
    if n == 0:
        return [], []
    observed = [v for (a, b), v in zscores.items() if a < b]
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


def _build_tree_from_zscores(zscores: dict, transform: str, exp_scale: float,
                             root: str, midpoint_root: bool) -> str:
    """Build Newick tree from Z-score dict."""
    labels, matrix = _zscores_to_distance_matrix(zscores, transform, exp_scale)
    if len(labels) < 2:
        return "();"
    return _build_tree_newick(labels, matrix, root, midpoint_root)


def _rank_structures(zscores: dict) -> list:
    """Rank by average Z-score (descending)."""
    labels = sorted(set(a for ab in zscores for a in ab))
    out = []
    for a in labels:
        vals = [zscores.get((a, b)) or zscores.get((b, a)) for b in labels if b != a]
        vals = [float(v) for v in vals if v is not None]
        out.append(
            {
                "label": a,
                "n_pairs": len(vals),
                "avg_z": sum(vals) / len(vals) if vals else 0.0,
                "max_z": max(vals) if vals else 0.0,
            }
        )
    out.sort(key=lambda r: (r["avg_z"], r["max_z"]), reverse=True)
    return out


def _plot_tree(newick: str, out_path: str, title: str) -> bool:
    """Generate tree plot. Returns True on success."""
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


def main():
    ap = argparse.ArgumentParser(
        description='Compute Dali-like structural similarity score. Pairwise or all-vs-all with tree output.',
        epilog='''
Pairwise: dali_score.py pdb_a pdb_b [OPTIONS]
All-vs-all: dali_score.py --all-vs-all dir_or_file [dir_or_file ...] [OPTIONS]

Alignment sources: DaliLite -> alignment file -> biotite -> sequence-order.
  DaliLite: --dalilite-path DIR or DALILITE_HOME (import.pl needs mkdssp; see bin/mpidali.pm $DSSP_EXE)

Tree/dendrogram (all-vs-all only): --output-tree, --output-plot, --output-matrix, --output-ranking.
''')
    ap.add_argument('paths', nargs='*', metavar='PATH',
                    help='Pairwise: pdb_a pdb_b. All-vs-all: dir(s) or file(s)')
    ap.add_argument('--all-vs-all', action='store_true',
                    help='Compare all structures (from dirs/files) pairwise; enable tree outputs')
    ap.add_argument('--filter', metavar='PATTERN',
                    help='Filter files by name (substring or glob) in all-vs-all mode')
    ap.add_argument('-a', '--alignment', metavar='FILE',
                    help='Alignment file with residue equivalences (TSV/CSV) [pairwise only]')
    ap.add_argument('--chain-a', metavar='ID', help='Chain ID for structure A (optional)')
    ap.add_argument('--chain-b', metavar='ID', help='Chain ID for structure B (optional)')
    ap.add_argument('--no-dalilite', action='store_true',
                    help='Skip DaliLite (use biotite/alignment/sequence-order)')
    ap.add_argument('--dalilite-path', metavar='DIR',
                    help='DaliLite installation directory (or set DALILITE_HOME); mkdssp must match bin/mpidali.pm')
    ap.add_argument('--no-biotite', action='store_true',
                    help='Skip biotite structural alignment')
    ap.add_argument('--no-sequence-order', action='store_true',
                    help='Do not fall back to sequence-order matching')
    ap.add_argument('-o', '--output', metavar='FILE',
                    help='Pairwise: write pair result CSV. All-vs-all: write pairwise table CSV')
    ap.add_argument('--output-tree', metavar='FILE', default='dali_tree.nwk',
                    help='All-vs-all: write Newick tree (default: dali_tree.nwk)')
    ap.add_argument('--output-plot', metavar='FILE',
                    help='All-vs-all: write dendrogram plot (requires ete3 or biopython+matplotlib)')
    ap.add_argument('--output-matrix', metavar='FILE',
                    help='All-vs-all: write Z-score matrix CSV')
    ap.add_argument('--output-ranking', metavar='FILE', default='dali_ranking.csv',
                    help='All-vs-all: write ranking CSV (default: dali_ranking.csv)')
    ap.add_argument('--transform', choices=['inv', 'maxminus', 'exp'], default='inv',
                    help='Z-score to distance transform for tree (default: inv)')
    ap.add_argument('--exp-scale', type=float, default=10.0,
                    help='Scale for exp transform: d=exp(-Z/scale)')
    ap.add_argument('--root', metavar='LABEL', help='Outgroup label to root tree on')
    ap.add_argument('--no-midpoint-root', action='store_true',
                    help='Disable midpoint rooting (when no --root specified)')
    ap.add_argument('-q', '--quiet', action='store_true', help='Reduce progress output (all-vs-all)')

    args = ap.parse_args()

    all_vs_all = args.all_vs_all
    if all_vs_all:
        if len(args.paths) < 2:
            ap.error('--all-vs-all requires at least 2 paths (dirs or PDB files)')
        files = _collect_pdb_files(args.paths, args.filter)
        if len(files) < 2:
            ap.error(f'Found {len(files)} structure file(s); need at least 2 for all-vs-all')
        print(f"All-vs-all: {len(files)} structures, {len(files)*(len(files)-1)//2} pairs", file=sys.stderr)

        data = run_all_vs_all(
            files,
            chain_a=args.chain_a, chain_b=args.chain_b,
            use_dalilite=not args.no_dalilite,
            use_biotite=not args.no_biotite,
            use_sequence_order=not args.no_sequence_order,
            dalilite_path=args.dalilite_path,
            verbose=not args.quiet,
        )
        zscores = data['zscores']
        labels = data['labels']

        if not zscores:
            print("Error: No pairwise scores computed.", file=sys.stderr)
            sys.exit(1)

        # Pairwise table CSV
        if args.output:
            os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
            with open(args.output, 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['pdb_a', 'pdb_b', 'raw_score', 'z_score', 'n_core'])
                for (la, lb), z in sorted(zscores.items(), key=lambda x: (x[0][0], x[0][1])):
                    if la < lb:
                        r = data['raw_scores'].get((la, lb), '')
                        nc = data['n_core'].get((la, lb), '')
                        w.writerow([la, lb, r, z, nc])
            print(f"Wrote pairwise table: {args.output}")

        # Ranking CSV
        ranking = _rank_structures(zscores)
        os.makedirs(os.path.dirname(args.output_ranking) or '.', exist_ok=True)
        with open(args.output_ranking, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['label', 'n_pairs', 'avg_z', 'max_z'])
            w.writeheader()
            w.writerows(ranking)
        print(f"Wrote ranking: {args.output_ranking}")

        # Z-score matrix CSV
        if args.output_matrix:
            os.makedirs(os.path.dirname(args.output_matrix) or '.', exist_ok=True)
            with open(args.output_matrix, 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Model'] + labels)
                for la in labels:
                    row = [la]
                    for lb in labels:
                        z = zscores.get((la, lb), zscores.get((lb, la)))
                        row.append(f"{z:.4f}" if z is not None else '-')
                    w.writerow(row)
            print(f"Wrote Z-score matrix: {args.output_matrix}")

        # Newick tree
        newick = _build_tree_from_zscores(
            zscores, args.transform, args.exp_scale,
            args.root, not args.no_midpoint_root,
        )
        os.makedirs(os.path.dirname(args.output_tree) or '.', exist_ok=True)
        with open(args.output_tree, 'w') as f:
            f.write(newick)
        print(f"Wrote tree: {args.output_tree}")

        # Plot
        if args.output_plot:
            if _plot_tree(newick, args.output_plot, "Dali structural similarity tree"):
                print(f"Wrote plot: {args.output_plot}")
            else:
                print("Plot failed: install ete3 or biopython+matplotlib", file=sys.stderr)

        sys.exit(0)

    # Pairwise mode
    if len(args.paths) != 2:
        ap.error('Pairwise mode requires exactly 2 paths (pdb_a pdb_b)')
    pdb_a, pdb_b = args.paths[0], args.paths[1]

    result = run(
        pdb_a, pdb_b,
        alignment_file=args.alignment,
        chain_a=args.chain_a, chain_b=args.chain_b,
        use_dalilite=not args.no_dalilite,
        use_biotite=not args.no_biotite,
        use_sequence_order=not args.no_sequence_order,
        dalilite_path=args.dalilite_path,
    )

    if 'error' in result:
        print(f"Error: {result['error']}", file=sys.stderr)
        if 'hint' in result:
            print(f"Hint: {result['hint']}", file=sys.stderr)
        sys.exit(1)

    print(f"raw_score: {result['raw_score']:.2f}")
    print(f"z_score:   {result['z_score']:.2f}")
    if result.get('dalilite_rmsd') is not None:
        print(f"dalilite_rmsd: {result['dalilite_rmsd']:.2f} Å")
    print(f"n_core:    {result['n_core']}")
    print(f"L_A, L_B:  {result['L_A']}, {result['L_B']}")
    print(f"source:    {result['alignment_source']}")

    if args.output:
        with open(args.output, 'w', newline='') as f:
            w = csv.writer(f)
            row_names = ['pdb_a', 'pdb_b', 'raw_score', 'z_score', 'n_core', 'alignment_source']
            row_vals = [
                pdb_a, pdb_b,
                f"{result['raw_score']:.4f}", f"{result['z_score']:.4f}",
                result['n_core'], result['alignment_source'],
            ]
            if result.get('dalilite_rmsd') is not None:
                row_names.extend(['dalilite_rmsd'])
                row_vals.append(f"{result['dalilite_rmsd']:.4f}")
            w.writerow(row_names)
            w.writerow(row_vals)
        print(f"Wrote {args.output}")


if __name__ == '__main__':
    main()
