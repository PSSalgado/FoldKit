"""Safe basenames for aligned PDB outputs (Linux NAME_MAX is 255 bytes)."""

from __future__ import annotations

import hashlib

MAX_ALIGNED_PDB_BASENAME = 255

# Fallback when input stems are still too long for pairwise output names.
COOT_ALIGNED_PDB_HELPERS_PY = """
import hashlib

def foldkit_aligned_pdb_basename(model_name, ref_name, tag):
    full = model_name + tag + ref_name + ".pdb"
    if len(full) <= 255:
        return full
    mh = hashlib.sha1(model_name.encode("utf-8")).hexdigest()[:10]
    rh = hashlib.sha1(ref_name.encode("utf-8")).hexdigest()[:10]
    overhead = len(tag) + len(".pdb") + len(mh) + len(rh) + 2
    budget = 255 - overhead
    if budget < 20:
        return mh + tag + rh + ".pdb"
    mlen = min(len(model_name), budget // 2)
    rlen = min(len(ref_name), budget - mlen)
    return model_name[:mlen] + "_" + mh + tag + ref_name[:rlen] + "_" + rh + ".pdb"
"""


def aligned_pdb_basename(model_stem: str, ref_stem: str, tag: str) -> str:
    """Return a filesystem-safe aligned PDB basename (<= 255 chars)."""
    full = f"{model_stem}{tag}{ref_stem}.pdb"
    if len(full) <= MAX_ALIGNED_PDB_BASENAME:
        return full
    mh = hashlib.sha1(model_stem.encode("utf-8")).hexdigest()[:10]
    rh = hashlib.sha1(ref_stem.encode("utf-8")).hexdigest()[:10]
    overhead = len(tag) + len(".pdb") + len(mh) + len(rh) + 2
    budget = MAX_ALIGNED_PDB_BASENAME - overhead
    if budget < 20:
        return f"{mh}{tag}{rh}.pdb"
    mlen = min(len(model_stem), budget // 2)
    rlen = min(len(ref_stem), budget - mlen)
    return f"{model_stem[:mlen]}_{mh}{tag}{ref_stem[:rlen]}_{rh}.pdb"
