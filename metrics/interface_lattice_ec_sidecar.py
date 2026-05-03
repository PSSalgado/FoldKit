#!/usr/bin/env python3
"""
JSON sidecar for two-phase lattice EC: SASA/BSA phase then EC-only phase.

See README (lattice EC, two-phase runs) for usage.
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from typing import Any, Mapping, Sequence

SCHEMA_VERSION = 1


def sha256_file(path: str, chunk_size: int = 1 << 20) -> str:
    """SHA-256 of file bytes."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _json_safe(x: Any) -> Any:
    """Recursively convert to JSON-serializable types (numpy-safe)."""
    if x is None or isinstance(x, (str, bool)):
        return x
    if isinstance(x, int) and not isinstance(x, bool):
        return int(x)
    if isinstance(x, float):
        return float(x)
    try:
        import numpy as np  # type: ignore

        if isinstance(x, (np.integer,)):
            return int(x)
        if isinstance(x, (np.floating,)):
            return float(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
    except Exception:
        pass
    if isinstance(x, dict):
        return {str(k): _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_safe(v) for v in x]
    if isinstance(x, set):
        return sorted(_json_safe(v) for v in x)
    return str(x)


def build_cli_signature(
    *,
    reference_chain_id: str,
    focus_chains: Sequence[str] | None,
    contact_distance: float,
    skip_accessibility_sasa: bool,
    ec_max_contact_points: int | None,
) -> dict[str, Any]:
    fc = sorted([str(c).strip() for c in (focus_chains or []) if str(c).strip()])
    return {
        "reference_chain_id": str(reference_chain_id).strip(),
        "focus_chains": fc,
        "contact_distance": float(contact_distance),
        "skip_accessibility_sasa": bool(skip_accessibility_sasa),
        "ec_max_contact_points": int(ec_max_contact_points)
        if ec_max_contact_points is not None
        else None,
    }


def validate_cli_signature(sidecar_sig: Mapping[str, Any], current: Mapping[str, Any]) -> None:
    keys = (
        "reference_chain_id",
        "focus_chains",
        "contact_distance",
        "skip_accessibility_sasa",
        "ec_max_contact_points",
    )
    for k in keys:
        a, b = sidecar_sig.get(k), current.get(k)
        if a != b:
            raise ValueError(
                f"Sidecar cli_signature mismatch for {k!r}: sidecar={a!r} current={b!r}. "
                "Use the same flags as the phase-sasa run."
            )


def validate_input_fingerprint(sidecar: Mapping[str, Any], pdb_path: str) -> None:
    inp = sidecar.get("input") or {}
    expected = inp.get("sha256")
    if not expected:
        raise ValueError("Sidecar missing input.sha256")
    got = sha256_file(pdb_path)
    if got.lower() != str(expected).lower():
        raise ValueError(
            "Input structure bytes do not match sidecar fingerprint (SHA-256). "
            f"Expected {expected}, got {got} for {pdb_path!r}. Use the identical PDB file."
        )
    recorded = inp.get("path_resolved")
    if recorded and os.path.abspath(pdb_path) != str(recorded):
        # Non-fatal warning could be logged; strict check is hash only
        pass


def write_sidecar(
    path: str,
    *,
    pdb_path: str,
    results: Mapping[str, Any],
    cli_signature: Mapping[str, Any],
) -> None:
    """Write sidecar JSON after phase sasa."""
    pdb_abs = os.path.abspath(pdb_path)
    payload = {
        "foldkit_sidecar_version": SCHEMA_VERSION,
        "schema_version": SCHEMA_VERSION,
        "generator": "interface_analyser_lattice_ec.py",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input": {
            "path": pdb_path,
            "path_resolved": pdb_abs,
            "sha256": sha256_file(pdb_abs),
            "size_bytes": os.path.getsize(pdb_abs),
        },
        "cli_signature": _json_safe(dict(cli_signature)),
        "summary": _json_safe(dict(results.get("summary") or {})),
        "interfaces": _json_safe(list(results.get("interfaces") or [])),
    }
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def read_sidecar(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("Sidecar root must be a JSON object")
    ver = data.get("foldkit_sidecar_version") or data.get("schema_version")
    if int(ver or 0) != SCHEMA_VERSION:
        raise ValueError(f"Unsupported sidecar schema version {ver!r} (expected {SCHEMA_VERSION})")
    return data


def results_from_sidecar(sidecar: Mapping[str, Any]) -> dict[str, Any]:
    """Rebuild minimal results dict for apply_ec_phase (copies so phase ec can mutate safely)."""
    raw_if = sidecar.get("interfaces") or []
    return {
        "interfaces": [dict(x) for x in raw_if],
        "summary": dict(sidecar.get("summary") or {}),
    }


def validate_interfaces_for_ec(interfaces: Sequence[Mapping[str, Any]]) -> None:
    for i, iface in enumerate(interfaces):
        for key in ("chain1_id", "chain2_id"):
            if iface.get(key) is None:
                raise ValueError(f"interfaces[{i}] missing {key!r}")
        if "buried_surface_area" not in iface:
            raise ValueError(f"interfaces[{i}] missing buried_surface_area")
