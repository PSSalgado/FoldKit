# Step-by-step guide: SSM superimposition for subfolder structure

Superimpose CIF/PDB models in nested subfolders to a reference using SSM in Coot, record the log, and extract RMSD to a separate file using the scripts in this folder.

**Example structure:**
```
<project>/                          e.g. condition_1/, condition_2/, ...
  ├── set_a/
  │     └── prefix_set_a_model_*.cif
  └── set_b/
        └── (same naming pattern with set_b)
```

---

## Prerequisites

- **Coot** installed and on `PATH` (scripts call `coot --script ...`).
- **Reference structure**: one PDB or CIF file to which all models will be aligned.
- **Root folder**: the directory that contains structure subfolders and CIF/PDB files.

---

## Step 1: Set paths (run from the repository root)

Command lines below assume the **FoldKit repository root** as the working directory (so `python superimposition/...` resolves). Alternatively, `cd superimposition` and run `python superimpose_coot_SSM.py` (no directory prefix in the command).

Set the reference and root folder (adjust paths to local locations):

```bash
REFERENCE="/path/to/project/reference.pdb"   # or .cif
ROOT_FOLDER="/path/to/folder/containing/models"
```

---

## Step 2: Run SSM superposition (collects all CIF/PDB in subfolders)

The script **recursively** finds all `.pdb` and `.cif` under `ROOT_FOLDER` (including all subfolders), then superimposes each to the reference using SSM in Coot.

**Option A – Default (first chain of each structure):**

```bash
python superimposition/superimpose_coot_SSM.py "$REFERENCE" "$ROOT_FOLDER"
```

**Option B – Only files matching the naming scheme (recommended):**

To restrict to files whose names contain a substring (e.g. `pattern`):

```bash
python superimposition/superimpose_coot_SSM.py --filter=set_a "$REFERENCE" "$ROOT_FOLDER"
```

Or only models that have `_model_` in the name:

```bash
python superimposition/superimpose_coot_SSM.py --filter=model_ "$REFERENCE" "$ROOT_FOLDER"
```

**Note:** `--filter` is a **substring** match (no glob). Use e.g. `--filter=tag` to include only matching models; don’t use `*` in the filter.

**Option C – Specific chains (reference and/or model):**

Use **`superimpose_coot_SSM.py`** with **`--ref-chain`** and/or **`--model-chain`** (the other defaults to `A` if only one is set):

```bash
python superimposition/superimpose_coot_SSM.py --ref-chain=B --model-chain=B --filter=set_a "$REFERENCE" "$ROOT_FOLDER"
```

- A confirmation prompt appears when `SSMaligned2_<ref_name>/` already exists (files may be overwritten).
- Coot will open, run the script, and the log is written automatically.

**Result of Step 2:**

- Output directory: `SSMaligned2_<ref_name>/` (e.g. `SSMaligned2_ref/`), created in the **current working directory**.
- Aligned models: `SSMaligned2_<ref_name>/<model>_SSMaligned2_<ref_name>.pdb`
- **Coot log:** `SSMaligned2_<ref_name>/coot_log.txt`

---

## Step 3: Log is already written

No extra step. The log is created in Step 2 at:

```text
SSMaligned2_<ref_name>/coot_log.txt
```

It contains Coot stdout/stderr, including SSM messages (e.g. “Aligning … to …”, “INFO: core rmsd”, “number of residues”, etc.).

---

## Step 4: Extract RMSD and write to a separate file

Use **`extract_rmsd.py`**: one script for both SSM and LSQ logs. Pass **`--format ssm`** for SSM (or **`--format auto`**, the default, to detect from the log).

### For SSM logs (this workflow)

Parses **SSM** Coot logs: recognises “Superposing &lt;model&gt; (chain X) onto reference (chain Y)” and, if present, Coot’s “Aligning … to …”, then writes that alignment line plus the RMSD block (core rmsd, number of residues, etc.) for each pair. Output: **`SSMaligned2_<ref_name>/rmsd_SSM_values.txt`**.

```bash
python ranking/extract_rmsd.py --format ssm "SSMaligned2_<ref_name>/coot_log.txt"
# or: python ranking/extract_rmsd.py --format auto "SSMaligned2_<ref_name>/coot_log.txt"
```

### For LSQ logs

Use **`--format lsq`** (or auto on an LSQ-only log). Optional: `--aligned=PATTERN`, `--reference=PATTERN`, `--debug`.

```bash
python ranking/extract_rmsd.py --format lsq "LSQaligned2_<ref_name>/coot_log.txt"
```

---

## Summary checklist

| Step | Action | Result |
|------|--------|--------|
| 1 | From repo root, set `REFERENCE` and `ROOT_FOLDER` | — |
| 2 | `python superimposition/superimpose_coot_SSM.py [--filter=...] [--ref-chain=...] [--model-chain=...] "$REFERENCE" "$ROOT_FOLDER"` | `SSMaligned2_<ref>/` with aligned PDBs + `coot_log.txt` |
| 3 | (automatic) | Log in `SSMaligned2_<ref>/coot_log.txt` |
| 4 | `python ranking/extract_rmsd.py --format ssm SSMaligned2_<ref>/coot_log.txt` | `SSMaligned2_<ref>/rmsd_SSM_values.txt` (alignment + RMSD per pair) |

---

## One-liner (after Step 2)

For **SSM** logs (this workflow):

```bash
python ranking/extract_rmsd.py --format ssm SSMaligned2_REFNAME/coot_log.txt
```

Replace `REFNAME` with the reference base name (no extension).

---

## Optional: Custom RMSD output path

To write RMSD output to a different file or directory:

```bash
python ranking/extract_rmsd.py file SSMaligned2_<ref_name>/coot_log.txt -o /path/to/rmsd_output.txt
```

For SSM and LSQ logs use **`extract_rmsd.py`** (single-log mode, optional **`file`** keyword, or **`--dir=...`** to scan logs under a tree) to get alignment and RMSD information into separate text files.

---

## Notes

- **Working directory:** Aligned outputs and the log are written under `SSMaligned2_<ref_name>/` in the directory from which the script is executed.
- **Reference path:** Use an absolute path for `REFERENCE` so Coot can resolve it regardless of the current working directory.
- **Many models:** Coot may take longer with large CIF sets; the log is updated during execution.
- **Chain choice:** Pass **`--ref-chain`** and **`--model-chain`** to **`superimpose_coot_SSM.py`** for explicit chains; omit both to use the first chain in each structure (default mode).
- **Coot window:** One-to-many and single-set all-vs-all keep Coot open by default; use `--not-interactive` to exit when finished. AxB (two sets with `--ref-filter` / `--model-filter`) is non-interactive by default (full matrix without pausing); use `--interactive` to keep Coot open after reload.
- **Which mode:** For **SSM** runs use **`extract_rmsd.py --format ssm`** (or `--format auto`). For **LSQ** use **`--format lsq`** with optional `--aligned` / `--reference` filters.
