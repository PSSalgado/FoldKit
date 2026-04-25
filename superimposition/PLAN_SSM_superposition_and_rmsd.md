# Step-by-step plan: SSM superimposition and RMSD extraction

Use the existing scripts in this folder to:
1. Superimpose all PDBs (in a folder and subfolders) to a reference PDB using SSM, with optional chain selection.
2. Record the Coot log.
3. Extract RMSD text from the log and write it to a separate file.

---

## Step 1: Superimpose all PDBs to a reference using SSM

Scripts already collect **all `.pdb` and `.cif`** files from the given directory(ies) **and their subdirectories** (recursive), then run Coot with SSM superposition. The Coot stdout/stderr is captured into a log file.

**Coot window:** For **one-to-many** and **single-set all-vs-all**, Coot **stays open** by default for inspection; pass **`--not-interactive`** to exit when finished. **AxB** (two sets via **`--ref-filter`** / **`--model-filter`**) is **non-interactive by default** so each reference-to-model pair runs in one session; pass **`--interactive`** to keep Coot open after reload.

### Option A: Use first chain of each structure (no chain option)

Use **`superimpose_coot_SSM.py`**:

```bash
# From the FoldKit repository root:
python superimposition/superimpose_coot_SSM.py reference.pdb /path/to/models
```

- **reference.pdb**: path to the reference structure (PDB or CIF).
- **/path/to/models**: one or more directories; all `.pdb`/`.cif` under them (including subfolders) are aligned to the reference.

Optional:

- **`--filter=PATTERN`**: only process files whose **filename** contains `PATTERN` (e.g. `--filter=set_a`).

Examples:

```bash
# All PDBs/CIFs under models/ and subfolders
python superimposition/superimpose_coot_SSM.py reference.pdb models/

# Multiple root folders
python superimposition/superimpose_coot_SSM.py reference.pdb dir1 dir2

# Only files matching a substring in the name
python superimposition/superimpose_coot_SSM.py reference.pdb models/ --filter=set_a
```

**Result:**

- Output directory: `SSMaligned2_<ref_name>/` (e.g. `SSMaligned2_ref/`).
- Aligned models: `SSMaligned2_<ref_name>/<model>_SSMaligned2_<ref_name>.pdb`.
- **Coot log**: `SSMaligned2_<ref_name>/coot_log.txt` (recorded automatically).

---

### Option B: Use a specific chain (reference and/or model)

Use **`superimpose_coot_SSM.py`** with the same flags:

```bash
python superimposition/superimpose_coot_SSM.py [--ref-chain=CHAIN] [--model-chain=CHAIN] [--filter=PATTERN] reference.pdb /path/to/models
```

- **`--ref-chain=CHAIN`**: chain ID in the **reference** (default: `A`).
- **`--model-chain=CHAIN`**: chain ID in **each model** (default: `A`).
- **`--filter=PATTERN`**: same as in Option A.

Examples:

```bash
# Reference chain A, model chain A (default)
python superimposition/superimpose_coot_SSM.py reference.pdb models/

# Reference chain B, all models use chain A
python superimposition/superimpose_coot_SSM.py --ref-chain=B reference.pdb models/

# Both chains specified
python superimposition/superimpose_coot_SSM.py --ref-chain=A --model-chain=A reference.pdb models/
```

**Result:** Same as Option A: `SSMaligned2_<ref_name>/` with aligned PDBs and **`coot_log.txt`**.

---

## Step 2: Coot log is already recorded

The script runs Coot with:

- `stdout` and `stderr` piped and appended to a log file.
- Log path: **`<output_dir>/coot_log.txt`** (i.e. `SSMaligned2_<ref_name>/coot_log.txt`).

No extra step is needed to “record” the log; it is created during Step 1.

---

## Step 3: Extract RMSD from the log and write to a separate file

Use one of the extractors on the Coot log produced in Step 1. They look for SSM-style lines (e.g. `INFO: core rmsd`, “number of residues”, etc.) and write a separate RMSD file.

### Option 1: `extract_rmsd.py` (SSM or LSQ)

For **SSM** logs, use **`--format ssm`** or **`--format auto`** (default). Writes **`rmsd_SSM_values.txt`** next to the log:

```bash
python ranking/extract_rmsd.py --format ssm SSMaligned2_<ref_name>/coot_log.txt
python ranking/extract_rmsd.py --format auto SSMaligned2_ref/coot_log.txt
```

For **LSQ** logs, use **`--format lsq`** with optional **`--aligned`**, **`--reference`**, **`--debug`**:

```bash
python ranking/extract_rmsd.py --format lsq LSQaligned2_ref/coot_log.txt --aligned=set_a --debug
```

---

### Option 2: `extract_rmsd.py` (single file or all logs under a dir)

- **Single log file** (output path optional):

  ```bash
  python ranking/extract_rmsd.py file SSMaligned2_ref/coot_log.txt
  python ranking/extract_rmsd.py file SSMaligned2_ref/coot_log.txt -o /path/to/rmsd_output.txt
  ```

- **All `coot_log.txt` under a directory** (e.g. after several runs):

  ```bash
  python ranking/extract_rmsd.py --dir=/path/to/base_dir --format ssm
  python ranking/extract_rmsd.py --dir /path/to/base_dir -o /path/to/rmsd_output_dir
  ```

Use **`extract_rmsd.py`** in default (single-log) mode, with the optional **`file`** keyword, or with **`--dir=DIR`**, for batch extraction and optional custom output paths.

---

## Summary: minimal workflow

1. **Superimpose (all PDBs in folder + subfolders, SSM, optional chain)**  
   - Default (first chain in each structure):  
     `python superimposition/superimpose_coot_SSM.py reference.pdb /path/to/models`  
   - Specific chains:  
     `python superimposition/superimpose_coot_SSM.py --ref-chain=A --model-chain=A reference.pdb /path/to/models`

2. **Log**  
   Already in `SSMaligned2_<ref_name>/coot_log.txt`.

3. **Extract RMSD to a separate file**  
   `python ranking/extract_rmsd.py --format ssm SSMaligned2_<ref_name>/coot_log.txt`  
   → RMSD text in `SSMaligned2_<ref_name>/rmsd_SSM_values.txt`.

---

## One-liner (after superimposition)

```bash
python ranking/extract_rmsd.py --format ssm SSMaligned2_REFNAME/coot_log.txt
```

Replace `REFNAME` with the reference filename (no extension).
