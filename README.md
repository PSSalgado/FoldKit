# FoldKit — protein structure analysis scripts

**FoldKit** is a collection of Python (and supporting R) tools for working with macromolecular **3D structures** in PDB or mmCIF format. The scripts cover **file preparation**, **pairwise and batch superimposition**, **similarity scoring and structure-based phylogenies** (with tabular and graphical outputs), and **quantitative crystal-packing and lattice metrics**. Some workflows use **Coot** for interactive superposition and log-based RMSD extraction; others use **DaliLite**, **TM-align**, or pure Python/R and do not require Coot.

Scripts and docs are grouped in top-level folders that mirror this README:

- `file_management/`
- `superimposition/`
- `ranking/`
- `metrics/`
- `appendix/`

From the repository root run `python <category>/script_name.py` (or `Rscript ranking/...` where noted).

**Example naming (this file, `metrics/*.md`, and script `--help` / docstrings):** `model_01`, `model_02`, … denote structure files or row/column labels; `set_a`, `set_b`, … denote filter substrings, `--sets` groups, or `--tags` batch tokens; `condition_1`, `condition_2`, … denote condition subfolders in batch drivers; `ref_id` / `model_id` are filename tokens in LSQ `--pattern` mode; `ref_01` illustrates a reference name in Coot log filters; `results.txt` denotes a merged interface text report (e.g. from `interface_analyser_asu_charge.py -o`); `contact_results.txt` denotes a merged report from `contact_analyser.py -o`; `./out` denotes a generic output directory; chain IDs in examples (`A`, `B`, …) are placeholders. Output name patterns such as `rmsd_table_<suffix>.csv` or `rmsd_heatmap_<suffix>.<ext>` (for example `.svg` from `rmsd_to_csv.py` batch mode) match the scripts that create them.

**Paths in examples:** `/path/to/project/...` marks a generic working tree (substitute your own directories or files; absolute paths are safest). Third-party installs may use placeholders such as `/path/to/DaliLite`. Command lines assume you run `python <category>/<script>.py` from the **repository root** unless stated otherwise.

## Contents


| Section                                                                                                | What it covers                                                                                      |
| ------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------- |
| [File management](#file-management)                                                                    | Renaming files; PDB chain ID replacement and merge; residue-range trimming (no Coot)                |
| [Superimposition](#superimposition)                                                                    | Coot SSM/LSQ and trim workflows; batch LSQ; DaliLite superposed coordinates; opening models in Coot |
| [Ranking, scoring, phylogeny, and graphical outputs](#ranking-scoring-phylogeny-and-graphical-outputs) | RMSD from logs; RMSD tables and heatmaps; Dali Z-scores; neighbour-joining trees and plots           |
| [Metrics (crystal packing and lattice)](#metrics-crystal-packing-and-lattice)                          | Packing density, interfaces, contacts; multi-copy / lattice workflows; optional batch JSON         |
| [Appendix](#appendix)                                                                                  | General notes, output layout, troubleshooting                                                       |


## Prerequisites

**Common:** Python 3; structure files in PDB or mmCIF where applicable.

**Depending on what you run:**

- **Coot** (command-line): for `superimpose_coot_*.py`, `trim_superimposeLSQ.py` (superposition modes), `open_models_in_coot.py`, and log parsing in `extract_rmsd.py`. Not needed for `trim_models.py` (trim-only).
- **R** (optional): `ranking/create_rmsd_heatmap.R` and other R utilities where documented; not required for the default `crystal_packing_analyser` CLI.
- **DaliLite** + `mkdssp`: `dalilite_superpose_scores.py` and DaliLite-backed modes of `dali_score.py` (see [Appendix — DaliLite](#dalilite-dali_scorepy-dalilite_superpose_scorespy)).
- **TM-align**: optional distance input for `structure_phylogeny.py --from-pdb`.
- `gemmi` (`pip install gemmi`): helpful for some CIF handling in superposition workflows.

---

## File management

### `rename_files.py`

**Legacy (two arguments):** in a single directory, renames PDBs matching `{prefix}_YYYY_MM_DD_HH_MM_*.pdb` to `{prefix}_*.pdb` (non-recursive).

```bash
python file_management/rename_files.py DIRECTORY PREFIX
```

**Regex mode:** remove or replace patterns in basenames; by default processes subdirectories and may rename directories unless `--files-only` / `--no-recursive`.

```bash
python file_management/rename_files.py DIRECTORY --remove='REGEX'
python file_management/rename_files.py DIRECTORY --replace='REGEX' --with='STRING'
```

Options (regex mode): `--remove`, `--replace`, `--with` (required with `--replace`), `--no-recursive`, `--files-only`

Examples:

```bash
python file_management/rename_files.py /path/to/project/workdir sample_prefix
python file_management/rename_files.py /path/to/project/workdir --remove='\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_'
python file_management/rename_files.py /path/to/project/workdir --replace='fold_' --with='protein_'
python file_management/rename_files.py /path/to/project/workdir --remove='temp_' --no-recursive
```

### `pdb_rechain.py`

Rewrite chain IDs on coordinate records (ATOM/HETATM/ANISOU/TER). Optionally merge one chain into another and renumber residues so the source chain continues after the last residue on the target chain (for example, merging B into A to form a single continuous chain).

- **Key options**
  - **`--merge-map 'FROM:TO,...'`**: Apply multiple merges left-to-right. Each merge moves all records from `FROM` to `TO` and renumbers the former `FROM` residues to continue after the last residue number on `TO` (PDB `resseq` field, max 9999).
  - **`--reorder-chains`**: Reorder coordinate records so each chain is contiguous in the output. Recommended after `--merge-map` to prevent interleaving due to the original record order.
  - **`--rename-sequential`**: Relabel all chains to `A..Z` after merges/renumbering.
  - **`--chain-order 'ID,ID,...'`**: Used with `--rename-sequential`. Specifies the preferred order of existing chain IDs (post-merge) that should become new `A,B,C,...` first; any remaining chains are appended in first-appearance order.
  - **`--renumber-per-chain`**: Renumber residues per chain starting at 1 (in residue first-appearance order). Keeps insertion codes unchanged.
  - **`-o/--output`**: For multiple input files, use either an output directory (one output per input) or a path template containing `{}` (replaced by the input stem), for example `-o '/path/to/project/out/{}_rechain.pdb'`.

```bash
# Merge chain B into chain A; residues from the former chain B continue after the last residue on chain A:
python file_management/pdb_rechain.py /path/to/project/input_dir --pattern '*_models.pdb' -f B -t A \
  --merge-renumber -o /path/to/project/output_dir/

# Rename only (target chain must not already exist in the file):
python file_management/pdb_rechain.py model_01.pdb -f X -t Y -o model_01_Y.pdb

# Multi-copy assemblies: merge several chain pairs, reorder so chains are contiguous, relabel to A..Z,
# and renumber residues per chain (starting at 1).
#
# Example merge list (pairs merged left-to-right):
#   B to A, D to C, F to E, ...
#
# Use --chain-order to control how the surviving chains are mapped onto new A,B,C,… after merging.
python file_management/pdb_rechain.py multicopy.pdb \
  --merge-map 'B:A,D:C,F:E,H:G' \
  --reorder-chains \
  --rename-sequential \
  --chain-order 'A,C,E,G' \
  --renumber-per-chain \
  -o multicopy_merged_ordered.pdb

# Batch processing (several files with identical chain conventions):
python file_management/pdb_rechain.py /path/to/project/models/ --pattern '*.pdb' \
  --merge-map 'B:A,D:C,F:E,H:G' \
  --reorder-chains --rename-sequential --chain-order 'A,C,E,G' \
  --renumber-per-chain \
  -o '/path/to/project/out/{}_merged_ordered.pdb'
```

`SEQRES` and `CRYST1` records are not rewritten; update these separately if required by downstream tools.

### `trim_models.py`

**Residue trimming only** (no superposition, no Coot): harmonise PDB/mmCIF models to the **shortest** residue span found among the inputs. **Self-contained** (stdlib + optional `gemmi` for mmCIF→PDB); it does not import other FoldKit scripts. `trim_superimposeLSQ.py` uses the same trimming implementation via `trim_models.py` for `--trim` / `--trim-only`.

```bash
python file_management/trim_models.py [--filter=set_a,set_b,...] directory1 [directory2 ...]
```

`--filter`: optional comma-separated patterns (basename substring or glob). If omitted, every `*.pdb` / `*.cif` in the given directories is included (output under `trimmed_all/`).

Examples:

```bash
python file_management/trim_models.py models/
python file_management/trim_models.py --filter=set_a models/
python superimposition/trim_superimposeLSQ.py --trim-only --filter=set_a models/
```

Output: `trimmed_<pattern>/` per filter, or `trimmed_all/` with no filter; trimmed PDBs. `gemmi` may be needed for mmCIF inputs.

#### References (file management)

- **wwPDB:** [PDBx/mmCIF and PDB archive documentation](https://www.wwpdb.org/documentation) — file layout and atom-record conventions.
- **Bernstein, F.C. et al.** (1977). The Protein Data Bank: a computer-based archival file for macromolecular structures. *J. Mol. Biol.* 112, 535–542. [DOI: 10.1016/0022-2836(77)90297-6](https://doi.org/10.1016/0022-2836(77)90297-6)

---

## Superimposition

These scripts align models to a reference or run all-vs-all jobs. **Coot-based** tools default to **keeping Coot open** after one-to-many or single-set all-vs-all runs (use `**--not-interactive`** to exit instead); **two-set (AxB)** runs default to **batch exit** unless `**--interactive`**. `**dalilite_superpose_scores.py**` writes superposed coordinates using **DaliLite** (separate install).

#### Main flows vs pattern mode (`--pattern`)

`superimpose_coot_SSM.py` and `superimpose_coot_LSQ.py` expose two **separate** command-line shapes. In a **single** run you use **one or the other**; they are **not** combinable.


| Goal                                                                                         | Invocation                                                                                                                                  |
| -------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| One fixed reference, many models                                                             | Leading `reference.pdb` (and dirs), or `--reference=` / `--ref=`                                                                            |
| All-vs-all on one pool of structures                                                         | `--all-vs-all` and one or more directories; optional `--filter=`                                                                            |
| Two disjoint sets (A×B)                                                                      | `--all-vs-all` with `--ref-filter=` and `--model-filter=` and two directory trees                                                           |
| Pair reference and model files by **filenames** (and optional subdir names) across two roots | `**--pattern` must be the first argument** after the script name, then `reference_dir model_dir ref_pattern model_pattern [target_pattern]` |


**Pattern mode rules**

- `**--pattern` must be `argv[1]`** (the first token after the script). If it is not first, the script uses the main flow instead and will not treat the rest as pattern mode.
- **Do not** combine `--pattern` with `--all-vs-all`, `--reference` / `--ref`, `--filter`, `--ref-filter`, `--model-filter`, or a leading reference positional in the same command.
- Pairing is implemented in `**superimposition/superimpose_pattern_match.py`** (`find_ref_model_matches`), shared by SSM and LSQ so neither script depends on the other.

**Examples — pattern mode (correct)**

```bash
python superimposition/superimpose_coot_SSM.py --pattern /path/to/project/refs /path/to/project/models proteinA fold1
python superimposition/superimpose_coot_LSQ.py --pattern /path/to/project/refs /path/to/project/models proteinA fold1
python superimposition/superimpose_coot_SSM.py --pattern --not-interactive /path/to/project/refs /path/to/project/models proteinA fold1
```

**Examples — invalid (do not use)**

```bash
# --pattern not first → not pattern mode; paths/options will be misread
python superimposition/superimpose_coot_SSM.py /path/to/project/refs /path/to/project/models --pattern proteinA fold1

# cannot mix pattern mode with all-vs-all or filters
python superimposition/superimpose_coot_SSM.py --all-vs-all --pattern refs/ models/ id_a id_b
python superimposition/superimpose_coot_SSM.py --pattern --all-vs-all refs/ models/ id_a id_b
```

### `superimpose_coot_SSM.py`

Secondary Structure Matching (SSM) superposition: one-to-many or all-vs-all.

**Default (no chain flags):** uses the **first chain** in each structure and SSM mode `1` (legacy Coot selection).

**Explicit chains:** pass `**--ref-chain`** and/or `**--model-chain**` to use mmdb-style selections (`/1/CHAIN/*`) and SSM with move flag `0`. If only one of the two is set, the other defaults to `**A**`.

```bash
# One-to-many: align all models to a reference (reference is first positional, or use --reference= / --ref=)
python superimposition/superimpose_coot_SSM.py [--filter=set_a] [--output-dir=DIR] reference.pdb dir1 [dir2 ...]
python superimposition/superimpose_coot_SSM.py [--filter=set_a] --reference=reference.pdb dir1 [dir2 ...]

# All-vs-all: each structure used as reference in turn
python superimposition/superimpose_coot_SSM.py [--filter=set_a] [--output-dir=DIR] --all-vs-all dir1 [dir2 ...]

# Explicit chains (one-to-many or all-vs-all)
python superimposition/superimpose_coot_SSM.py --ref-chain=B --model-chain=B reference.pdb models/
python superimposition/superimpose_coot_SSM.py --all-vs-all --ref-chain=B --model-chain=B models/

# All-vs-all between two sets (AxB), potentially from different trees
python superimposition/superimpose_coot_SSM.py --all-vs-all \
  --ref-filter='*set_a*' --model-filter='*set_b*' \
  --ref-chain=A --model-chain=B \
  dir_for_set_a/ dir_for_set_b/

# Same AxB run but keep the Coot window open after all superpositions (default is batch: Coot exits)
python superimposition/superimpose_coot_SSM.py --all-vs-all --interactive \
  --ref-filter='*set_a*' --model-filter='*set_b*' \
  --ref-chain=A --model-chain=B \
  dir_for_set_a/ dir_for_set_b/
```

**AxB mode:** one Coot session runs every reference in A against every structure in B (skipping a file paired with itself). The default is **non-interactive**: the generated Coot script ends with `coot_real_exit(0)  # AxB default: non-interactive batch (exit when done)` so the full matrix of alignments finishes without leaving a window open. Pass `**--interactive`** to reload inputs after the last alignment and keep Coot open (same idea as single-set all-vs-all, which stays open by default).

Examples:

```bash
python superimposition/superimpose_coot_SSM.py reference.pdb dir1 dir2 dir3
python superimposition/superimpose_coot_SSM.py --all-vs-all --filter=set_a models/
python superimposition/superimpose_coot_SSM.py --filter=set_a --output-dir=SSM_run reference.pdb models/
python superimposition/superimpose_coot_SSM.py --ref-chain=C --model-chain=B reference.pdb models/
python superimposition/superimpose_coot_SSM.py --all-vs-all --ref-filter=set_a --model-filter=set_b --ref-chain=A --model-chain=B models_a/ models_b/
```

Options: `--filter` (single-set substring or glob on basename, e.g. `set_a`), `--ref-filter` (reference set A), `--model-filter` (model set B; both require `--all-vs-all`), `--reference` / `--ref` (one-to-many: reference file instead of a leading positional), `--interactive` (AxB only: keep Coot open after reload), `--not-interactive` (one-to-many and single-set all-vs-all: exit Coot when done; default keeps Coot open), `--output-dir` / `--out-dir` (placeholders: `[reference_name]`, `[filter]`), `--ref-chain`, `--model-chain`, `--all-vs-all`

Output: `SSMaligned2_[reference_name]` or `SSMaligned_all_vs_all[*]`; `[model]_SSMaligned2_[reference].pdb`; `coot_log.txt` (AxB also writes `coot_log_AxB_*.txt` in the current working directory); `rmsd_SSM_values.txt` (via `extract_rmsd.py --format ssm` or `--format auto`). Prior outputs whose paths contain `SSMaligned` or `LSQaligned` are skipped when collecting models.

**Pattern mode (`--pattern`):** a **separate** entry point from one-to-many / all-vs-all / AxB (see **Main flows vs pattern mode** above). Uses the same pairing rules as LSQ pattern mode via `**superimpose_pattern_match.find_ref_model_matches`**. Each run writes `coot_log.txt` with the same header fields as LSQ pattern mode (`# Reference:`, `# Models directory:`, `# Number of models:`) plus an SSM-specific mode line; Coot output includes `**Superposing … onto …**` lines before each superposition so `extract_rmsd.py --format ssm` can pair RMSD blocks with models. Output directories are `[subdir]_SSMaligned_[reference_stem]`; aligned filenames use `_SSMaligned_` between model and reference stems.

```bash
# --pattern must be first; optional flags can appear before the directories
python superimposition/superimpose_coot_SSM.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]
```

Pattern-only options (same idea as LSQ): `--strict-position`, `--divider=`, `--ref-file-pattern=`, `--model-file-pattern=`, `--output-suffix=`, `--not-interactive`.

### `superimpose_coot_LSQ.py`

Least Squares (LSQ) superposition: one-to-many or all-vs-all. By default uses the **first chain** in each structure. `**--ref-chain` / `--model-chain`** apply to **one-to-many** runs; **AxB mode** still uses the first chain per structure (explicit-chain LSQ for AxB is not implemented yet).

```bash
# One-to-many: specify reference explicitly
python superimposition/superimpose_coot_LSQ.py [--output-dir=DIR] --reference=reference.pdb dir1 [dir2 ...]

# One-to-many: interactive reference selection from directories
python superimposition/superimpose_coot_LSQ.py [--output-dir=DIR] dir1 [dir2 ...]

# All-vs-all
python superimposition/superimpose_coot_LSQ.py --all-vs-all [--output-dir=DIR] dir1 [dir2 ...]

# All-vs-all between two sets (AxB), potentially from different trees (first chain per structure)
python superimposition/superimpose_coot_LSQ.py --all-vs-all \
  --ref-filter='*set_a*' --model-filter='*set_b*' \
  dir_for_set_a/ dir_for_set_b/

# AxB with Coot left open when finished
python superimposition/superimpose_coot_LSQ.py --all-vs-all --interactive \
  --ref-filter='*set_a*' --model-filter='*set_b*' \
  dir_for_set_a/ dir_for_set_b/
```

**AxB (LSQ):** same behaviour as SSM AxB: default **non-interactive** batch with explicit `coot_real_exit(0)` in the generated script; `**--interactive`** keeps Coot open after reload. Log file `coot_log_AxB_*.txt` in the current working directory.

Examples:

```bash
python superimposition/superimpose_coot_LSQ.py --reference=reference.cif --filter=set_a models/
python superimposition/superimpose_coot_LSQ.py --all-vs-all --filter=set_a models/
python superimposition/superimpose_coot_LSQ.py --all-vs-all --ref-filter=set_a --model-filter=set_b models_a/ models_b/
python superimposition/superimpose_coot_LSQ.py dir1 dir2
```

Options: `--reference` / `--ref`, `--filter` (single-set substring/glob match, e.g. `set_a`), `--ref-filter` (reference set A), `--model-filter` (model set B; both require `--all-vs-all`), `--interactive` (AxB only: keep Coot open after reload), `--not-interactive` (one-to-many and single-set all-vs-all: exit Coot when done; default keeps Coot open), `--ref-chain`, `--model-chain` (one-to-many explicit chains only), `--output-dir` / `--out-dir`

Output: `LSQaligned2_[reference_name]` or `LSQaligned_all_vs_all[*]`; `coot_log.txt` (plus `coot_log_AxB_*.txt` in the working directory for AxB); `rmsd_values.txt` (via `extract_rmsd.py`)

**Pattern mode (`--pattern`):** a **separate** entry point from one-to-many / all-vs-all / AxB (see **Main flows vs pattern mode** above). Pairs reference and model files by **filenames** (and optional subdirectory names) via `**superimpose_pattern_match.find_ref_model_matches`**, then runs the same LSQ Coot job per match. This replaces the old standalone `trim_superimposeLSQ_pattern.py` script.

```bash
# --pattern must be first; optional flags can appear before the directories
python superimposition/superimpose_coot_LSQ.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]
```

Options:

- `--strict-position`: ref_pattern before divider, target_pattern after
- `--divider=STRING`: divider string inside reference filenames
- `--ref-file-pattern=GLOB` (default: `"*.pdb"`)
- `--model-file-pattern=GLOB` (default: `"*.cif"`)
- `--output-suffix=STRING`: output directory infix (default: `"_LSQaligned_"`)
- `--not-interactive`: exit Coot after each alignment set (default: keep Coot open for inspection)

Examples:

```bash
python superimposition/superimpose_coot_LSQ.py --pattern /path/to/project/refs /path/to/project/models ref_id model_id
python superimposition/superimpose_coot_LSQ.py --pattern --ref-file-pattern="*final.pdb" --model-file-pattern="*.pdb" /path/to/project/refs /path/to/project/models ref_id model_id
python superimposition/superimpose_coot_LSQ.py --pattern --divider=LSQaligned2 /path/to/project/refs /path/to/project/models ref_id model_id
```

Output: `[subdir]_LSQaligned_[ref_name]/`; aligned PDBs named `[model]_LSQaligned_[ref].pdb` inside each output directory.

### `trim_superimposeLSQ.py`

Trims models to a **common residue span** (shortest model among the set, optional `--filter` groups), then performs LSQ superposition—**or** trim only without Coot. Trimming is implemented in `**trim_models.py`**; use that script for a standalone trim-only run, or `**--trim-only**` here for the same logic before superposition.

```bash
# Trim and superimpose (all-vs-all or one-to-many with interactive reference)
python superimposition/trim_superimposeLSQ.py --trim [--all-vs-all] [--filter=set_a,set_b] dir1 [dir2 ...]

# Only trim (no Coot) — same implementation as trim_models.py
python superimposition/trim_superimposeLSQ.py --trim-only [--filter=set_a,set_b,...] dir1 [dir2 ...]

# Superimpose without trimming
python superimposition/trim_superimposeLSQ.py [--all-vs-all] [--filter=set_a] dir1 [dir2 ...]
```

Examples:

```bash
python superimposition/trim_superimposeLSQ.py --trim --all-vs-all --filter=set_a,set_b models/
python superimposition/trim_superimposeLSQ.py --all-vs-all models/
python superimposition/trim_superimposeLSQ.py --trim dir1 dir2
python file_management/trim_models.py --filter=set_a dir1
python superimposition/trim_superimposeLSQ.py --trim-only --filter=set_a dir1
```

Options: `--trim`, `--trim-only` (implies `--trim`; incompatible with `--all-vs-all`), `--all-vs-all`, `--filter=set_a[,set_b,...]` (comma-separated substrings)

Output: With trimming, creates `trimmed_set_a/`, `trimmed_set_b/` (or `trimmed_all/`); superposition output in `LSQ_set_a_all_vs_all/`, `LSQ_all_vs_all/`, or `LSQaligned2_[ref]` (one-to-many).

### `run_all_superpositions.py`

Batch driver for **`superimpose_coot_LSQ.py --pattern`**: loops over **conditions** (subfolders under each base, e.g. `condition_1`) and **tags** (same role as `--filter=set_a` in SSM/LSQ). Expects `refs_base/<condition>/<REFERENCE_SUBDIR>/` for reference PDBs (default subdir `LSQaligned2_reference_m0` in the script) and `models_base/<condition>/*_fl_<set>*` for models. The script passes **`--pattern`**, **`--divider=LSQ_`**, and globs derived from **`CONDITION_PREFIX`** (e.g. `run1_sd_set_a` vs `run1_sd_reference`); edit **`CONDITION_PREFIX`**, **`REFERENCE_SUBDIR`**, **`REFERENCE_STEM`**, and **`TAGS`** in the script to match your layout.

```bash
python superimposition/run_all_superpositions.py --tags set_a set_b \
  --ref-base /path/to/project/refs --models-base /path/to/project/models

# Subset of conditions (default: all keys in CONDITION_PREFIX)
python superimposition/run_all_superpositions.py --tags set_a --conditions condition_1 condition_2 \
  --ref-base /path/to/project/refs --models-base /path/to/project/models
```

Options: `--tags` (required unless `TAGS` is set in the script), `--conditions`, `--ref-base`, `--models-base`.

### `dalilite_superpose_scores.py`

Pairwise or **all-vs-all** DaliLite runs with **superposed PDB** output (target on query frame), pairwise or merged CSV, optional Z-matrix CSV, Newick tree, optional dendrogram image, and an optional **Z-score heatmap** (matplotlib; layout and colour bar match **`rmsd_to_csv.py`** via shared **`ranking/foldkit_heatmap.py`**). Requires DaliLite and **mkdssp** (same constraints as `dali_score.py`; see [Appendix — DaliLite](#dalilite-dali_scorepy-dalilite_superpose_scorespy)). In **all-vs-all** mode, structure files are taken from the **root of each given directory only** unless **`--recursive`** is set (then subfolders are included, matching the **default** directory scan in `dali_score.py`; that script can match root-only behaviour with **`--root-only-dirs`**).

**Optional run summary log:** most scripts accept `--log [FILE]` to write a short summary of tasks and errors (no stdout/stderr tee). Logging is optional and summary-only.

**Heatmap (`--heatmap PATH`):** writes a square figure from pairwise Dali Z (colour bar label **Dali Z**). The file extension selects **PNG**, **PDF**, or **SVG**. Options aligned with `rmsd_to_csv.py` / `foldkit_heatmap.py` include `--heatmap-title`, `--cmap`, `--vmin`, `--vmax`, `--short-heatmap-labels`, `--heatmap-diverging-center` (`none` or `median`), `--heatmap-colorbar-orientation` (`vertical` or `horizontal`), and `--heatmap-y-axis-right`. Autoscale and median use **all finite off-diagonal** Z values (not RMSD-style “positive only”). All-vs-all produces the full matrix; pairwise mode produces a 2×2 plot. **Axis order** for the heatmap and **`--output-matrix`** matches **`rmsd_to_csv.py`**: structure labels are **natural-sorted** (alphanumeric, numeric-aware).

**Optional second channel (hatch):** **`--heatmap-n-core-patterns`** encodes **`n_core`** (aligned Cα count per pair) as **binned hatch patterns** on off-diagonal cells; **colour still encodes Dali Z**. Use **`--heatmap-n-core-bins N`** for equal-width bins over observed `n_core` (default **4**), or **`--heatmap-n-core-edges E0,E1,…`** for explicit boundaries (comma-separated; extended to the observed min/max if needed). Legend title uses the script’s label for that quantity (aligned Cα / Dali `n_core`). Hatch line density is controlled in **`foldkit_heatmap.py`** (`_HATCH_REPEAT_*` constants) if defaults look too sparse or too solid for your matrix size.

**Table outputs and Z-score column:** The **`z_score` field in CSV** is the **DaliLite summary Z** when DaliLite returns a value on the hit line. If that is missing, FoldKit fills **`z_score`** with the **empirical (Holm-style) Z** from `dali_score.compute_z_score`. The script states which applies: progress lines use **`Z=… (DaliLite)`** versus **`empirical Z-score=… (FoldKit dali-like)`**; in pairwise mode, **`empirical Z-score:`** or **`Z-score: … (DaliLite reported)`** is printed accordingly. The **pairs** CSV and **`--output-matrix`** are **natural-ordered** by label; the **ranking** CSV (`--output-ranking`, default `dali_ranking.csv`) is ordered **only by Z** (average then maximum), not by label.

**Residue core (optional):** The pairs CSV includes **`core_resnum_min_a` / `core_resnum_max_a` / `core_resnum_min_b` / `core_resnum_max_b`** (PDB residue number span of the aligned core on each side; non-contiguous alignments still show overall min–max). Pass **`--equivalences-dir /path/to/project/equivalences`** to write one TSV per pair (aligned Cα pairs used for **`n_core`** and the **raw** Dali sum in FoldKit, with a short comment header). **`n_core`** is the count of those pairs after mapping DaliLite’s block to your coordinates.

```bash
python ranking/dalilite_superpose_scores.py /path/to/project/models/model_01.pdb /path/to/project/models/model_02.pdb \
  -d /path/to/project/out --dalilite-path /path/to/DaliLite

python ranking/dalilite_superpose_scores.py --all-vs-all /path/to/project/models/ /path/to/project/models/ \
  -d /path/to/project/dali_run1/ --filter '*.pdb' --dalilite-path /path/to/DaliLite \
  -o /path/to/project/dali_run1/pairs.csv --output-matrix /path/to/project/dali_run1/z_matrix.csv \
  --output-ranking /path/to/project/dali_run1/dali_ranking.csv \
  --equivalences-dir /path/to/project/dali_run1/equivalences \
  --heatmap /path/to/project/dali_run1/z_scores.svg --heatmap-diverging-center median --cmap RdYlBu_r \
  --heatmap-n-core-patterns
```

### `open_models_in_coot.py`

Opens PDB/mmCIF in Coot with C-alpha representation. **Directory mode:** pass one or more directories (top-level `*.pdb`, `*.cif`, `*.mmcif`); optional `--filter SUBSTRING` on basename. **Pattern mode:** pass glob patterns, with optional directories to search (recursive).

```bash
python superimposition/open_models_in_coot.py [--filter SUBSTRING] dir1 [dir2 ...]
python superimposition/open_models_in_coot.py pattern [pattern ...] [dir1] [dir2 ...]
```

Examples:

```bash
python superimposition/open_models_in_coot.py models/
python superimposition/open_models_in_coot.py --filter set_a dir1 dir2
python superimposition/open_models_in_coot.py '*set_a*LSQaligned2*ref*pdb'
python superimposition/open_models_in_coot.py '*model_*.cif' models/
```

#### References (superimposition)

- **Emsley, P. et al.** (2010). Features and development of Coot. *Acta Crystallogr. D* **66**, 486–501. [DOI: 10.1107/S0907444910007493](https://doi.org/10.1107/S0907444910007493)
- **Holm, L. & Sander, C.** (1993). Protein structure comparison by alignment of distance matrices. *J. Mol. Biol.* **233**, 123–138. [DOI: 10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489) — Dali / DaliLite distance geometry (see also Holm 2020 below).
- **Holm, L.** (2020). Dali and the persistence of protein shape. *Protein Sci.* **29**, 128–140. [DOI: 10.1002/pro.3749](https://doi.org/10.1002/pro.3749)

---

## Ranking, scoring, phylogeny, and graphical outputs

Scripts here convert superposition logs or similarity tables into **tables, heatmaps, rankings, and trees** (Newick, and PDF, PNG, or SVG figures where supported).

### `extract_rmsd.py`

Extracts RMSD blocks from Coot logs. **LSQ** logs use “Aligning … to …” (optional `--aligned` / `--reference` filters). **SSM** logs use “Superposing … onto …”. Use `--format lsq`, `--format ssm`, or `--format auto` (default) to detect from the log.

```bash
python ranking/extract_rmsd.py /path/to/project/coot_log.txt [--format auto|lsq|ssm] [LSQ options]
```

LSQ-only options: `--aligned`, `-a`; `--reference`, `-r`; `--case-sensitive`, `-c`; `--debug`, `-d`.

Examples:

```bash
python ranking/extract_rmsd.py /path/to/project/LSQaligned2_ref/coot_log.txt --format lsq
python ranking/extract_rmsd.py /path/to/project/SSMaligned2_ref/coot_log.txt --format ssm
python ranking/extract_rmsd.py /path/to/project/coot_log.txt   # auto: picks SSM vs LSQ from log content
python ranking/extract_rmsd.py /path/to/project/coot_log.txt --aligned=set_a --reference=ref_01 --format lsq
python ranking/extract_rmsd.py /path/to/project/coot_log.txt --debug --format lsq
```

Output: `rmsd_values*.txt` (LSQ) or `rmsd_SSM_values*.txt` (SSM).

**Batch modes (same script):**

```bash
# File mode: single log with optional explicit output path
python ranking/extract_rmsd.py file /path/to/project/SSMaligned2_ref/coot_log.txt
python ranking/extract_rmsd.py file /path/to/project/SSMaligned2_ref/coot_log.txt -o /path/to/project/rmsd_output.txt

# Directory scan: sweep all Coot logs under a base directory (--dir=DIR or --dir DIR)
python ranking/extract_rmsd.py --dir=/path/to/project/run1 --format ssm
python ranking/extract_rmsd.py --dir /path/to/project/run1 -o /path/to/project/rmsd_output_dir
```

- **Single-log mode (default)**: same as before; pick SSM vs LSQ automatically or via `--format`, with LSQ filters (`--aligned`, `--reference`, `--case-sensitive`, `--debug`), writing standard `rmsd_SSM_values*.txt` / `rmsd_values*.txt` next to the log.
- **`file` keyword** (optional) and **`--dir`**: same extractors; optional `file` disambiguates a log literally named `file`. **`--dir`** scans recursively; **`-o`** with **`--dir`** sets a mirror root for outputs.

### `create_rmsd_heatmap.R`

Creates heatmap PDFs from RMSD tables.

```bash
Rscript ranking/create_rmsd_heatmap.R /path/to/project/base_directory
Rscript ranking/create_rmsd_heatmap.R /path/to/project/base_directory viridis
Rscript ranking/create_rmsd_heatmap.R /path/to/project/base_directory YlOrRd
```

Usage: `Rscript ranking/create_rmsd_heatmap.R <base_directory> [palette]` (use a directory that contains the RMSD tables this script expects).  
Palette: RdYlBu (default), RdYlGn, YlOrRd, viridis, plasma, etc. (RColorBrewer or viridis options)

Requirements: R with viridis, pheatmap, RColorBrewer

Output: `rmsd_heatmap_[subdomain].pdf`, `combined_rmsd_heatmap.pdf`, `combined_rmsd_heatmap_clustered.pdf`

### `rmsd_to_csv.py`

Converts pairwise RMSD inputs (all-vs-all SSM or LSQ text, or an existing square RMSD CSV) to a square CSV table. Optional **matplotlib** heatmaps require no R. Figures are drawn by **`ranking/foldkit_heatmap.py`** (same module as **`dalilite_superpose_scores.py`** heatmaps; `dalilite` does not import this script).

#### Single-file mode (default)

Pass a single path as `INPUT`; do not use `--scan-dir`. Suitable for one Coot or `extract_rmsd` log, or one RMSD table to reformat.

1. **Read and detect format** — LSQ-style text (`Alignment: Aligning A to B` plus a `core rmsd` line), SSM-style text (`Superposing` / `Aligning` plus `core rmsd`), or a square CSV with a `Model` column and one column per structure (same rules as `structure_phylogeny.py`).
2. **Build one symmetric matrix** — All structure names from pairwise lines are collected; missing pairs remain empty. Row and column order default to a natural sort unless ordering options override it.
3. **Write one CSV** — Square table: column `Model`, then one column per structure, diagonal `-`, off-diagonals RMSD in Å (four decimal places). Default output is `rmsd_table_<suffix>.csv` beside the input (`rmsd_values_run1.txt` → `run1`; generic names fall back to `rmsd`). Use `-o` / `--output` for an explicit path.
4. **Optional heatmap** — `--plot PATH` writes one figure; the extension selects **PNG**, **PDF**, or **SVG**. Use `--title`, `--cmap`, `--vmin`, `--vmax`, and the options in [Heatmap colours](#heatmap-colours).
5. **Optional row/column order** — `--order` (comma-separated list or path to one label per line) or `--order-dir` with `--order-glob` (recursive scan of PDB/mmCIF stems) reorders the matrix before the CSV and plot. Labels must match names in the RMSD file.

For a neighbour-joining tree from one file, run `structure_phylogeny.py` on the same log or on the written CSV.

#### Heatmap colours

Heatmaps share the following controls for `--plot`, for figures under `--heatmap-dir`, and for `--combined-heatmap`. The **hatch** rows apply only to **square** RMSD heatmaps (single-file `--plot` and per-file plots under **`--heatmap-dir`**), not to the wide **combined** figure (see row below the table).

| Option | Meaning |
| ------ | ------- |
| `--cmap` | Matplotlib colour map (`cmap`). Default `viridis_r`. Other registered maps include `plasma_r`, `inferno_r`, `RdYlBu_r`, `coolwarm`, `YlOrRd`. |
| `--vmin`, `--vmax` | Colour scale limits in Å. If omitted: square matrices use positive off-diagonal RMSDs only; the combined heatmap uses positive finite values in the merged matrix. |
| `--short-heatmap-labels` | Shortens axis labels on figures only; CSV data unchanged. On the combined heatmap, row ticks use the shortened model name when this flag is set. |
| `--heatmap-diverging-center` | `none` (default): linear normalisation. `median`: `TwoSlopeNorm` with the centre at the median RMSD (suited to diverging maps). In batch mode, the median is taken from the merged CSV when that file is written; with `--no-combined`, each subdomain matrix uses its own median. |
| `--heatmap-colorbar-orientation` | `vertical` (default) or `horizontal` (colour bar below the matrix). |
| `--heatmap-y-axis-right` | Row labels on the right. With a vertical colour bar, the bar is drawn on the left to avoid overlap. |
| `--plot-hatch` | **Single-file** (`--plot`) or **batch** with **`--heatmap-dir`** only: add a second channel—**hatch patterns** on off-diagonal cells binned from **pair counts parsed from the log** (not from CSV). **LSQ** logs: **matched atoms** (Coot); **SSM** logs: **aligned residues**. **CSV** input has no per-pair counts in-file → hatch is skipped (warning if this flag is set). Legend labels state the quantity. |
| `--heatmap-hatch-bins N` | With `--plot-hatch`: number of equal-width bins over observed counts (default **4**). |
| `--heatmap-hatch-edges E0,E1,…` | With `--plot-hatch`: optional comma-separated bin edges; extended to data min/max when needed (at least two values). |
| Missing / diagonal | No value, same-structure diagonal, or combined-table gaps map to white (below the scale). The colour bar is a plain rectangle (no triangular extension). Tick values in Å are chosen for both linear and median-centred scales. |

**Combined heatmap** (`combined_rmsd_heatmap.*`): RMSD only; the second hatch channel is **not** applied there.

**Figure format (single file):** set the extension on `--plot` (`.png`, `.pdf`, or `.svg`). PNG is often the most robust on headless systems.

**Batch figure format:** `--heatmap-format` selects `png`, `pdf`, or `svg` for auto-named outputs under `--heatmap-dir` (default **svg**). The file extension on `--combined-heatmap`, when present, selects the format for that path. SVG and PDF heatmaps use vector cells with an embedded raster colour bar; PNG uses a raster matrix. Paths may use `~`; they are expanded to absolute paths internally.

`--title` sets the heatmap title. Batch mode appends the subdomain suffix on per-file plots and “(combined)” on the merged figure when using the default title pattern.

#### Batch mode (`--scan-dir`)

For a directory tree of RMSD logs where each match produces its own square CSV and optionally one merged wide table.

1. **Recursive search** — Under `--scan-dir`, every file matching `--scan-glob` is collected (default `rmsd_values_*.txt` for LSQ logs; use `rmsd_SSM_values*.txt` for SSM logs).
2. **One square CSV per match** — Each file is parsed and written as `rmsd_table_<suffix>.csv` in the same directory. The suffix is the subdomain label for that block in the combined table.
3. **Combined table** — Unless `--no-combined` is set, `combined_rmsd_table.csv` is written (default: `<scan-dir>/combined_rmsd_table.csv`, or `--combined-csv`). Rows are sorted by `Subdomain` then `Model` (natural sort), unless `--combined-subdomain-order` fixes block order (comma-separated tokens or a file of tokens; each token matches a subdomain label with an underscore boundary, e.g. `D1` matches `…_D1s` but not `…_D10s`; unmatched subdomains follow in natural order). Structure columns follow the merged row order (first-seen model order, then any remaining names). Cells are empty where a subdomain has no value for that column.

If two logs in different directories share the same basename, they receive the same `Subdomain` string; use distinct log names when groups must stay separate.

**Row/column order in batch** — `--order`, `--order-dir`, and `--order-glob` apply to every matched log, with the same semantics as single-file mode (`reorder_matrix` skips missing labels per file, then appends the rest).

**Heatmaps in batch** — With `--heatmap-dir`, the script writes `rmsd_heatmap_<suffix>.<fmt>` for each log and, by default, `combined_rmsd_heatmap.<fmt>` from the merged CSV (`<fmt>` from `--heatmap-format`). Use `--no-combined-heatmap` to skip the combined figure. Use `--combined-heatmap /path/to/project/figs/combined_rmsd_heatmap.svg` (or `.png` / `.pdf`) to set the combined path explicitly; that option can be used without `--heatmap-dir` to emit only the combined figure. Batch mode does not use `--plot`.

Unless `--no-shared-heatmap-scale` is set, per-subdomain heatmaps use the same vmin/vmax as the merged table (or combined limits when the merged file is absent). `--vmin` / `--vmax` still bound the scale when supplied.

```bash
python ranking/rmsd_to_csv.py /path/to/project/rmsd_SSM_values.txt
python ranking/rmsd_to_csv.py /path/to/project/rmsd_values.txt -o /path/to/project/rmsd_table.csv
python ranking/rmsd_to_csv.py /path/to/project/rmsd_SSM_values.txt --plot /path/to/project/heatmap.svg --order model_01,model_02,model_03
python ranking/rmsd_to_csv.py /path/to/project/rmsd_SSM_values.txt --plot /path/to/project/heatmap.pdf --cmap RdYlBu_r --heatmap-diverging-center median
python ranking/rmsd_to_csv.py /path/to/project/rmsd_values.txt --plot /path/to/project/heatmap.png --vmin 0 --vmax 5

python ranking/rmsd_to_csv.py /path/to/project/rmsd_values.txt --order-dir /path/to/project/models --order-glob 'model_*.pdb'
python ranking/rmsd_to_csv.py /path/to/project/rmsd_values.txt --order-dir /path/to/project/models

python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/
python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --scan-glob 'rmsd_SSM_values*.txt'
python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --heatmap-dir /path/to/project/figs --heatmap-format png --combined-csv /path/to/project/combined_rmsd_table.csv
python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --heatmap-dir /path/to/project/figs --combined-subdomain-order 'token_1,token_2,token_3'
python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --no-combined

python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --order-dir /path/to/project/models --order-glob 'model_*.pdb'
```

**Single-file options:** `-o` / `--output`; `--plot`; `--title`; `--cmap`, `--vmin`, `--vmax`; `--short-heatmap-labels`; `--heatmap-diverging-center`; `--heatmap-colorbar-orientation`; `--heatmap-y-axis-right`; **`--plot-hatch`**, **`--heatmap-hatch-bins`**, **`--heatmap-hatch-edges`**; `--order` / `--order-dir` / `--order-glob` (see [Heatmap colours](#heatmap-colours)).

**Batch (`--scan-dir`) options:** `--scan-glob`; `--combined-csv`; `--no-combined`; `--combined-subdomain-order`; `--heatmap-dir`; `--heatmap-format`; `--combined-heatmap`; `--no-combined-heatmap`; `--no-shared-heatmap-scale`; `--cmap`, `--vmin`, `--vmax`, `--title`; `--short-heatmap-labels`; `--heatmap-diverging-center`; `--heatmap-colorbar-orientation`; `--heatmap-y-axis-right`; **`--plot-hatch`** (with **`--heatmap-dir`** only; parses counts from each **LSQ/SSM** log), **`--heatmap-hatch-bins`**, **`--heatmap-hatch-edges`**; `--order` / `--order-dir` / `--order-glob`. For a tree from one RMSD input, run `structure_phylogeny.py` on that log or CSV after conversion.

### `structure_phylogeny.py`

Structure-based phylogenetic tree from pairwise RMSD (Coot LSQ/SSM or CSV). Uses neighbour-joining; midpoint or outgroup rooting. Optionally computes distances from PDBs via TM-align.

```bash
python ranking/structure_phylogeny.py /path/to/project/rmsd_values.txt -o /path/to/project/structure_tree.nwk
python ranking/structure_phylogeny.py /path/to/project/rmsd_values.txt -o /path/to/project/tree.nwk --plot /path/to/project/tree.pdf

# Root on outgroup
python ranking/structure_phylogeny.py /path/to/project/rmsd_values.txt -o /path/to/project/tree.nwk --root "outgroup"

# Unrooted tree
python ranking/structure_phylogeny.py /path/to/project/rmsd_values.txt -o /path/to/project/tree.nwk --unrooted --plot /path/to/project/tree.pdf

# Compute distances from PDBs (TM-align required)
python ranking/structure_phylogeny.py --from-pdb /path/to/project/models -o /path/to/project/tree.nwk --plot /path/to/project/tree.pdf
```

Options: `-o`, `--output`; `--plot PATH`; `--root LABEL`; `--unrooted`; `--from-pdb [DIR]`; `--format auto|lsq_txt|ssm_txt|csv`

Input: LSQ `rmsd_values_*.txt`, SSM `rmsd_SSM_values*.txt`, or CSV `rmsd_table_*.csv`

Dependencies: scikit-bio (recommended) or Biopython; ete3 or matplotlib for `--plot`; TM-align for `--from-pdb`

#### DALI-like mode: per-query pseudo Z-scores (`structure_phylogeny.py`)

`structure_phylogeny.py` can optionally build the tree from a DALI-like, query-centric pseudo Z-score pipeline instead of using RMSD directly.

Enable it with `--pseudo-z per_query`. Internally the mode does: RMSD → similarity → per-query z-scores → symmetrise → clamp `z+ = max(0,z)` → convert z-scores to distances (`--zdist`) for NJ/UPGMA.

```bash
# Default pseudo-Z settings (similarity = exp(-RMSD/tau), zdist = inv, tau = 2.0)
python ranking/structure_phylogeny.py rmsd_values.txt -o dali_like_tree.nwk --pseudo-z per_query

# Tune similarity kernel and z→distance transform
python ranking/structure_phylogeny.py rmsd_values.txt -o dali_like_tree.nwk --pseudo-z per_query \
  --similarity exp --tau 2.5 --zdist exp --zscale 1.0

# Also write the per-query pseudo-Z ranking table
python ranking/structure_phylogeny.py rmsd_values.txt -o dali_like_tree.nwk --pseudo-z per_query \
  --ranking-csv dali_like_ranking.csv
```

Related options:

- `--pseudo-z {off,per_query}`
- `--similarity {exp,neg_rmsd}` and `--tau` (used with `exp`)
- `--zdist {inv,exp,maxminus}` and `--zscale` (used with `exp`)
- `--ranking-csv PATH` (optional CSV output)

### `dali_phylogeny.py`

If you already have pairwise Dali Z-scores (for example from DaliLite or `dali_score.py`), `dali_phylogeny.py` can rank structures, convert Z to distances, and build a Newick tree.

Input format: a delimited text table (TSV/CSV/space-delimited) with 3 fields per row:  
`model_01  model_02  zscore`

Comment lines starting with `#` are ignored.

```bash
# Basic usage
python ranking/dali_phylogeny.py pairwise_zscores.tsv -o dali_tree.nwk

# Write ranking CSV
python ranking/dali_phylogeny.py pairwise_zscores.tsv -o dali_tree.nwk --output-ranking dali_ranking.csv

# Choose how to convert Z → distance
python ranking/dali_phylogeny.py pairwise_zscores.tsv -o dali_tree.nwk --transform exp --exp-scale 10.0
```

Options:

- `-o, --output-tree PATH`
- `--output-ranking PATH`
- `--transform {inv,maxminus,exp}`
- `--exp-scale VALUE` (only for `exp`)
- `--root LABEL` and `--no-midpoint-root`
- `--plot PATH`
- `--log [FILE]` (optional summary log)

### `dali_score.py`

Computes Dali-equivalent structural similarity scores from two or more structures. It can optionally run **DaliLite** locally (if the user specifies its path) or use biotite/alignment file for residue equivalences. Supports pairwise and all-vs-all modes with tree/dendrogram output.

With DaliLite, the script runs `import.pl` then `dali.pl --cd1/--cd2` from a short temporary work directory (not bare `--pdbfile1/--pdbfile2`). `import.pl` requires `mkdssp` at the path set in DaliLite’s `bin/mpidali.pm` (`$DSSP_EXE`; symlink your site’s `mkdssp` there if the default path is wrong). Set `DALILITE_HOME` or `--dalilite-path` to the DaliLite root. Optional `MKDSSP` can point to your `mkdssp` binary for clearer error messages.

**Pairwise:**

```bash
python ranking/dali_score.py /path/to/project/model_01.pdb /path/to/project/model_02.pdb \
  -o /path/to/project/result.csv
# With DaliLite (canonical Z-score): set DALILITE_HOME or --dalilite-path
export DALILITE_HOME=/path/to/DaliLite
python ranking/dali_score.py /path/to/project/model_01.pdb /path/to/project/model_02.pdb
```

**All-vs-all with tree:**

```bash
python ranking/dali_score.py --all-vs-all /path/to/project/models/ -o /path/to/project/pairs.csv \
  --output-tree /path/to/project/dali_tree.nwk \
  --output-plot /path/to/project/dendrogram.png \
  --output-ranking /path/to/project/ranking.csv
```

**All-vs-all directory scan:** by default, each directory argument is searched **recursively** for `*.pdb`, `*.cif`, and `*.ent`. Pass **`--root-only-dirs`** to use only the **top level** of each directory (aligned with `dalilite_superpose_scores.py` default). Structures are **natural-sorted** before pairwise comparison; the **pairwise table**, **Z-matrix**, and **tree** inputs use that order. The **ranking** CSV is ordered **by Z-score** (average then maximum), not alphabetically. Where DaliLite does not report a summary Z, **`z_score`** uses the **empirical** FoldKit normalisation (same as in `dali_score.compute_z_score`).

Options: `--dalilite-path DIR`, `--no-dalilite`, `--filter` (substring or glob on basename for directory inputs), `--root-only-dirs` (all-vs-all: no subdirectory scan); tree: `--output-tree`, `--output-plot`, `--output-matrix`, `--transform`, `--root`. **Optional summary log:** `--log [FILE]`.

#### References (ranking, scoring, phylogeny)

- **Holm, L. & Sander, C.** (1993). Protein structure comparison by alignment of distance matrices. *J. Mol. Biol.* **233**, 123–138. [DOI: 10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489)
- **Holm, L.** (2020). Dali and the persistence of protein shape. *Protein Sci.* **29**, 128–140. [DOI: 10.1002/pro.3749](https://doi.org/10.1002/pro.3749)
- **Saitou, N. & Nei, M.** (1987). The neighbour-joining method: a new method for reconstructing phylogenetic trees. *Mol. Biol. Evol.* **4**, 406–425. [DOI: 10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)
- **Zhang, Y. & Skolnick, J.** (2005). TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Res.* **33**, 2302–2309. [DOI: 10.1093/nar/gki524](https://doi.org/10.1093/nar/gki524) — TM-align (`structure_phylogeny.py --from-pdb`).

---

## Metrics (crystal packing and lattice)

A Python toolkit for analysing differences in **X-ray crystal lattice packing**, including interfaces and crystal contacts. It complements FoldKit’s **structural similarity** tools by quantifying **how molecules pack in the unit cell**, not only how similar their isolated coordinates are.

### Overview

The pipeline implements established ways to quantify crystal packing differences that are visually apparent but hard to measure by eye. Individual molecules may have low pairwise RMSD while packing differs in interfaces and symmetry contacts.

### Key features

- **Quantitative metrics:** packing density (Matthews coefficient, solvent content), interface and contact analysis.
- **Batch export:** `--compare` writes `batch_analysis_results.json` with all per-structure outputs in one file.

### Scientific background

#### Packing density and void volume

- **Matthews coefficient (Vm):** unit cell volume per molecular weight (Å³/Da); typical proteins ~1.7–3.5 Å³/Da.
- **Solvent content** and **packing efficiency** as fractions of the cell.

#### Interface analysis

Buried surface area, contact area, interface complementarity, hydrogen-bond summaries.

#### Crystal contacts

Classification (H-bond, hydrophobic, electrostatic, van der Waals), contact density, symmetry-related vs biological interfaces.

### Installation

#### Prerequisites

```bash
# Create virtual environment (recommended)
python -m venv crystal_analysis_env
source crystal_analysis_env/bin/activate  # Linux/Mac
# or
crystal_analysis_env\Scripts\activate     # Windows
```

#### Install dependencies

Install packages needed for the crystal-packing stack (for example `biopython`, `numpy`, `scipy`, `pandas`; see script imports and error messages). If you maintain a project `requirements.txt`, use `pip install -r requirements.txt`.

#### Core dependencies

- **BioPython:** structure parsing
- **NumPy, SciPy (optional), pandas**
- **R** (optional; only for documented `ranking/*.R` scripts such as `create_rmsd_heatmap.R`)

### Crystal packing scripts

#### `crystal_packing_analyser.py`

Runs packing metrics, interface analysis, and contact analysis per structure. With `--compare`, also writes `batch_analysis_results.json`.

```bash
python metrics/crystal_packing_analyser.py --input model_01.pdb
python metrics/crystal_packing_analyser.py --input dir/
python metrics/crystal_packing_analyser.py --input *.pdb --compare --output analysis_run
python metrics/crystal_packing_analyser.py --input *.pdb --sets set_a set_b --output "crystal_analysis_{}"
python metrics/crystal_packing_analyser.py --input *.pdb --sets set_a set_b --output analysis_output --dry-run
python metrics/crystal_packing_analyser.py --input model_01.pdb --chains A --output analysis_chainA
python metrics/crystal_packing_analyser.py --input expanded_assembly.pdb --reference-chain A --output ./analysis_expanded
```

Options: `--input`, `--output` (use `{}` for one dir per set), `--set`, `--sets`, `--compare` (combined JSON), `--verbose`, `--dry-run`, `--chains`, `--reference-chain` (multi-copy lattice metrics in the interface stage)

#### `packing_metrics.py`

Matthews coefficient, solvent content, packing efficiency.

```bash
python metrics/packing_metrics.py model_01.pdb
python metrics/packing_metrics.py dir/
python metrics/packing_metrics.py *.pdb -o results.txt
python metrics/packing_metrics.py *.pdb -q -o summary.txt
python metrics/packing_metrics.py *.pdb --sets set_a set_b -o "analysis_{}.txt"
python metrics/packing_metrics.py *.pdb --per-structure -o "{}_metrics.txt"
python metrics/packing_metrics.py *.pdb --per-structure --sets set_a set_b -o "{}_metrics.txt"
```

Options: `-o`, `-q` (quiet), `--per-structure`, `--set`, `--sets`

#### `lattice_packing_analyser.py`

Packing-style metrics for **symmetry-expanded / multi-copy** coordinate models (a “supercell” or lattice fragment in one Cartesian frame, typically P1). Unlike `packing_metrics.py` (unit-cell-based), this script computes lattice-wide **density** and **packing fraction** using either the PDB `CRYST1` volume or a coordinate **bounding box**.

Key options:

- **`--volume-source {auto,cryst1,bbox}`**: choose the volume definition (default: `auto`).
- **`--bbox-pad A`**: padding (Å) added to the bounding box when `--volume-source=bbox`.
- **`--expected-chains N`**: sanity check: require exactly N chains (useful for fixed supercell expansions).
- **`--allow-cryst1-mismatch`**: override the safety check that avoids using `CRYST1` volume when `--expected-chains > 1` (often `CRYST1` is the *original* unit cell, not the supercell).
- **Set selection**: same pattern set controls as other analysers: `--set` / `--sets`.

```bash
# Single supercell (stdout JSON by default)
python metrics/lattice_packing_analyser.py supercell.pdb --volume-source bbox --bbox-pad 2.0

# Require expected number of chains (helpful for fixed expansion protocols)
python metrics/lattice_packing_analyser.py supercell.pdb --expected-chains 12 --volume-source bbox

# Multiple structures: write per-structure JSON/TXT under an output directory
python metrics/lattice_packing_analyser.py supercells/ --sets set_a set_b --output-dir "lattice_pack_{}" --output-txt "lattice_{}.txt"

# Multiple structures: write one combined JSON
python metrics/lattice_packing_analyser.py supercells/*.pdb --output-json lattice_packing_combined.json
```

Output fields (high level): chosen volume metadata, atom/mass densities, packing density fraction/percent, Matthews-like ratio (Å³/Da) and heuristic solvent content (see `metrics/metrics_details.md`, Section 1.6).

#### Interface analysis decision table

| Goal | Recommended script | Notes / outputs |
| --- | --- | --- |
| Pairwise interfaces in the ASU + **charge-tag** complementarity | `interface_analyser_asu_charge.py` | Reports `charge_complementarity*` per interface. |
| Pairwise interfaces in the ASU + **electrostatic complementarity (EC; McCoy)** | `interface_analyser_asu_ec.py` | Reports `ec_*` per interface; charge-tag metrics are suppressed. |
| Lattice interfaces for a reference chain (multi-copy models) + **charge-tag** metrics | `interface_analyser_lattice_charge.py` | Requires `--reference-chain`; reports `lattice_charge_*`. |
| Lattice interfaces for a reference chain (multi-copy models) + **EC (McCoy)** | `interface_analyser_lattice_ec.py` | Requires `--reference-chain`; reports `lattice_ec_*`. |

#### `interface_analyser_asu_charge.py`

Protein–protein interface analysis (buried surface area, contacts, complementarity).

```bash
python metrics/interface_analyser_asu_charge.py model_01.pdb
python metrics/interface_analyser_asu_charge.py model_01.pdb --chains A -o chainA_interfaces.txt
```

Options: `-o` (use `{}` for per-set or per-structure), `--per-structure`, `--set`, `--sets`, `--dry-run`, `--chains`, `--reference-chain` (focal chain: isolated vs embedded SASA, lattice burial fraction, cross-chain contact-residue fraction)

The text report summary lists **isolated SASA** (Shrake–Rupley, each chain alone) for every chain that appears in any reported interface, plus the **sum** of those values (e.g. chains A, B, C when interfaces A–B and A–C are reported). With **`--reference-chain`**, the summary adds **`sasa_reference_isolated`**, **`sasa_reference_in_cluster`**, **`lattice_burial_fraction`** (approximately `1 − SASA_cluster/SASA_iso`, a geometric occlusion index, **not** thermodynamic **ΔSASA**), and **`lattice_contact_residue_fraction`** (reference residues with a cross-chain neighbour within the contact cut-off).

**Complementarity metrics (naming and methods).**

- **Charge-tag complementarity (contact-based metrics)**: reported by `interface_analyser_asu_charge.py` and `interface_analyser_lattice_charge.py` as `charge_complementarity*` (and lattice fields `lattice_charge_*`). These are derived from charged residue identities within the filtered atom–atom contact set.
- **Electrostatic complementarity (EC; McCoy method)**: reported by `interface_analyser_asu_ec.py` and `interface_analyser_lattice_ec.py` as `ec_*` (per interface) and `lattice_ec_*` (lattice summaries). EC is a Pearson correlation of electrostatic potentials on facing surface points, following McCoy, Epa & Colman (1997). Definitions and reported fields: `metrics/EC_details.md` and `metrics/metrics_details.md` (Section 2.5.4).

#### `interface_analyser_asu_ec.py`

Electrostatic complementarity (EC) for interfaces using the McCoy method (per-interface `ec_r`, `ec_n_pairs`, and density; plus lattice-weighted summaries in lattice mode).

```bash
python metrics/interface_analyser_asu_ec.py model_01.pdb -o ec_results.txt
```

#### `interface_analyser_lattice_charge.py`

```bash
python metrics/interface_analyser_lattice_charge.py expanded_assembly.pdb --reference-chain A -o lattice_charge.txt
```

#### `interface_analyser_lattice_ec.py`

```bash
python metrics/interface_analyser_lattice_ec.py expanded_assembly.pdb --reference-chain A -o lattice_ec.txt
```

#### Lattice entrypoints (reference-chain reports)

For multi-copy coordinate models, use the lattice interface entrypoints:

```bash
python metrics/interface_analyser_lattice_ec.py expanded_assembly.pdb --reference-chain A -o lattice_ec.txt
python metrics/interface_analyser_lattice_charge.py expanded_assembly.pdb --reference-chain A -o lattice_charge.txt
```

### Multi-copy models (crystallographic interfaces)

Applicable when the coordinate file contains **several chains** representing **repeated molecules** in one frame (symmetry expansion, expanded lattice fragment, oligomer with one chain per subunit). Space-group operators are **not** applied; outputs describe the supplied Cartesian geometry only. Pairwise interface metrics (BSA, contacts, complementarity) are reported **per chain pair**. Each **BSA** value uses the usual formula with a **two-chain-only** complex (isolated SASA for each chain minus SASA of those two chains together), so pairs such as A–B, A–C, and A–D are each well-defined. Summed BSA over many pairs still aggregates separate pairwise terms and is **not** a single global burial measure for one molecule in the full assembly. The focal-copy **`lattice_burial_fraction`** is a geometric occlusion index, **not** a thermodynamic **ΔSASA**. Definitions, limitations, and interpretation: **`metrics/metrics_details.md`** (Sections 2.7.1–2.7.2).

**Typical workflow**

1. Build or export a multi-chain structure with distinct chain IDs per copy.
2. Run **`interface_analyser_lattice_charge.py`** (charge-tag) or **`interface_analyser_lattice_ec.py`** (EC) on that file; use **`--chains`** to restrict to interfaces involving a focal molecule, and **`--reference-chain`** for the focal-copy occlusion summary in the report header.
3. Optionally run **`crystal_packing_analyser.py`** with the same flags so **`interface_analysis.summary`** in each `*_analysis.json` includes the reference-chain fields.

**Examples** (repository root; substitute paths and chain IDs)

```bash
# All pairwise interfaces in an expanded multi-chain model
python metrics/interface_analyser_lattice_charge.py /path/to/project/expanded_assembly.pdb --reference-chain A -o assembly_interfaces.txt

# Only interfaces where the focal chain participates (partner chains still listed per pair)
python metrics/interface_analyser_lattice_charge.py /path/to/project/expanded_assembly.pdb --reference-chain A --chains A -o interfaces_focal_A.txt

# Focal-copy summary: isolated vs embedded SASA, lattice_burial_fraction, contact-residue fraction
python metrics/interface_analyser_lattice_charge.py /path/to/project/expanded_assembly.pdb --reference-chain A --chains A -o lattice_focal_A.txt

# Full pipeline; per-structure JSON under the output directory includes the interface summary
python metrics/crystal_packing_analyser.py --input /path/to/project/expanded_assembly.pdb \
  --chains A --reference-chain A --interface-metrics charge --output ./analysis_expanded_assembly

# Full pipeline with EC (McCoy) interface metrics
python metrics/crystal_packing_analyser.py --input /path/to/project/expanded_assembly.pdb \
  --chains A --reference-chain A --interface-metrics ec --output ./analysis_expanded_assembly_ec

# Batch with combined JSON (--compare); same focal chain for every file matched by the glob
python metrics/crystal_packing_analyser.py --input '/path/to/project/models/*_assembly.pdb' \
  --reference-chain A --interface-metrics charge --compare --output ./batch_lattice_analysis
```

**Outputs.** Text: one block per chain pair (BSA, contacts, complementarity, etc.) plus summary lines (isolated SASA for chains in reported interfaces; with **`--reference-chain`**, **`sasa_reference_isolated`**, **`sasa_reference_in_cluster`**, **`lattice_burial_fraction`**, **`lattice_contact_residue_fraction`**). JSON (pipeline): nested under `interface_analysis.summary` in `*_analysis.json`.

#### `interface_molecule_report_csv.py`

Filter an interface text report into CSV (optional filters by PDB basename and chain). Input is typically the file from `interface_analyser_asu_charge.py … -o results.txt` (or the EC analogue). Writes one file per structure when using `--output-dir` or `-o` with `{}`.

```bash
python metrics/interface_molecule_report_csv.py results.txt -m A -m B --output-dir ./out
python metrics/interface_molecule_report_csv.py results.txt --chains A,B --group-by-chain --output-dir ./out
python metrics/interface_molecule_report_csv.py results.txt --pdbs model_01.pdb,model_02.pdb --chains A,B -o 'out/{}.csv'
python metrics/interface_molecule_report_csv.py results.txt --chains A,B --combine-regex '^(model\d+[^_]*)_' --output-dir ./out
python metrics/interface_molecule_report_csv.py results.txt --chains A,B --combine-glob 'run_a*' --combine-glob 'run_b*' --output-dir ./out
python metrics/interface_molecule_report_csv.py results.txt --chains A,B \
  --combine-glob 'run_a*' --combine-glob 'run_b*' --combine-glob 'run_*' --output-dir ./out
```

Options: `--pdb`, `--pdbs`, `-m`, `--molecule`, `--chains`, `-o`, `--output-dir`, `--group-by-chain`, `--combine-regex` (Python regex; **first capturing group** = merge key / output stem), `--combine-glob` (fnmatch on stem). **Order matters** for globs: list more specific patterns before broader ones (e.g. `run_a`* before `run_*`), because a single broad glob can match stems you intended for a narrower rule.

The regex above captures `model`, digits, then **any suffix without an underscore** up to the next `_` (e.g. a replicate or batch token after the first `_`). It maps stems like `model_01_rep1`, `model_01_rep2`, `model_02_rep1` to `model_01.csv`, `model_02.csv`. If stems are only `modelN_…` with no letters after the digits, `^(model\d+)_` is enough. To restrict allowed suffixes after the digits, narrow the regex (character class or alternation). To list variants explicitly, enumerate them in the alternation (e.g. `'^(variant_a|variant_b)_'`). Invalid: `model\d+`* — use `\d+` plus a separate suffix rule (`[^_]*`, etc.).

#### `contact_analyser.py`

Crystal contact and symmetry interaction analysis.

```bash
python metrics/contact_analyser.py model_01.pdb
python metrics/contact_analyser.py *.pdb -o contact_results.txt
python metrics/contact_analyser.py *.pdb --sets set_a set_b -o "contact_{}.txt"
python metrics/contact_analyser.py *.pdb --per-structure -o "{}_contact.txt"
python metrics/contact_analyser.py *.pdb --per-structure --sets set_a set_b -o "{}_contact.txt"
python metrics/contact_analyser.py model_01.pdb --chains A -o chainA_contacts.txt
```

Options: same core options as the interface analysers (including `--chains`).

#### `contact_molecule_report_csv.py`

Filter a `contact_analyser.py` text report into CSV rows (one per atom–atom contact). Columns include `chain1`, `chain2`, `res1`, `atom1`, `res2`, `atom2`, `distance_A`, `contact_type`, plus `set_label` and `structure_basename`. Same filtering and merge options as `interface_molecule_report_csv.py` (`--pdb`, `--pdbs`, `-m`, `--chains`, `-o`, `--output-dir`, `--combine-regex`, `--combine-glob`). For `*_asu_contacts.txt` sidecars (no progress lines in file), use `--structure-basename model_01.pdb` or rely on the filename heuristic (`…_<stem>_asu_contacts.txt`).

```bash
python metrics/contact_molecule_report_csv.py contact_results.txt -m A -m B --output-dir ./out
python metrics/contact_molecule_report_csv.py contact_results.txt --pdbs model_01.pdb,model_02.pdb --chains A,B -o 'out/{}.csv'
python metrics/contact_molecule_report_csv.py contact_results_model_01_asu_contacts.txt --structure-basename model_01.pdb --chains A,B -o model_01_contacts.csv
```

### Python API usage

Run Python with **`metrics/`** on the module path (for example `PYTHONPATH=metrics` from the repository root, or `cd metrics`).

```python
from crystal_packing_analyser import CrystalPackingAnalyser

analyser = CrystalPackingAnalyser(output_dir="analysis_output")
results = analyser.analyse_single_structure("model_01.pdb")

model_paths = ["model_01.pdb", "model_02.pdb", "model_03.pdb"]
all_results = [analyser.analyse_single_structure(p) for p in model_paths]
batch = analyser.compare_structures(all_results)  # writes batch_analysis_results.json
```

### Individual module usage

#### Packing metrics

```python
from packing_metrics import PackingMetricsCalculator

calc = PackingMetricsCalculator()
metrics = calc.calculate_metrics("model_01.pdb")
print(f"Matthews coefficient: {metrics['matthews_coefficient']:.2f}")
print(f"Solvent content: {metrics['solvent_content_percent']:.1f}%")
```

#### Interface analysis

```python
from interface_analyser_base import InterfaceAnalyser

analyser = InterfaceAnalyser(contact_distance=5.0)
interfaces = analyser.analyse_interfaces("model_01.pdb")
print(f"Total interfaces: {interfaces['summary']['total_interfaces']}")
```

### Output files

#### Analysis results (typical layout)

```
crystal_analysis_output/
├── model_01_analysis.json
├── model_02_analysis.json
└── batch_analysis_results.json    # when using --compare (or compare_structures in API)
```

#### JSON output (excerpt)

```json
{
  "structure_id": "model_01",
  "packing_metrics": {
    "matthews_coefficient": 2.15,
    "solvent_content_percent": 43.2,
    "packing_efficiency_percent": 74.8
  },
  "interface_analysis": {
    "summary": {
      "total_interfaces": 6,
      "total_buried_surface_area": 1250.5
    }
  },
  "contact_analysis": {
    "contact_summary": {
      "asu_stats": {
        "total_contacts": 145,
        "average_distance": 3.8
      }
    }
  }
}
```

### Recommended workflow

**With Coot:** prepare or superpose structures, export PDBs, then run `crystal_packing_analyser` on those files.

**With ChimeraX or other viewers:** align if needed, export PDBs, run the same pipeline.

**After FoldKit superposition:**

```bash
python superimposition/superimpose_coot_SSM.py reference.pdb dir1 dir2 dir3
python metrics/crystal_packing_analyser.py --input SSMaligned2_ref/*.pdb --compare --output packing_batch
```

### Scientific applications

- Polymorphism and alternate crystal forms
- Crystal engineering and contact optimisation
- Ligand effects on packing (apo vs holo)
- Biological vs crystal interfaces

### Performance considerations

- Single-structure analysis is roughly O(n²) in atom count for contact-heavy steps.
- Typical memory ~100–500 MB per structure; runtime on the order of tens of seconds to a few minutes.

### Quick start

```bash
python appendix/quick_start.py
python metrics/crystal_packing_analyser.py --input model_01.pdb
python metrics/crystal_packing_analyser.py --input model_01.pdb model_02.pdb model_03.pdb --compare
```

#### References (metrics / crystal packing)

1. **Matthews, B.W.** (1968). Solvent content of protein crystals. *J. Mol. Biol.* **33**, 491–497. [DOI: 10.1016/0022-2836(68)90205-2](https://doi.org/10.1016/0022-2836(68)90205-2)
2. **Janin, J.** (1997). Specific versus non-specific contacts in protein crystals. *Nat. Struct. Biol.* **4**, 973–974. [DOI: 10.1038/nsb1297-973](https://doi.org/10.1038/nsb1297-973)
3. **Carugo, O. & Argos, P.** (1997). Protein–protein crystal-packing contacts. *Protein Sci.* **6**, 2261–2263. [DOI: 10.1002/pro.5560061021](https://doi.org/10.1002/pro.5560061021)

*These DOIs were verified against the Crossref REST API.*

#### Troubleshooting (metrics pipeline)

1. **BioPython import errors:** `pip install biopython>=1.79`
2. **Memory:** reduce contact distance cutoffs for large structures
3. **R:** install from [r-project.org](https://www.r-project.org/) only if you run documented `ranking/*.R` scripts; on macOS ensure R is on `PATH` (e.g. `export PATH="/usr/local/bin:$PATH"`).

#### Features summary

- Modular analysers usable alone or via `crystal_packing_analyser`
- Metrics aligned with common crystallographic practice
- Optional R utilities (e.g. RMSD heatmaps) where listed under superposition / ranking sections
- Fits naturally after superposition or alongside external refinement workflows

---

## Appendix

### General

- Scripts accept **PDB and mmCIF** where not otherwise noted.
- **Coot** superposition: single-set one-to-many and single-set all-vs-all **keep Coot open** by default; use `**--not-interactive`** to exit when done. **AxB** (two-set) mode is **non-interactive by default**; use `**--interactive`** to keep Coot open after reload.
- Run `**extract_rmsd.py**` on `coot_log.txt` with `**--format auto**` to classify SSM vs LSQ logs.

### Logging (optional summary logs)

Most command-line scripts accept **`--log [FILE]`** to write a short **summary log** (tasks run + errors). This is **opt-in** and it does **not** mirror stdout/stderr.

- **`--log`** (no value): write `<script_stem>.log` in the current working directory.
- **`--log myrun.log`**: write to that file (directories are created if needed).

For Coot superposition scripts, `--log` is independent of the established Coot logs (`coot_log*.txt`) and RMSD outputs.

**Coot superposition scripts** (`superimpose_coot_SSM.py`, `superimpose_coot_LSQ.py`, `trim_superimposeLSQ.py`) also write their established Coot logs (`coot_log*.txt`) and RMSD outputs; the standard log options above are supported in a compatible way and do not replace `coot_log*.txt`.

### Output directory structure (Coot superposition)

```
working_directory/
├── SSMaligned2_[reference_name]/     # SSM one-to-many
│   ├── [model]_SSMaligned2_[reference].pdb
│   ├── coot_log.txt
│   └── rmsd_SSM_values.txt
├── SSMaligned_all_vs_all/            # SSM all-vs-all
├── LSQaligned2_[reference_name]/     # LSQ one-to-many
├── LSQaligned_all_vs_all/            # LSQ all-vs-all (superimpose_coot_LSQ)
├── LSQ_all_vs_all/                   # LSQ all-vs-all (trim_superimposeLSQ)
├── LSQ_set_a_all_vs_all/            # LSQ all-vs-all with filter
├── trimmed_set_a/                    # Trimmed models (trim_models.py / trim_superimposeLSQ --trim)
├── trimmed_set_b/
├── trimmed_all/                      # Trimmed models, no --filter (same)
└── [subdir]_LSQaligned_[ref]/        # superimpose_coot_LSQ.py --pattern
```

### Troubleshooting

- **Coot:** installed and on `PATH`; reference and model files exist; working directory writable.
- **CIF:** install `gemmi` if conversion issues appear.

#### DaliLite (`dali_score.py`, `dalilite_superpose_scores.py`)

- **`DALILITE_HOME`** or **`--dalilite-path`** must point at the DaliLite install root (`bin/dali.pl`).
- **`mkdssp`** must exist at the path in **`bin/mpidali.pm`** (`$DSSP_EXE`). If imports produce empty **`DAT/*.dat`**, install or symlink **`mkdssp`** or edit **`$DSSP_EXE`**. Perl files under **`bin/`** must be readable (fix permissions if you see `Can't locate FSSP.pm: Permission denied`).
- **Path length:** DaliLite Fortran limits full path length (~80 characters); these scripts stage under the system temporary directory automatically.
- **All-vs-all file discovery:** **`dalilite_superpose_scores.py`** collects structures from the **root** of each given directory unless **`--recursive`** is set. **`dali_score.py`** defaults to a **recursive** directory scan; use **`--root-only-dirs`** for root-level only.
- **Z in tables:** When DaliLite writes a hit line, the **`z_score`** column is **DaliLite’s reported Z**; otherwise it is **FoldKit’s empirical Z** (the script labels **empirical Z-score** in the log when the latter is used). **`raw_score`** in FoldKit is always the **recomputed** Dali sum over the mapped core; **`lali` / `nres` / `%id`** in CSV come from DaliLite’s summary when present.
- **Table order (dalilite_superpose_scores):** pairs CSV and Z-matrix use **natural-sorted** labels; the ranking file is **sorted by Z** only. **`--equivalences-dir`** writes per-pair TSVs of aligned residues (see the **`dalilite_superpose_scores.py`** subsection).
- **Z heatmap:** **`dalilite_superpose_scores.py`** can write a matplotlib Z-score heatmap with **`--heatmap PATH`** using **`ranking/foldkit_heatmap.py`** (same style options as **`rmsd_to_csv.py`**: `--cmap`, `--vmin`, `--vmax`, `--heatmap-diverging-center`, etc.). Optional second channel: **`--heatmap-n-core-patterns`** (with **`--heatmap-n-core-bins`** / **`--heatmap-n-core-edges`**) encodes **`n_core`** as hatch while colour remains **Dali Z**; see [Superimposition — `dalilite_superpose_scores.py`](#dalilite_superpose_scorespy).

### Notes

- SSM superposition reports RMSD and alignment statistics in the log.
- LSQ superposition is standard least-squares Cα fitting.
- Trimming aligns residue ranges before LSQ for fair comparison.
- Original input files are not modified by superposition scripts; outputs go to new directories.
- RMSD extraction parses alignment and RMSD lines from Coot logs.

## License

FoldKit is licensed under the [Apache License, Version 2.0](LICENSE).