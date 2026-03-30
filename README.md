# FoldKit — protein structure analysis scripts

**FoldKit** is a collection of Python (and supporting R) tools for working with macromolecular **3D structures** in PDB or mmCIF format. The scripts cover **file preparation**, **pairwise and batch superimposition**, **similarity scoring and structure-based phylogenies** (with tabular and graphical outputs), and **quantitative crystal-packing and lattice metrics**. Some workflows use **Coot** for interactive superposition and log-based RMSD extraction; others use **DaliLite**, **TM-align**, or pure Python/R and do not require Coot.

All scripts live under `FoldKit/`. From the repository root run `python FoldKit/script_name.py` (or `Rscript FoldKit/...` where noted).

**Example naming (this file, `FoldKit/*.md`, and script `--help` / docstrings):** `model_01`, `model_02`, … denote structure files or row/column labels; `set_a`, `set_b`, … denote filter substrings, `--sets` groups, or `--tags` batch tokens; `condition_1`, `condition_2`, … denote condition subfolders in batch drivers; `ref_id` / `model_id` are filename tokens in LSQ `--pattern` mode; `ref_01` illustrates a reference name in Coot log filters; `results.txt` denotes a merged text report from tools such as `interface_analyzer.py -o` or `packing_metrics.py -o`; `./out` denotes a generic output directory; chain IDs in examples (`A`, `B`, …) are placeholders.

## Contents

| Section | What it covers |
|--------|----------------|
| [File management](#file-management) | Renaming files; PDB chain ID replacement and merge; residue-range trimming (no Coot) |
| [Superimposition](#superimposition) | Coot SSM/LSQ and trim workflows; batch LSQ; DaliLite superposed coordinates; opening models in Coot |
| [Ranking, scoring, phylogeny, and graphical outputs](#ranking-scoring-phylogeny-and-graphical-outputs) | RMSD from logs; RMSD tables and heatmaps; Dali Z-scores; neighbor-joining trees and plots |
| [Metrics (crystal packing and lattice)](#metrics-crystal-packing-and-lattice) | Packing density, interfaces, contacts, channels, graphs, comparative analysis, R visualization |
| [Appendix](#appendix) | General notes, output layout, troubleshooting |

## Prerequisites

**Common:** Python 3; structure files in PDB or mmCIF where applicable.

**Depending on what you run:**

- **Coot** (command-line): for `superimpose_coot_*.py`, `trim_superimposeLSQ.py` (superposition modes), `open_models_in_coot.py`, and log parsing in `extract_rmsd.py`. Not needed for **`trim_models.py`** (trim-only).
- **R** (optional): `create_rmsd_heatmap.R` and the crystal-packing visualization pipeline.
- **DaliLite** + **`mkdssp`**: `dalilite_superpose_scores.py` and DaliLite-backed modes of `dali_score.py` (see [Appendix — DaliLite](#dalilite-dali_scorepy-dalilite_superpose_scorespy)).
- **TM-align**: optional distance input for `structure_phylogeny.py --from-pdb`.
- **`gemmi`** (`pip install gemmi`): helpful for some CIF handling in superposition workflows.

---

## File management

### `rename_files.py`

**Legacy (two arguments):** in a single directory, renames PDBs matching `{prefix}_YYYY_MM_DD_HH_MM_*.pdb` to `{prefix}_*.pdb` (non-recursive).

```bash
python FoldKit/rename_files.py DIRECTORY PREFIX
```

**Regex mode:** remove or replace patterns in basenames; by default processes subdirectories and may rename directories unless `--files-only` / `--no-recursive`.

```bash
python FoldKit/rename_files.py DIRECTORY --remove='REGEX'
python FoldKit/rename_files.py DIRECTORY --replace='REGEX' --with='STRING'
```

Options (regex mode): `--remove`, `--replace`, `--with` (required with `--replace`), `--no-recursive`, `--files-only`

Examples:

```bash
python FoldKit/rename_files.py /path/to/dir sample_prefix
python FoldKit/rename_files.py /path/to/dir --remove='\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_'
python FoldKit/rename_files.py /path/to/dir --replace='fold_' --with='protein_'
python FoldKit/rename_files.py /path/to/dir --remove='temp_' --no-recursive
```

### `pdb_rechain.py`

Replace chain IDs on ATOM/HETATM/TER lines, or merge one chain into another and renumber residues so the source chain continues after the last residue on the target chain (e.g. B → A for a single continuous chain).

```bash
# Copy all matching PDBs into out/, chain B → A with residue numbers after last A:
python FoldKit/pdb_rechain.py /path/to/input_dir --pattern '*_models.pdb' -f B -t A \
  --merge-renumber -o /path/to/output_dir/

# Simple rename (only if target chain is empty in the file):
python FoldKit/pdb_rechain.py model_01.pdb -f X -t Y -o model_01_Y.pdb
```

`SEQRES` / `CRYST1` lines are not rewritten; adjust those separately if your downstream tool requires them.

### `trim_models.py`

**Residue trimming only** (no superposition, no Coot): harmonize PDB/mmCIF models to the **shortest** residue span found among the inputs. **Self-contained** (stdlib + optional **`gemmi`** for mmCIF→PDB); it does not import other FoldKit scripts. **`trim_superimposeLSQ.py`** uses the same trimming implementation via **`trim_models.py`** for `--trim` / `--trim-only`.

```bash
python FoldKit/trim_models.py [--filter=set_a,set_b,...] directory1 [directory2 ...]
```

**`--filter`:** optional comma-separated patterns (basename substring or glob). If omitted, every `*.pdb` / `*.cif` in the given directories is included (output under `trimmed_all/`).

Examples:

```bash
python FoldKit/trim_models.py models/
python FoldKit/trim_models.py --filter=set_a models/
python FoldKit/trim_superimposeLSQ.py --trim-only --filter=set_a models/
```

Output: `trimmed_<pattern>/` per filter, or `trimmed_all/` with no filter; trimmed PDBs. **`gemmi`** may be needed for mmCIF inputs.

#### References (file management)

- **wwPDB:** [PDBx/mmCIF and PDB archive documentation](https://www.wwpdb.org/documentation) — file layout and atom-record conventions.
- **Bernstein, F.C. et al.** (1977). The Protein Data Bank: a computer-based archival file for macromolecular structures. *J. Mol. Biol.* 112, 535–542. [DOI: 10.1016/0022-2836(77)90297-6](https://doi.org/10.1016/0022-2836(77)90297-6)

---

## Superimposition

These scripts align models to a reference or run all-vs-all jobs. **Coot-based** tools leave Coot running for inspection; **`dalilite_superpose_scores.py`** writes superposed coordinates using **DaliLite** (separate install).

### `superimpose_coot_SSM.py`

Secondary Structure Matching (SSM) superposition: one-to-many or all-vs-all.

**Default (no chain flags):** uses the **first chain** in each structure and SSM mode `1` (legacy Coot selection).

**Explicit chains:** pass **`--ref-chain`** and/or **`--model-chain`** to use mmdb-style selections (`/1/CHAIN/*`) and SSM with move flag `0`. If only one of the two is set, the other defaults to **`A`**.

```bash
# One-to-many: align all models to a reference
python FoldKit/superimpose_coot_SSM.py [--filter=set_a] [--output-dir=DIR] reference.pdb dir1 [dir2 ...]

# All-vs-all: each structure used as reference in turn
python FoldKit/superimpose_coot_SSM.py [--filter=set_a] [--output-dir=DIR] --all-vs-all dir1 [dir2 ...]

# Explicit chains (one-to-many or all-vs-all)
python FoldKit/superimpose_coot_SSM.py --ref-chain=B --model-chain=B reference.pdb models/
python FoldKit/superimpose_coot_SSM.py --all-vs-all --ref-chain=B --model-chain=B models/
```

Examples:

```bash
python FoldKit/superimpose_coot_SSM.py reference.pdb dir1 dir2 dir3
python FoldKit/superimpose_coot_SSM.py --all-vs-all --filter=set_a models/
python FoldKit/superimpose_coot_SSM.py --filter=set_a --output-dir=SSM_run reference.pdb models/
python FoldKit/superimpose_coot_SSM.py --ref-chain=C --model-chain=B reference.pdb models/
```

Options: `--filter` (substring or glob on basename, e.g. `set_a`), `--output-dir=DIR` (placeholders: `[reference_name]`, `[filter]`), `--ref-chain`, `--model-chain`, `--all-vs-all`

Output: `SSMaligned2_[reference_name]` or `SSMaligned_all_vs_all[*]`; `[model]_SSMaligned2_[reference].pdb`; `coot_log.txt`; `rmsd_SSM_values.txt` (via `extract_rmsd.py --format ssm` or `--format auto`). Prior outputs whose paths contain `SSMaligned` or `LSQaligned` are skipped when collecting models.

### `superimpose_coot_LSQ.py`

Least Squares (LSQ) superposition: one-to-many or all-vs-all. Uses first chain only.

```bash
# One-to-many: specify reference explicitly
python FoldKit/superimpose_coot_LSQ.py [--output-dir=DIR] --reference=reference.pdb dir1 [dir2 ...]

# One-to-many: interactive reference selection from directories
python FoldKit/superimpose_coot_LSQ.py [--output-dir=DIR] dir1 [dir2 ...]

# All-vs-all
python FoldKit/superimpose_coot_LSQ.py --all-vs-all [--output-dir=DIR] dir1 [dir2 ...]
```

Examples:

```bash
python FoldKit/superimpose_coot_LSQ.py --reference=reference.cif --filter=set_a models/
python FoldKit/superimpose_coot_LSQ.py --all-vs-all --filter=set_a models/
python FoldKit/superimpose_coot_LSQ.py dir1 dir2
```

Options: `--reference=FILE`, `--filter` (substring match, e.g. `set_a`), `--output-dir=DIR`

Output: `LSQaligned2_[reference_name]` or `LSQaligned_all_vs_all[*]`; `coot_log.txt`; `rmsd_values.txt` (via `extract_rmsd.py`)

**Pattern mode (`--pattern`):** pair reference and model files by **filenames** (and optional subdirectory names), then run the same LSQ Coot job per match. This replaces the old standalone `trim_superimposeLSQ_pattern.py` script.

```bash
python FoldKit/superimpose_coot_LSQ.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]
```

Options:

- `--strict-position`: ref_pattern before divider, target_pattern after
- `--divider=STRING`: divider string inside reference filenames
- `--ref-file-pattern=GLOB` (default: `"*.pdb"`)
- `--model-file-pattern=GLOB` (default: `"*.cif"`)
- `--output-suffix=STRING`: output directory infix (default: `"_LSQaligned_"`)

Examples:

```bash
python FoldKit/superimpose_coot_LSQ.py --pattern /path/to/refs /path/to/models ref_id model_id
python FoldKit/superimpose_coot_LSQ.py --pattern --ref-file-pattern="*final.pdb" --model-file-pattern="*.pdb" /path/to/refs /path/to/models ref_id model_id
python FoldKit/superimpose_coot_LSQ.py --pattern --divider=LSQaligned2 /path/to/refs /path/to/models ref_id model_id
```

Output: `[subdir]_LSQaligned_[ref_name]/`; aligned PDBs named `[model]_LSQaligned_[ref].pdb` inside each output directory.

### `trim_superimposeLSQ.py`

Trims models to a **common residue span** (shortest model among the set, optional `--filter` groups), then performs LSQ superposition—**or** trim only without Coot. Trimming is implemented in **`trim_models.py`**; use that script for a standalone trim-only run, or **`--trim-only`** here for the same logic before superposition.

```bash
# Trim and superimpose (all-vs-all or one-to-many with interactive reference)
python FoldKit/trim_superimposeLSQ.py --trim [--all-vs-all] [--filter=set_a,set_b] dir1 [dir2 ...]

# Only trim (no Coot) — same implementation as trim_models.py
python FoldKit/trim_superimposeLSQ.py --trim-only [--filter=set_a,set_b,...] dir1 [dir2 ...]

# Superimpose without trimming
python FoldKit/trim_superimposeLSQ.py [--all-vs-all] [--filter=set_a] dir1 [dir2 ...]
```

Examples:

```bash
python FoldKit/trim_superimposeLSQ.py --trim --all-vs-all --filter=set_a,set_b models/
python FoldKit/trim_superimposeLSQ.py --all-vs-all models/
python FoldKit/trim_superimposeLSQ.py --trim dir1 dir2
python FoldKit/trim_models.py --filter=set_a dir1
python FoldKit/trim_superimposeLSQ.py --trim-only --filter=set_a dir1
```

Options: `--trim`, `--trim-only` (implies `--trim`; incompatible with `--all-vs-all`), `--all-vs-all`, `--filter=set_a[,set_b,...]` (comma-separated substrings)

Output: With trimming, creates `trimmed_set_a/`, `trimmed_set_b/` (or `trimmed_all/`); superposition output in `LSQ_set_a_all_vs_all/`, `LSQ_all_vs_all/`, or `LSQaligned2_[ref]` (one-to-many).

### `run_all_superpositions.py`

Batch driver for **`superimpose_coot_LSQ.py --pattern`**: loops over **conditions** (subfolders under each base, e.g. `condition_1`) and **tags** (same role as `--filter=set_a` in SSM/LSQ). Expects `refs_base/<condition>/<REFERENCE_SUBDIR>/` for reference PDBs (default subdir `LSQaligned2_reference_m0` in the script) and `models_base/<condition>/*_fl_<set>*` for models. The script passes **`--pattern`**, **`--divider=LSQ_`**, and globs derived from **`CONDITION_PREFIX`** (e.g. `run1_sd_set_a` vs `run1_sd_reference`); edit **`CONDITION_PREFIX`**, **`REFERENCE_SUBDIR`**, **`REFERENCE_STEM`**, and **`TAGS`** in the script to match your layout.

```bash
python FoldKit/run_all_superpositions.py --tags set_a set_b \
  --ref-base /path/to/refs --models-base /path/to/models

# Subset of conditions (default: all keys in CONDITION_PREFIX)
python FoldKit/run_all_superpositions.py --tags set_a --conditions condition_1 condition_2 \
  --ref-base /path/to/refs --models-base /path/to/models
```

Options: `--tags` (required unless `TAGS` is set in the script), `--conditions`, `--ref-base`, `--models-base`.

### `dalilite_superpose_scores.py`

Pairwise or **all-vs-all** DaliLite runs with **superposed PDB** output (target on query frame) plus CSV / Z-matrix / Newick tree. Requires DaliLite and **mkdssp** (same constraints as `dali_score.py`; see [Appendix — DaliLite](#dalilite-dali_scorepy-dalilite_superpose_scorespy)).

**Session log:** by default, stdout and stderr are also copied to `foldkit_dalilite.log` under `-d` (override with `--log PATH`, or pass `--no-log` to disable).

```bash
python FoldKit/dalilite_superpose_scores.py model_01.pdb model_02.pdb -d ./out \
  --dalilite-path /path/to/DaliLite
```

### `open_models_in_coot.py`

Opens PDB/mmCIF in Coot with C-alpha representation. **Directory mode:** pass one or more directories (top-level `*.pdb`, `*.cif`, `*.mmcif`); optional `--filter SUBSTRING` on basename. **Pattern mode:** pass glob patterns, with optional directories to search (recursive).

```bash
python FoldKit/open_models_in_coot.py [--filter SUBSTRING] dir1 [dir2 ...]
python FoldKit/open_models_in_coot.py pattern [pattern ...] [dir1] [dir2 ...]
```

Examples:

```bash
python FoldKit/open_models_in_coot.py models/
python FoldKit/open_models_in_coot.py --filter set_a dir1 dir2
python FoldKit/open_models_in_coot.py '*set_a*LSQaligned2*ref*pdb'
python FoldKit/open_models_in_coot.py '*model_*.cif' models/
```

#### References (superimposition)

- **Emsley, P. et al.** (2010). Features and development of Coot. *Acta Crystallogr. D* **66**, 486–501. [DOI: 10.1107/S0907444910007493](https://doi.org/10.1107/S0907444910007493)
- **Holm, L. & Sander, C.** (1993). Protein structure comparison by alignment of distance matrices. *J. Mol. Biol.* **233**, 123–138. [DOI: 10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489) — Dali / DaliLite distance geometry (see also Holm 2020 below).
- **Holm, L.** (2020). Dali and the persistence of protein shape. *Protein Sci.* **29**, 128–140. [DOI: 10.1002/pro.3749](https://doi.org/10.1002/pro.3749)

---

## Ranking, scoring, phylogeny, and graphical outputs

Turn superposition logs or similarity tables into **tables, heatmaps, rankings, and trees** (Newick, PDF/PNG plots where supported).

### `extract_rmsd.py`

Extracts RMSD blocks from Coot logs. **LSQ** logs use “Aligning … to …” (optional `--aligned` / `--reference` filters). **SSM** logs use “Superposing … onto …”. Use `--format lsq`, `--format ssm`, or **`--format auto`** (default) to detect from the log.

```bash
python FoldKit/extract_rmsd.py path/to/coot_log.txt [--format auto|lsq|ssm] [LSQ options]
```

LSQ-only options: `--aligned`, `-a`; `--reference`, `-r`; `--case-sensitive`, `-c`; `--debug`, `-d`.

Examples:

```bash
python FoldKit/extract_rmsd.py LSQaligned2_ref/coot_log.txt --format lsq
python FoldKit/extract_rmsd.py SSMaligned2_ref/coot_log.txt --format ssm
python FoldKit/extract_rmsd.py coot_log.txt   # auto: picks SSM vs LSQ from log content
python FoldKit/extract_rmsd.py coot_log.txt --aligned=set_a --reference=ref_01 --format lsq
python FoldKit/extract_rmsd.py coot_log.txt --debug --format lsq
```

Output: `rmsd_values*.txt` (LSQ) or `rmsd_SSM_values*.txt` (SSM).

### `create_rmsd_table.py`

Creates CSV tables from RMSD value files in a directory.

```bash
python FoldKit/create_rmsd_table.py base_directory
```

Examples:

```bash
python FoldKit/create_rmsd_table.py SSMaligned_all_vs_all/
```

Output: `rmsd_table_[subdomain].csv`, `combined_rmsd_table.csv`

### `create_rmsd_heatmap.R`

Creates heatmap PDFs from RMSD tables.

```bash
Rscript FoldKit/create_rmsd_heatmap.R base_directory
Rscript FoldKit/create_rmsd_heatmap.R base_directory viridis
Rscript FoldKit/create_rmsd_heatmap.R base_directory YlOrRd
```

Usage: `Rscript FoldKit/create_rmsd_heatmap.R <base_directory> [palette]`  
Palette: RdYlBu (default), RdYlGn, YlOrRd, viridis, plasma, etc. (RColorBrewer or viridis options)

Requirements: R with viridis, pheatmap, RColorBrewer

Output: `rmsd_heatmap_[subdomain].pdf`, `combined_rmsd_heatmap.pdf`, `combined_rmsd_heatmap_clustered.pdf`

### `rmsd_to_csv.py`

Converts pairwise RMSD file (all-vs-all SSM/LSQ) to square CSV table. Optionally plots heatmap (matplotlib, no R).

```bash
python FoldKit/rmsd_to_csv.py rmsd_SSM_values.txt
python FoldKit/rmsd_to_csv.py rmsd_values.txt -o rmsd_table.csv
python FoldKit/rmsd_to_csv.py rmsd_SSM_values.txt --plot heatmap.pdf --order model_01,model_02,model_03
python FoldKit/rmsd_to_csv.py rmsd_SSM_values.txt --plot heatmap.pdf --cmap RdYlBu_r
python FoldKit/rmsd_to_csv.py rmsd_values.txt --plot heatmap.png --vmin 0 --vmax 5

# Row/column order from a directory (recursive); each label is the file stem and must match the RMSD file
python FoldKit/rmsd_to_csv.py rmsd_values.txt --order-dir /path/to/root --order-glob 'model_*.pdb'
python FoldKit/rmsd_to_csv.py rmsd_values.txt --order-dir /path/to/root   # glob defaults to *.pdb
```

Options: `-o`, `--output`; `--plot PATH`; `--title`; `--vmin` / `--vmax` (optional heatmap color scale in Å; defaults: min/max of **positive off-diagonal** pairwise RMSDs, so diagonal zeros do not compress the scale); `--order FILE_OR_LIST` (comma-separated list or one label per line in a file); `--order-dir DIR` with optional `--order-glob PATTERN` (default `*.pdb`) to build order from matching files under `DIR` recursively—sorted by relative path with natural numeric order; `--order` and `--order-dir` cannot be used together; `--cmap` (matplotlib colormap, default: viridis_r; e.g. plasma_r, RdYlBu_r, coolwarm, YlOrRd)

### `structure_phylogeny.py`

Structure-based phylogenetic tree from pairwise RMSD (Coot LSQ/SSM or CSV). Uses neighbor-joining; midpoint or outgroup rooting. Optionally computes distances from PDBs via TM-align.

```bash
python FoldKit/structure_phylogeny.py rmsd_values.txt -o structure_tree.nwk
python FoldKit/structure_phylogeny.py rmsd_values.txt -o tree.nwk --plot tree.pdf

# Root on outgroup
python FoldKit/structure_phylogeny.py rmsd_values.txt -o tree.nwk --root "outgroup"

# Unrooted tree
python FoldKit/structure_phylogeny.py rmsd_values.txt -o tree.nwk --unrooted --plot tree.pdf

# Compute distances from PDBs (TM-align required)
python FoldKit/structure_phylogeny.py --from-pdb /path/to/models -o tree.nwk --plot tree.pdf
```

Options: `-o`, `--output`; `--plot PATH`; `--root LABEL`; `--unrooted`; `--from-pdb [DIR]`; `--format auto|lsq_txt|ssm_txt|csv`

Input: LSQ `rmsd_values_*.txt`, SSM `rmsd_SSM_values*.txt`, or CSV `rmsd_table_*.csv`

Dependencies: scikit-bio (recommended) or Biopython; ete3 or matplotlib for `--plot`; TM-align for `--from-pdb`

#### DALI-like mode: per-query pseudo Z-scores (`structure_phylogeny.py`)

`structure_phylogeny.py` can optionally build the tree from a DALI-like, query-centric pseudo Z-score pipeline instead of using RMSD directly.

Enable it with `--pseudo-z per_query`. Internally the mode does: RMSD → similarity → per-query z-scores → symmetrize → clamp `z+ = max(0,z)` → convert z-scores to distances (`--zdist`) for NJ/UPGMA.

```bash
# Default pseudo-Z settings (similarity = exp(-RMSD/tau), zdist = inv, tau = 2.0)
python FoldKit/structure_phylogeny.py rmsd_values.txt -o dali_like_tree.nwk --pseudo-z per_query

# Tune similarity kernel and z→distance transform
python FoldKit/structure_phylogeny.py rmsd_values.txt -o dali_like_tree.nwk --pseudo-z per_query \
  --similarity exp --tau 2.5 --zdist exp --zscale 1.0

# Also write the per-query pseudo-Z ranking table
python FoldKit/structure_phylogeny.py rmsd_values.txt -o dali_like_tree.nwk --pseudo-z per_query \
  --ranking-csv dali_like_ranking.csv
```

Related options:

- `--pseudo-z {off,per_query}`
- `--similarity {exp,neg_rmsd}` and `--tau` (used with `exp`)
- `--zdist {inv,exp,maxminus}` and `--zscale` (used with `exp`)
- `--ranking-csv PATH` (optional CSV output)

### `dali_phylogeny.py`

If you already have pairwise DALI Z-scores (for example from DaliLite or `dali_score.py`), `dali_phylogeny.py` can rank structures, convert Z to distances, and build a Newick tree.

Input format: a delimited text table (TSV/CSV/space-delimited) with 3 fields per row:  
`model_01  model_02  zscore`

Comment lines starting with `#` are ignored.

```bash
# Basic usage
python FoldKit/dali_phylogeny.py pairwise_zscores.tsv -o dali_tree.nwk

# Write ranking CSV
python FoldKit/dali_phylogeny.py pairwise_zscores.tsv -o dali_tree.nwk --output-ranking dali_ranking.csv

# Choose how to convert Z → distance
python FoldKit/dali_phylogeny.py pairwise_zscores.tsv -o dali_tree.nwk --transform exp --exp-scale 10.0
```

Options:

- `-o, --output-tree PATH`
- `--output-ranking PATH`
- `--transform {inv,maxminus,exp}`
- `--exp-scale VALUE` (only for `exp`)
- `--root LABEL` and `--no-midpoint-root`
- `--plot PATH`
- `--log FILE` (default: `foldkit_dali_phylogeny.log` next to `-o`; `--no-log` disables)

### `dali_score.py`

Computes Dali-equivalent structural similarity scores from two or more structures. It can optionally run **DaliLite** locally (if the user specifies its path) or use biotite/alignment file for residue equivalences. Supports pairwise and all-vs-all modes with tree/dendrogram output.

With DaliLite, the script runs **`import.pl`** then **`dali.pl --cd1/--cd2`** from a short temporary work directory (not bare `--pdbfile1/--pdbfile2`). **`import.pl` requires `mkdssp`** at the path set in DaliLite’s **`bin/mpidali.pm`** (`$DSSP_EXE`; symlink your site’s `mkdssp` there if the default path is wrong). Set **`DALILITE_HOME`** or **`--dalilite-path`** to the DaliLite root. Optional **`MKDSSP`** can point to your `mkdssp` binary for clearer error messages.

**Pairwise:**

```bash
python FoldKit/dali_score.py model_01.pdb model_02.pdb -o result.csv
# With DaliLite (canonical Z-score): set DALILITE_HOME or --dalilite-path
export DALILITE_HOME=/path/to/DaliLite
python FoldKit/dali_score.py model_01.pdb model_02.pdb
```

**All-vs-all with tree:**

```bash
python FoldKit/dali_score.py --all-vs-all models_dir/ -o pairs.csv \
  --output-tree dali_tree.nwk --output-plot dendrogram.png --output-ranking ranking.csv
```

Options: `--dalilite-path DIR`, `--no-dalilite`, `--filter` (substring or glob on basename for directory inputs); tree: `--output-tree`, `--output-plot`, `--output-matrix`, `--transform`, `--root`. **Session log:** default `foldkit_dali_score.log` in the current directory (`--log FILE`, `--no-log`).

#### References (ranking, scoring, phylogeny)

- **Holm, L. & Sander, C.** (1993). Protein structure comparison by alignment of distance matrices. *J. Mol. Biol.* **233**, 123–138. [DOI: 10.1006/jmbi.1993.1489](https://doi.org/10.1006/jmbi.1993.1489)
- **Holm, L.** (2020). Dali and the persistence of protein shape. *Protein Sci.* **29**, 128–140. [DOI: 10.1002/pro.3749](https://doi.org/10.1002/pro.3749)
- **Saitou, N. & Nei, M.** (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Mol. Biol. Evol.* **4**, 406–425. [DOI: 10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)
- **Zhang, Y. & Skolnick, J.** (2005). TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Res.* **33**, 2302–2309. [DOI: 10.1093/nar/gki524](https://doi.org/10.1093/nar/gki524) — TM-align (`structure_phylogeny.py --from-pdb`).

---

## Metrics (crystal packing and lattice)

A Python toolkit for analyzing differences in **X-ray crystal lattice packing**, including interfaces, contacts, void spaces, and graph-theoretical descriptors. It complements FoldKit’s **structural similarity** tools by quantifying **how molecules pack in the unit cell**, not only how similar their isolated coordinates are.

### Overview

The pipeline implements established ways to quantify crystal packing differences that are visually apparent but hard to measure by eye. Individual molecules may have low pairwise RMSD while packing differs in interfaces, solvent channels, and symmetry contacts.

### Key features

- **Quantitative metrics:** packing density (Matthews coefficient, solvent content), interface and contact analysis, solvent channels, graph-based descriptors.
- **Comparative analysis:** PCA, clustering, correlations, outlier detection across structures.
- **Visualization:** R-generated comparative plots, heatmaps, PCA and clustering figures (via `visualization.py` and generated R scripts).

### Scientific background

#### Packing density and void volume

- **Matthews coefficient (Vm):** unit cell volume per molecular weight (Å³/Da); typical proteins ~1.7–3.5 Å³/Da.
- **Solvent content** and **packing efficiency** as fractions of the cell.

#### Interface analysis

Buried surface area, contact area, interface complementarity, hydrogen-bond summaries.

#### Crystal contacts

Classification (H-bond, hydrophobic, electrostatic, van der Waals), contact density, symmetry-related vs biological interfaces.

#### Solvent channels

Connectivity, bottlenecks, accessible void volume.

#### Graph-theoretical descriptors

Residue contact networks, modularity, small-world metrics, centrality.

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

Install packages needed for the crystal-packing stack (for example `biopython`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `networkx`; see script imports and error messages). If you maintain a project `requirements.txt`, use `pip install -r requirements.txt`.

#### R visualization setup

The pipeline generates R scripts for publication-quality plots. Install R and required packages:

1. **Install R:** [https://www.r-project.org/](https://www.r-project.org/)
2. **R packages** (auto-installed when running scripts): ggplot2, viridis, pheatmap, corrplot, RColorBrewer, reshape2, gridExtra, jsonlite.
3. **Generate and run visualizations:**

```bash
# Python pipeline generates R scripts and data (use your files, directory, or glob)
python FoldKit/crystal_packing_analyzer.py --input model_01.pdb model_02.pdb --compare
# or: --input dir/   or: --input *.pdb

# Run R visualizations
cd comparison_plots/
Rscript master_visualization.R
```

#### Core dependencies

- **BioPython:** structure parsing
- **NumPy/SciPy, pandas, scikit-learn, NetworkX**
- **R** (optional but recommended for full plots)

### Crystal packing scripts

#### `crystal_packing_analyzer.py`

Full pipeline: packing metrics, interface, contact, channel, graph, optional comparative analysis.

```bash
python FoldKit/crystal_packing_analyzer.py --input model_01.pdb
python FoldKit/crystal_packing_analyzer.py --input dir/
python FoldKit/crystal_packing_analyzer.py --input *.pdb --compare --output comparison_results
python FoldKit/crystal_packing_analyzer.py --input *.pdb --sets set_a set_b --output "crystal_analysis_{}"
python FoldKit/crystal_packing_analyzer.py --input *.pdb --sets set_a set_b --output analysis_output --dry-run
```

Options: `--input`, `--output` (use `{}` for one dir per set), `--set`, `--sets`, `--compare`, `--verbose`, `--dry-run`

#### `packing_metrics.py`

Matthews coefficient, solvent content, packing efficiency.

```bash
python FoldKit/packing_metrics.py model_01.pdb
python FoldKit/packing_metrics.py dir/
python FoldKit/packing_metrics.py *.pdb -o results.txt
python FoldKit/packing_metrics.py *.pdb -q -o summary.txt
python FoldKit/packing_metrics.py *.pdb --sets set_a set_b -o "analysis_{}.txt"
python FoldKit/packing_metrics.py *.pdb --per-structure -o "{}_metrics.txt"
python FoldKit/packing_metrics.py *.pdb --per-structure --sets set_a set_b -o "{}_metrics.txt"
```

Options: `-o`, `-q` (quiet), `--per-structure`, `--set`, `--sets`

#### `interface_analyzer.py`

Protein–protein interface analysis (buried surface area, contacts, complementarity).

```bash
python FoldKit/interface_analyzer.py model_01.pdb
python FoldKit/interface_analyzer.py *.pdb -o results.txt
python FoldKit/interface_analyzer.py *.pdb --set set_a,set_b -o by_set.txt
python FoldKit/interface_analyzer.py *.pdb --sets set_a set_b -o "interface_{}.txt"
python FoldKit/interface_analyzer.py *.pdb --per-structure -o "{}_interface.txt"
python FoldKit/interface_analyzer.py *.pdb --per-structure --sets set_a set_b -o "{}_interface.txt"
python FoldKit/interface_analyzer.py *.pdb --sets set_a set_b -o results.txt --dry-run
```

Options: `-o` (use `{}` for per-set or per-structure), `--per-structure`, `--set`, `--sets`, `--dry-run`

#### `interface_molecule_report_csv.py`

Filter an `interface_analyzer.py` text report into CSV (optional filters by PDB basename and chain). Input is typically the file from `interface_analyzer.py … -o results.txt`. Writes one file per structure when using `--output-dir` or `-o` with `{}`.

```bash
python FoldKit/interface_molecule_report_csv.py results.txt -m A -m B --output-dir ./out
python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --group-by-chain --output-dir ./out
python FoldKit/interface_molecule_report_csv.py results.txt --pdbs model_01.pdb,model_02.pdb --chains A,B -o 'out/{}.csv'
python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-regex '^(model\d+[^_]*)_' --output-dir ./out
python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-glob 'model1*' --combine-glob 'model2*' --output-dir ./out
python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B \
  --combine-glob 'model1a*' --combine-glob 'model1del*' --combine-glob 'model1_*' --output-dir ./out
```

Options: `--pdb`, `--pdbs`, `-m`, `--molecule`, `--chains`, `-o`, `--output-dir`, `--group-by-chain`, `--combine-regex` (Python regex; **first capturing group** = merge key / output stem), `--combine-glob` (fnmatch on stem). **Order matters** for globs: use more specific patterns first (e.g. `model1a*` before `model1_*`), because `model1*` matches `model1a_*` and `model1del_*` as well.

The regex above captures `model`, digits, then **any suffix without an underscore** up to the next `_` (e.g. time/lane token `2m`). It maps stems like `model1_2m`, `model1a_7m`, `model1del_12m`, `model2_2m` to `model1.csv`, `model1a.csv`, `model1del.csv`, `model2.csv`. If stems are only `modelN_…` with no letters after the digits, `^(model\d+)_` is enough. To allow only specific suffixes, use `'^(model\d+(?:a|del)?)_'`. To list variants explicitly: `'^(model1a|model1del|model1|model2)_'`. Invalid: `model\d+*` — use `\d+` plus a separate suffix rule (`[^_]*`, `(?:a|del)?`, etc.).

#### `contact_analyzer.py`

Crystal contact and symmetry interaction analysis.

```bash
python FoldKit/contact_analyzer.py model_01.pdb
python FoldKit/contact_analyzer.py *.pdb -o contact_results.txt
python FoldKit/contact_analyzer.py *.pdb --sets set_a set_b -o "contact_{}.txt"
python FoldKit/contact_analyzer.py *.pdb --per-structure -o "{}_contact.txt"
python FoldKit/contact_analyzer.py *.pdb --per-structure --sets set_a set_b -o "{}_contact.txt"
```

Options: same as `interface_analyzer.py`

#### `channel_analyzer.py`

Solvent channel and void space analysis.

```bash
python FoldKit/channel_analyzer.py model_01.pdb
python FoldKit/channel_analyzer.py *.pdb --sets set_a set_b -o "channel_{}.txt"
python FoldKit/channel_analyzer.py *.pdb --per-structure -o "{}_channel.txt"
```

Options: same as `interface_analyzer.py`

#### `graph_analyzer.py`

Graph-theoretical analysis (residue networks, modularity, centrality).

```bash
python FoldKit/graph_analyzer.py model_01.pdb
python FoldKit/graph_analyzer.py *.pdb --sets set_a set_b -o "graph_{}.txt"
python FoldKit/graph_analyzer.py *.pdb --per-structure -o "{}_graph.txt"
```

Options: same as `interface_analyzer.py`

#### `comparative_analyzer.py`

Multi-structure comparison (PCA, clustering, correlation). Requires ≥2 structures per set.

```bash
python FoldKit/comparative_analyzer.py model_01.pdb model_02.pdb
python FoldKit/comparative_analyzer.py *.pdb --sets set_a set_b -o "comparison_{}.txt"
```

Options: `-o`, `--set`, `--sets`, `--dry-run` (no `--per-structure`)

#### `visualization.py`

R script generation for publication-quality plots (used by `crystal_packing_analyzer`).

### Python API usage

Run Python with `FoldKit` on the module path (for example `cd FoldKit` or `PYTHONPATH=FoldKit` from the repo root).

```python
from crystal_packing_analyzer import CrystalPackingAnalyzer

analyzer = CrystalPackingAnalyzer(output_dir="analysis_output")
results = analyzer.analyze_single_structure("model_01.pdb")

model_paths = ["model_01.pdb", "model_02.pdb", "model_03.pdb"]
all_results = [analyzer.analyze_single_structure(p) for p in model_paths]
comparison = analyzer.compare_structures(all_results)
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
from interface_analyzer import InterfaceAnalyzer

analyzer = InterfaceAnalyzer(contact_distance=5.0)
interfaces = analyzer.analyze_interfaces("model_01.pdb")
print(f"Total interfaces: {interfaces['summary']['total_interfaces']}")
```

### Output files

#### Analysis results (typical layout)

```
crystal_analysis_output/
├── model_01_analysis.json
├── model_02_analysis.json
├── comparison_report.txt
├── comparison_plots/
│   ├── summary_statistics.csv
│   ├── correlations.csv
│   ├── pca_results.csv
│   ├── clustering_results.csv
│   ├── master_visualization.R
│   ├── summary_statistics.R
│   ├── correlations.R
│   ├── pca_analysis.R
│   ├── clustering_analysis.R
│   ├── summary_statistics.pdf
│   ├── correlations.pdf
│   ├── pca_analysis.pdf
│   ├── clustering_results.pdf
│   └── cluster_assignments.txt
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

**With Coot:** prepare or superpose structures, export PDBs, then run `crystal_packing_analyzer` on those files.

**With ChimeraX or other viewers:** align if needed, export PDBs, run the same pipeline.

**After FoldKit superposition:**

```bash
python FoldKit/superimpose_coot_SSM.py reference.pdb dir1 dir2 dir3
python FoldKit/crystal_packing_analyzer.py --input SSMaligned2_ref/*.pdb --compare
```

### Scientific applications

- Polymorphism and alternate crystal forms
- Crystal engineering and contact optimization
- Ligand effects on packing (apo vs holo)
- Biological vs crystal interfaces

### Performance considerations

- Single-structure analysis is roughly O(n²) in atom count for contact-heavy steps.
- Typical memory ~100–500 MB per structure; runtime on the order of tens of seconds to a few minutes.

### Quick start

```bash
python FoldKit/quick_start.py
python FoldKit/crystal_packing_analyzer.py --input model_01.pdb
python FoldKit/crystal_packing_analyzer.py --input model_01.pdb model_02.pdb model_03.pdb --compare
cd comparison_plots/
Rscript master_visualization.R
```

#### References (metrics / crystal packing)

1. **Matthews, B.W.** (1968). Solvent content of protein crystals. *J. Mol. Biol.* **33**, 491–497. [DOI: 10.1016/0022-2836(68)90205-2](https://doi.org/10.1016/0022-2836(68)90205-2)
2. **Janin, J.** (1997). Specific versus non-specific contacts in protein crystals. *Nat. Struct. Biol.* **4**, 973–974. [DOI: 10.1038/nsb1297-973](https://doi.org/10.1038/nsb1297-973)
3. **Carugo, O. & Argos, P.** (1997). Protein–protein crystal-packing contacts. *Protein Sci.* **6**, 2261–2263. [DOI: 10.1002/pro.5560061021](https://doi.org/10.1002/pro.5560061021)

*These DOIs were verified against the Crossref REST API.*

#### Troubleshooting (metrics pipeline)

1. **BioPython import errors:** `pip install biopython>=1.79`
2. **Memory:** reduce contact distance cutoffs for large structures
3. **NetworkX:** `pip install networkx>=2.6`
4. **R:** install from [r-project.org](https://www.r-project.org/); on macOS ensure R is on `PATH` (e.g. `export PATH="/usr/local/bin:$PATH"`).

#### Features summary

- Modular analyzers usable alone or via `crystal_packing_analyzer`
- Metrics aligned with common crystallographic practice
- Optional R layer for publication-style figures
- Fits naturally after superposition or alongside external refinement workflows

---

## Appendix

### General

- Scripts accept **PDB and mmCIF** where not otherwise noted.
- **Coot** superposition scripts keep Coot open for interactive review.
- Run **`extract_rmsd.py`** on `coot_log.txt` with **`--format auto`** to classify SSM vs LSQ logs.

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

- **`DALILITE_HOME` / `--dalilite-path`** must point at the DaliLite install root (`bin/dali.pl`).
- **`mkdssp`** must exist at the path in **`bin/mpidali.pm`** (`$DSSP_EXE`). If imports produce empty **`DAT/*.dat`**, install or symlink **`mkdssp`** or edit **`$DSSP_EXE`**. Perl files under **`bin/`** must be readable (fix permissions if you see `Can't locate FSSP.pm: Permission denied`).
- **Path length:** DaliLite Fortran limits full path length (~80 characters); these scripts stage under the system temporary directory automatically.

### Notes

- SSM superposition reports RMSD and alignment statistics in the log.
- LSQ superposition is standard least-squares Cα fitting.
- Trimming aligns residue ranges before LSQ for fair comparison.
- Original input files are not modified by superposition scripts; outputs go to new directories.
- RMSD extraction parses alignment and RMSD lines from Coot logs.

## License

FoldKit is licensed under the [Apache License, Version 2.0](LICENSE).
