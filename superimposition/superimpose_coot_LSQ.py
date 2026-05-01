import hashlib
import os
import sys
from subprocess import Popen, PIPE, STDOUT
import glob
import itertools
import re
import fnmatch

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from superimposition.superimpose_pattern_match import find_ref_model_matches

from utils.cli_log import setup_log_from_argv

_argv_no_log, _ = setup_log_from_argv(
    script_path=__file__,
    argv=sys.argv[1:],
    inputs=[sys.argv[1]] if len(sys.argv) > 1 else [],
    pattern=None,
)
sys.argv = [sys.argv[0]] + _argv_no_log

if _ is not None:
    _.task("Coot LSQ superposition (superimpose_coot_LSQ.py)")
    try:
        _.kv("argv", _argv_no_log)
    except Exception:
        pass

def sanitize_pattern_for_filename(pattern):
    """Make a filter pattern safe for use in log/rmsd filenames."""
    if not pattern:
        return ""
    s = re.sub(r'[\*?\ /\\]', '_', pattern)
    return s[:64]

def expand_output_dir_pattern(template, ref_name=None, filter_pattern=None):
    """
    Expand --output-dir pattern: [reference_name] -> ref stem, [pattern] or [filter] -> sanitised filter.
    Also replaces literal *filter* with the sanitised filter (convenience placeholder).
    """
    if not template:
        return template
    s = template
    ref = ref_name if ref_name is not None else ""
    pat = sanitize_pattern_for_filename(filter_pattern) if filter_pattern else ""
    s = s.replace("[reference_name]", ref).replace("[pattern]", pat).replace("[filter]", pat)
    if pat:
        s = s.replace("*filter*", pat)
    return s


def _announce_superposition_finished(
    log_file, exit_code, summary="", rmsd_format=None
):
    """
    Announce completion on stderr and append the same text to the Coot log.

    While Coot runs, its stdout is often only copied to the log file. Using stderr
    avoids losing the message when Python stdout is fully buffered or redirected
    (non-TTY, some IDEs, wrappers).
    """
    lines = [
        "",
        "=" * 72,
    ]
    if exit_code == 0:
        lines.append(
            "Finished: all superpositions for this run are complete (Coot exit code 0)."
        )
    else:
        lines.append(
            "Coot process ended (exit code {}). Check the log if outputs look wrong.".format(
                exit_code
            )
        )
    if summary:
        lines.append(summary)
    lines.append("Log: {}".format(log_file))
    lines.append("=" * 72)
    if rmsd_format:
        lines.append("")
        lines.append("To extract RMSD values, run:")
        lines.append(
            "python ranking/extract_rmsd.py --format {} {}".format(rmsd_format, log_file)
        )
    text = "\n".join(lines) + "\n"
    print(text, file=sys.stderr, end="", flush=True)
    try:
        with open(log_file, "a", encoding="utf-8", errors="replace") as logf:
            logf.write(text)
    except OSError:
        pass


COOT_STDOUT_ALIGNMENTS_DONE_MARKER = "FOLDKIT_ALIGNMENTS_DONE"


def _pipe_coot_stdout_line(line, log_f, echo_all_stdout=False):
    """Append Coot stdout to the log; echo marker lines to stderr; optionally echo every line to stdout."""
    log_f.write(line)
    log_f.flush()
    stripped = line.rstrip("\n\r")
    if COOT_STDOUT_ALIGNMENTS_DONE_MARKER in line:
        print(stripped, file=sys.stderr, flush=True)
    if echo_all_stdout:
        print(stripped)


def _drain_coot_stdout_and_announce(
    process,
    log_file,
    stays_open_after_run,
    summary_on_exit,
    rmsd_format="lsq",
    echo_all_stdout=False,
):
    """
    Copy Coot stdout into the log. If Coot is left open for inspection, print the same
    completion banner as batch mode as soon as FOLDKIT_ALIGNMENTS_DONE appears, so the user
    does not wait until the window is closed (process.wait would otherwise block indefinitely).
    """
    announced = False
    with open(log_file, "a") as log:
        if process.stdout is None:
            print("Warning: process.stdout is None, cannot capture output.")
        else:
            while True:
                output = process.stdout.readline()
                if output == "" and process.poll() is not None:
                    break
                if output:
                    _pipe_coot_stdout_line(
                        output, log, echo_all_stdout=echo_all_stdout
                    )
                    if (
                        stays_open_after_run
                        and not announced
                        and COOT_STDOUT_ALIGNMENTS_DONE_MARKER in output
                    ):
                        _announce_superposition_finished(
                            log_file,
                            0,
                            summary_on_exit
                            + " Coot is still open for inspection; close the window after inspection.",
                            rmsd_format=rmsd_format,
                        )
                        announced = True
    rc = process.wait()
    if not announced:
        _announce_superposition_finished(
            log_file, rc, summary_on_exit, rmsd_format=rmsd_format
        )
    return rc


# Prepended to every generated Coot LSQ script. Matches Coot User Manual (Least-Squares Fitting) and
# upstream python/coot_lsq.py: clear_lsq_matches + add_lsq_match + apply_lsq_matches.
COOT_LSQ_HELPERS_PY = """
def _foldkit_chain_residue_seq_min_max(imol, chain_id):
    n = chain_n_residues(chain_id, imol)
    if n <= 0:
        return None, None
    lo = None
    hi = None
    for i in range(n):
        seq = seqnum_from_serial_number(imol, chain_id, i)
        if lo is None or seq < lo:
            lo = seq
        if hi is None or seq > hi:
            hi = seq
    return lo, hi

def _foldkit_lsq_match_type_symbol(m):
    if m in ("CA", "ca", "Ca"):
        return 2
    if m in ("main", "Main", "mainchain", "Mainchain"):
        return 1
    if m in ("all", "ALL", "All"):
        return 0
    return 1

def foldkit_least_squares_superpose(imol_ref, imol_mov, ref_chain_id, mov_chain_id, match_type_str):
    ref_lo, ref_hi = _foldkit_chain_residue_seq_min_max(imol_ref, ref_chain_id)
    mov_lo, mov_hi = _foldkit_chain_residue_seq_min_max(imol_mov, mov_chain_id)
    if ref_lo is None or ref_hi is None:
        raise ValueError("empty or invalid reference chain for LSQ: %r" % (ref_chain_id,))
    if mov_lo is None or mov_hi is None:
        raise ValueError("empty or invalid mobile chain for LSQ: %r" % (mov_chain_id,))
    mt = _foldkit_lsq_match_type_symbol(match_type_str)
    clear_lsq_matches()
    add_lsq_match(ref_lo, ref_hi, ref_chain_id, mov_lo, mov_hi, mov_chain_id, mt)
    try:
        return apply_lsq_matches(imol_ref, imol_mov)
    except NameError:
        return apply_lsq_matches_simple(imol_ref, imol_mov)
"""


def _normalise_lsq_match_type(s):
    """Return one of ca, main, all (default main)."""
    if not s:
        return "main"
    x = str(s).strip().lower()
    if x in ("ca", "main", "all"):
        return x
    return "main"


# Legacy alternate template (reload-from-disk display). Same LSQ API as create_lsq_script.
def create_coot_script(
    reference_file,
    model_files,
    output_dir,
    keep_coot_open=True,
    lsq_match_type="main",
):
    """Legacy one-to-many template. If keep_coot_open is False, skip reload-for-display and exit Coot."""
    mt = repr(_normalise_lsq_match_type(lsq_match_type))
    script_content = f"""{COOT_LSQ_HELPERS_PY}
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

reference_mol = read_pdb({reference_file!r})
reference_chain = chain_ids(reference_mol)[0]
ref_name = os.path.splitext(os.path.basename({reference_file!r}))[0]

if not os.path.exists({output_dir!r}):
    os.makedirs({output_dir!r})

for model_path in {model_files!r}:
    model_mol = read_pdb(model_path)
    model_chain = chain_ids(model_mol)[0]
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    print("Aligning " + model_name + " to " + ref_name)
    try:
        foldkit_least_squares_superpose(reference_mol, model_mol, reference_chain, model_chain, {mt})
    except Exception as e:
        print("Error during superposition:", e)
        continue
    output_name = os.path.join({output_dir!r}, model_name + "_LSQaligned2_" + ref_name + ".pdb")
    write_pdb_file(model_mol, output_name)

for i in range(graphics_n_molecules()):
    close_molecule(i)

__POST_CLOSE__
"""
    post_open = f"""
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Reloading structures in Coot for inspection...")
handle_read_draw_molecule_with_recentre({reference_file!r}, 0)
ref_mol = graphics_n_molecules() - 1
set_molecule_bonds_colour_map_rotation(ref_mol, 0)
graphics_to_ca_representation(int(ref_mol))

for model_path in {model_files!r}:
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    aligned_path = os.path.join({output_dir!r}, model_name + "_LSQaligned2_" + ref_name + ".pdb")
    handle_read_draw_molecule_with_recentre(aligned_path, 0)
    mol = graphics_n_molecules() - 1
    set_molecule_bonds_colour_map_rotation(mol, 21 * (mol - ref_mol))
    graphics_to_ca_representation(int(mol))

x, y, z = molecule_centre(ref_mol)
set_rotation_centre(x, y, z)

# coot_real_exit(0)  # Comment out to keep Coot window open
"""
    post_batch = f"""
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Exiting Coot.")
print("Aligned structures saved to: {output_dir}")
coot_real_exit(0)
"""
    post_close = post_open if keep_coot_open else post_batch
    return script_content.replace("__POST_CLOSE__", post_close)

def create_lsq_script(
    reference_file,
    model_files,
    output_dir,
    aligned_tag="_LSQaligned2_",
    keep_coot_open=True,
    ref_chain=None,
    model_chain=None,
    lsq_match_type="main",
):
    """Create Coot script for LSQ superposition (Coot manual: clear_lsq_matches / add_lsq_match / apply_lsq_matches).

    aligned_tag: substring between model basename and reference basename in output PDB names
    (default _LSQaligned2_; use _LSQaligned_ for --pattern mode compatibility).
    """
    mt = _normalise_lsq_match_type(lsq_match_type)
    mt_lit = repr(mt)
    ref_chain_assign = (
        "reference_chain = chain_ids(reference_mol)[0]"
        if ref_chain is None
        else f"reference_chain = {ref_chain!r}"
    )
    model_chain_assign = (
        "model_chain = chain_ids(model_mol)[0]"
        if model_chain is None
        else f"model_chain = {model_chain!r}"
    )
    script_content = f"""{COOT_LSQ_HELPERS_PY}
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

reference_mol = read_pdb({reference_file!r})
{ref_chain_assign}
ref_name = os.path.splitext(os.path.basename({reference_file!r}))[0]

graphics_to_ca_representation(reference_mol)

if not os.path.exists({output_dir!r}):
    os.makedirs({output_dir!r})

loaded_molecules = [reference_mol]

for model_path in {model_files!r}:
    model_mol = read_pdb(model_path)
    {model_chain_assign}
    graphics_to_ca_representation(model_mol)
    loaded_molecules.append(model_mol)
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    print("Aligning " + model_name + " to " + ref_name)
    try:
        foldkit_least_squares_superpose(reference_mol, model_mol, reference_chain, model_chain, {mt_lit})
    except Exception as e:
        print("Error during superposition:", e)
        continue
    output_name = os.path.join({output_dir!r}, model_name + {aligned_tag!r} + ref_name + ".pdb")
    write_pdb_file(model_mol, output_name)

__TAIL__
"""
    tail_open = f"""
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Applying colours in Coot — structures are ready for inspection.")
for i, mol in enumerate(loaded_molecules):
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)

print("Superposition complete. All structures loaded in Coot for visual inspection.")
print("Reference: " + ref_name)
print("Aligned structures saved to: {output_dir}")

# Keep Coot open for manual visualisation
# coot_real_exit(0)  # Commented out to keep Coot window open
"""
    tail_batch = f"""
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Exiting Coot.")
for i in range(graphics_n_molecules()):
    close_molecule(i)
print("Superposition complete (batch). Aligned structures saved to: {output_dir}")
coot_real_exit(0)
"""
    tail = tail_open if keep_coot_open else tail_batch
    return script_content.replace("__TAIL__", tail)

def create_all_vs_all_lsq_script(
    model_files,
    output_dir,
    keep_coot_open=True,
    ref_chain=None,
    model_chain=None,
    lsq_match_type="main",
):
    """Create Coot script for all-vs-all LSQ superposition (Coot manual LSQ API)."""
    mt = repr(_normalise_lsq_match_type(lsq_match_type))
    ref_assign = (
        "reference_chain = chain_ids(reference_mol)[0]"
        if ref_chain is None
        else f"reference_chain = {ref_chain!r}"
    )
    mov_assign = (
        "model_chain = chain_ids(model_mol)[0]"
        if model_chain is None
        else f"model_chain = {model_chain!r}"
    )
    script_content = f"""{COOT_LSQ_HELPERS_PY}
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

if not os.path.exists({output_dir!r}):
    os.makedirs({output_dir!r})

model_files = {model_files!r}
final_molecules = []

for ref_path in model_files:
    reference_mol = read_pdb(ref_path)
    {ref_assign}
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    graphics_to_ca_representation(reference_mol)
    print("Using reference: " + ref_name)

    for model_path in model_files:
        if model_path == ref_path:
            continue
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        print("  Aligning " + model_name + " to " + ref_name)
        model_mol = read_pdb(model_path)
        {mov_assign}
        graphics_to_ca_representation(model_mol)
        try:
            foldkit_least_squares_superpose(reference_mol, model_mol, reference_chain, model_chain, {mt})
        except Exception as e:
            print("Error during superposition of " + model_name + " to " + ref_name + ": " + str(e))
            close_molecule(model_mol)
            continue
        output_name = os.path.join({output_dir!r}, model_name + "_LSQaligned2_" + ref_name + ".pdb")
        write_pdb_file(model_mol, output_name)
        if model_path not in [mol[1] for mol in final_molecules]:
            final_molecules.append((model_mol, model_path))
        else:
            close_molecule(model_mol)

    if ref_path not in [mol[1] for mol in final_molecules]:
        final_molecules.append((reference_mol, ref_path))
    else:
        close_molecule(reference_mol)

__ALL_VS_ALL_TAIL__
"""
    tail_open = f"""
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Reloading structures in Coot for inspection...")
print("\\nReloading structures for visual inspection...")
loaded_for_display = []
unique_files = list(set(model_files))

for i, file_path in enumerate(unique_files):
    mol = read_pdb(file_path)
    graphics_to_ca_representation(mol)
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)
    loaded_for_display.append(mol)
    print("Loaded for display: " + os.path.basename(file_path))

print("\\nAll-vs-all superposition complete. All structures loaded in Coot for visual inspection.")
print("Aligned structures saved to: {output_dir}")

# Keep Coot open for manual visualisation
# coot_real_exit(0)  # Commented out to keep Coot window open
"""
    tail_batch = f"""
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Exiting Coot.")
for i in range(graphics_n_molecules()):
    close_molecule(i)
print("\\nAll-vs-all superposition complete (batch). Aligned structures saved to: {output_dir}")
coot_real_exit(0)
"""
    tail = tail_open if keep_coot_open else tail_batch
    return script_content.replace("__ALL_VS_ALL_TAIL__", tail)


def create_axb_lsq_script(
    ref_files,
    model_files_B,
    output_dirs_per_ref,
    keep_coot_open=False,
    ref_chain=None,
    model_chain=None,
    lsq_match_type="main",
):
    """
    AxB LSQ: for each reference in ref_files, least-squares align each model in model_files_B
    (skipping same path). Chains: --ref-chain/--model-chain when set, else first chain per structure.
    keep_coot_open False (default): call coot_real_exit(0) after all alignments (batch).
    keep_coot_open True (--interactive): reload inputs and leave Coot open.
    """
    exit_line = "# coot_real_exit(0)  # AxB --interactive: keep Coot open after reload"
    if not keep_coot_open:
        exit_line = "coot_real_exit(0)  # AxB default: non-interactive batch (exit when done)"

    mt = repr(_normalise_lsq_match_type(lsq_match_type))
    ref_assign = (
        "reference_chain = chain_ids(reference_mol)[0]"
        if ref_chain is None
        else f"reference_chain = {ref_chain!r}"
    )
    mov_assign = (
        "model_chain = chain_ids(model_mol)[0]"
        if model_chain is None
        else f"model_chain = {model_chain!r}"
    )
    ref_configs = list(zip(ref_files, output_dirs_per_ref))
    script_content = f"""{COOT_LSQ_HELPERS_PY}
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

ref_configs = {repr(ref_configs)}
model_paths = {repr(model_files_B)}
keep_coot_open = {repr(bool(keep_coot_open))}

for ref_path, out_dir in ref_configs:
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    reference_mol = read_pdb(ref_path)
    {ref_assign}
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    graphics_to_ca_representation(reference_mol)
    print("LSQ AxB: reference " + ref_name)
    for model_path in model_paths:
        if model_path == ref_path:
            continue
        model_mol = read_pdb(model_path)
        {mov_assign}
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        print("  Aligning " + model_name + " to " + ref_name)
        graphics_to_ca_representation(model_mol)
        try:
            foldkit_least_squares_superpose(reference_mol, model_mol, reference_chain, model_chain, {mt})
        except Exception as e:
            print("Error LSQ " + model_name + " to " + ref_name + ": " + str(e))
            close_molecule(model_mol)
            continue
        output_name = os.path.join(out_dir, model_name + "_LSQaligned2_" + ref_name + ".pdb")
        write_pdb_file(model_mol, output_name)
        close_molecule(model_mol)
    close_molecule(reference_mol)

for i in range(graphics_n_molecules()):
    close_molecule(i)

if keep_coot_open:
    print("FOLDKIT_ALIGNMENTS_DONE: All AxB superpositions written to disk. Reloading structures in Coot for inspection...")
    print("\\nReloading structures for visual inspection (LSQ AxB)...")
    unique_paths = []
    seen = set()
    for p, _ in ref_configs:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            unique_paths.append(p)
    for p in model_paths:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            unique_paths.append(p)
    for i, file_path in enumerate(unique_paths):
        mol = read_pdb(file_path)
        graphics_to_ca_representation(mol)
        set_molecule_bonds_colour_map_rotation(mol, 20 * i)
        print("Loaded for display: " + os.path.basename(file_path))
    print("LSQ AxB complete.")
else:
    print("FOLDKIT_ALIGNMENTS_DONE: All AxB superpositions written to disk. Exiting Coot.")

__EXIT_LINE__
""".replace("__EXIT_LINE__", exit_line)

    return script_content


def run_pattern_lsq_superposition(
    ref_dir,
    model_dir,
    ref_pattern,
    model_pattern,
    target_pattern=None,
    divider=None,
    strict_position=False,
    ref_file_pattern="*.pdb",
    model_file_pattern="*.cif",
    output_suffix="_LSQaligned_",
    keep_coot_open=True,
    lsq_match_type="main",
):
    """Pair reference and model files by patterns, then LSQ-superpose each match in Coot (legacy trim_superimposeLSQ_pattern behaviour)."""
    matches = find_ref_model_matches(
        ref_dir,
        model_dir,
        ref_pattern,
        model_pattern,
        target_pattern,
        divider,
        strict_position,
        ref_file_pattern,
        model_file_pattern,
    )
    if not matches:
        print("No matching files found to process.")
        return

    for ref_file, subdir, model_files in matches:
        subdir_name = os.path.basename(subdir)
        ref_name = os.path.splitext(os.path.basename(ref_file))[0]
        output_dir = f"{subdir_name}{output_suffix}{ref_name}"

        if os.path.exists(output_dir):
            print(f"Directory {output_dir} already exists. Skipping this set.")
            continue

        os.makedirs(output_dir)
        print(f"\nProcessing {subdir_name}...")
        print(f"Reference: {os.path.basename(ref_file)}")
        print(f"Models to align: {len(model_files)}")

        script_file = "temp_coot_script.py"
        with open(script_file, "w") as f:
            f.write(
                create_lsq_script(
                    ref_file,
                    model_files,
                    output_dir,
                    aligned_tag="_LSQaligned_",
                    keep_coot_open=keep_coot_open,
                    lsq_match_type=lsq_match_type,
                )
            )

        try:
            log_file = os.path.join(output_dir, "coot_log.txt")
            with open(log_file, "w") as log:
                log.write("# LSQ alignment log\n")
                log.write(
                    "# Mode: pattern (pair reference and model files by name patterns; same pairing rules as SSM --pattern)\n"
                )
                log.write(f"# Reference: {ref_file}\n")
                log.write(f"# Models directory: {subdir}\n")
                log.write(f"# Number of models: {len(model_files)}\n")
                log.write(f"# LSQ match type: {lsq_match_type}\n\n")

            process = Popen(
                ["coot", "--script", script_file],
                stdout=PIPE,
                stderr=STDOUT,
                universal_newlines=True,
                bufsize=1,
            )
            _drain_coot_stdout_and_announce(
                process,
                log_file,
                keep_coot_open,
                "LSQ pattern mode — subdirectory: {}.".format(subdir_name),
                rmsd_format="lsq",
            )
        except Exception as e:
            print(f"Error during superposition: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)


def main_pattern_mode():
    """CLI: --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"""
    args = sys.argv[2:]
    strict_position = False
    divider = None
    target_pattern = None
    ref_file_pattern = "*.pdb"
    model_file_pattern = "*.cif"
    output_suffix = "_LSQaligned_"
    keep_coot_open = True
    lsq_match_type = "main"

    i = 0
    while i < len(args):
        if args[i] == "--strict-position":
            strict_position = True
            args.pop(i)
        elif args[i] == "--not-interactive":
            keep_coot_open = False
            args.pop(i)
        elif args[i].startswith("--divider="):
            divider = args[i].split("=", 1)[1]
            args.pop(i)
        elif args[i].startswith("--ref-file-pattern="):
            ref_file_pattern = args[i].split("=", 1)[1]
            args.pop(i)
        elif args[i].startswith("--model-file-pattern="):
            model_file_pattern = args[i].split("=", 1)[1]
            args.pop(i)
        elif args[i].startswith("--output-suffix="):
            output_suffix = args[i].split("=", 1)[1]
            args.pop(i)
        elif args[i].startswith("--lsq-match-type="):
            lsq_match_type = _normalise_lsq_match_type(args[i].split("=", 1)[1])
            args.pop(i)
        else:
            i += 1

    if len(args) < 4:
        print("Pattern mode: pair reference and model files by name patterns (superimpose_pattern_match).")
        print("This mode cannot be combined with --all-vs-all, --reference, --filter, AxB filters, etc.")
        print("")
        print(
            "Usage: python superimposition/superimpose_coot_LSQ.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"
        )
        print("\nOptions:")
        print("  --strict-position        ref_pattern before divider and target_pattern after")
        print("  --divider=STRING")
        print("  --ref-file-pattern=GLOB  (default: \"*.pdb\")")
        print("  --model-file-pattern=GLOB (default: \"*.cif\")")
        print("  --output-suffix=STRING   output dir suffix (default: \"_LSQaligned_\")")
        print("  --not-interactive        exit Coot after each alignment set (default: keep open)")
        print("  --lsq-match-type=TYPE    ca | main | all (default: main)")
        print("\nExamples:")
        print("  python superimposition/superimpose_coot_LSQ.py --pattern /path/to/refs /path/to/models ref_id model_id")
        print(
            "  python superimposition/superimpose_coot_LSQ.py --pattern --ref-file-pattern=\"*final.pdb\" --model-file-pattern=\"*.pdb\" "
            "/path/to/refs /path/to/models ref_id model_id"
        )
        sys.exit(1)

    ref_dir = os.path.abspath(args[0])
    model_dir = os.path.abspath(args[1])
    ref_pattern = args[2]
    model_pattern = args[3]
    if len(args) >= 5:
        target_pattern = args[4]

    if not os.path.exists(ref_dir):
        print(f"Error: Reference directory '{ref_dir}' does not exist.")
        sys.exit(1)
    if not os.path.exists(model_dir):
        print(f"Error: Model directory '{model_dir}' does not exist.")
        sys.exit(1)

    run_pattern_lsq_superposition(
        ref_dir,
        model_dir,
        ref_pattern,
        model_pattern,
        target_pattern,
        divider,
        strict_position,
        ref_file_pattern,
        model_file_pattern,
        output_suffix,
        keep_coot_open=keep_coot_open,
        lsq_match_type=lsq_match_type,
    )


def _pattern_has_glob_chars(pattern: str) -> bool:
    return any(ch in pattern for ch in ["*", "?", "["])

def _matches_filter(basename: str, pattern: str) -> bool:
    """
    Match `pattern` against a file basename.
    - If pattern contains glob chars (* ? [), use shell-style wildcard matching.
    - Otherwise, preserve legacy behaviour: substring match.
    Also matches against the filename stem (without extension) for convenience.
    """
    stem, _ = os.path.splitext(basename)
    if _pattern_has_glob_chars(pattern):
        return (
            fnmatch.fnmatchcase(basename, pattern)
            or fnmatch.fnmatchcase(stem, pattern)
        )
    return (pattern in basename) or (pattern in stem)

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--pattern":
        main_pattern_mode()
        return

    # Check for minimum arguments
    if len(sys.argv) < 3:
        print("Usage: python superimposition/superimpose_coot_LSQ.py [--output-dir=DIR] --all-vs-all --filter=TAG dir1 [dir2 ...]")
        print("       python superimposition/superimpose_coot_LSQ.py [--output-dir=DIR] --reference=REF_FILE --filter=TAG dir1 [dir2 ...]")
        print("       python superimposition/superimpose_coot_LSQ.py [--output-dir=DIR] dir1 [dir2 ...]")
        print(
            "       python superimposition/superimpose_coot_LSQ.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"
        )
        print("\nAdditional options:")
        print("  --filter=TAG         Substring or glob on basename / stem (single-set all-vs-all)")
        print("  --ref-filter=TAG     Like --filter, but selects reference set (A) for AxB mode")
        print("  --model-filter=TAG   Like --filter, but selects model set (B) for AxB mode")
        print("  --reference=FILE     Reference structure (one-to-many); alternative to interactive pick")
        print("  --ref=FILE           Same as --reference=")
        print("  --output-dir=DIR     Output dir; placeholders [reference_name], [filter], *filter*")
        print("  --out-dir=DIR        Alias for --output-dir=")
        print("  --ref-chain=CHAIN    Explicit reference chain (one-to-many, all-vs-all, AxB)")
        print("  --model-chain=CHAIN  Explicit model chain (same)")
        print("  --lsq-match-type=TYPE  LSQ only: ca | main | all (default: main; Coot manual match types)")
        print("  --interactive        AxB only: after all alignments, reload structures and keep Coot open")
        print("  --not-interactive    One-to-many and single-set all-vs-all: exit Coot when done (default: keep open)")
        print("\nExamples:")
        print("  python superimposition/superimpose_coot_LSQ.py --reference=/path/to/ref.cif --filter=set_a /path/to/models/")
        print("  python superimposition/superimpose_coot_LSQ.py --all-vs-all --filter=set_a /path/to/models/")
        print("  python superimposition/superimpose_coot_LSQ.py --all-vs-all --ref-filter=set_a --model-filter=set_b models_dir1 models_dir2")
        print("  python superimposition/superimpose_coot_LSQ.py --pattern /path/to/refs /path/to/models ref_id model_id")
        print("")
        print("Pattern mode (--pattern) is a separate entry point, not combinable with the lines above:")
        print("  --pattern must be the FIRST argument after the script name.")
        print("  Do not use --pattern together with --all-vs-all, --reference/--ref, --filter,")
        print("  --ref-filter, --model-filter, or interactive directory-only pick in the same command.")
        sys.exit(1)
    
    # Parse arguments
    all_vs_all = False
    filter_pattern = None
    ref_filter = None
    model_filter = None
    reference_file = None
    output_dir_arg = None
    directories = []
    ref_chain = None
    model_chain = None
    axb_keep_coot_open = False
    legacy_keep_coot_open = True
    lsq_match_type = "main"

    # Allow --output-dir anywhere (e.g. after positionals); only this pass sets it
    for a in sys.argv[1:]:
        if a.startswith("--output-dir=") or a.startswith("--out-dir="):
            output_dir_arg = a.split("=", 1)[1]
    
    for arg in sys.argv[1:]:
        if arg == "--all-vs-all":
            all_vs_all = True
        elif arg == "--interactive":
            axb_keep_coot_open = True
        elif arg == "--not-interactive":
            legacy_keep_coot_open = False
        elif arg.startswith("--filter="):
            filter_pattern = arg.split("=", 1)[1]
        elif arg.startswith("--ref-filter="):
            ref_filter = arg.split("=", 1)[1]
        elif arg.startswith("--model-filter="):
            model_filter = arg.split("=", 1)[1]
        elif arg.startswith("--reference=") or arg.startswith("--ref="):
            reference_file = arg.split("=", 1)[1]
            reference_file = os.path.abspath(reference_file)
            if not os.path.exists(reference_file):
                print("Error: Reference file '{}' does not exist.".format(reference_file))
                sys.exit(1)
        elif arg.startswith("--ref-chain="):
            ref_chain = arg.split("=", 1)[1]
        elif arg.startswith("--model-chain="):
            model_chain = arg.split("=", 1)[1]
        elif arg.startswith("--lsq-match-type="):
            lsq_match_type = _normalise_lsq_match_type(arg.split("=", 1)[1])
        elif arg.startswith("--output-dir=") or arg.startswith("--out-dir="):
            pass  # already collected above
        elif not arg.startswith("--"):
            directories.append(os.path.abspath(arg))
    
    if not directories:
        print("Error: No directories specified.")
        sys.exit(1)
    
    # Check for conflicting options
    if all_vs_all and reference_file:
        print(
            "Error: Cannot use both --all-vs-all and --reference/--ref together.\n"
            "  One-to-many: use --reference=FILE without --all-vs-all.\n"
            "  AxB (two sets): use --all-vs-all --ref-filter=... --model-filter=... dir_A dir_B; "
            "each structure in set A is a reference in turn (no single --reference)."
        )
        sys.exit(1)

    explicit_chains = ref_chain is not None or model_chain is not None
    if explicit_chains:
        if ref_chain is None:
            ref_chain = "A"
        if model_chain is None:
            model_chain = "A"
    
    # Collect model files from all specified directories
    model_files = []
    for directory in directories:
        if not os.path.exists(directory):
            print("Warning: Directory '{}' does not exist. Skipping.".format(directory))
            continue
        
        # Get all PDB and CIF files from the directory and subdirectories
        dir_models = []
        # Search directly in the directory
        dir_models.extend(glob.glob(os.path.join(directory, "*.pdb")))
        dir_models.extend(glob.glob(os.path.join(directory, "*.cif")))
        
        # Search recursively in subdirectories
        dir_models.extend(glob.glob(os.path.join(directory, "**", "*.pdb"), recursive=True))
        dir_models.extend(glob.glob(os.path.join(directory, "**", "*.cif"), recursive=True))
        
        # Remove duplicates that might occur from searching both directly and recursively
        dir_models = list(set(dir_models))
        
        if not dir_models:
            print("Warning: No PDB or CIF files found in '{}' or its subdirectories. Skipping.".format(directory))
            continue
        
        model_files.extend(dir_models)
    
    if not model_files:
        print("Error: No PDB or CIF files found in any of the specified directories.")
        sys.exit(1)

    # Two-set AxB mode: --ref-filter / --model-filter (requires --all-vs-all)
    use_two_sets = bool(ref_filter or model_filter)
    if use_two_sets and not all_vs_all:
        print("Error: --ref-filter/--model-filter currently require --all-vs-all (AxB mode).")
        sys.exit(1)

    if use_two_sets and filter_pattern:
        print("Warning: --filter is ignored when --ref-filter/--model-filter are used.")
        filter_pattern = None

    ref_files = None
    model_files_B = None

    if use_two_sets:
        # Build reference set A
        if ref_filter:
            ref_files = [
                f for f in model_files
                if _matches_filter(os.path.basename(f), ref_filter)
            ]
        else:
            ref_files = list(model_files)

        # Build model set B
        if model_filter:
            model_files_B = [
                f for f in model_files
                if _matches_filter(os.path.basename(f), model_filter)
            ]
        else:
            model_files_B = list(model_files)

        if not ref_files:
            if ref_filter is not None:
                print(f"Error: No structures match --ref-filter='{ref_filter}'.")
            else:
                print(
                    "Error: Reference set (A) is empty (no structure files available; "
                    "check input directories)."
                )
            sys.exit(1)
        if not model_files_B:
            if model_filter is not None:
                print(f"Error: No structures match --model-filter='{model_filter}'.")
            else:
                print(
                    "Error: Model set (B) is empty (no structure files available; "
                    "check input directories)."
                )
            sys.exit(1)

        print(f"AxB mode: {len(ref_files)} references (A), {len(model_files_B)} models (B).")

    # Single-set filter (legacy --filter)
    if not use_two_sets and filter_pattern:
        filtered_models = []
        for model in model_files:
            if _matches_filter(os.path.basename(model), filter_pattern):
                filtered_models.append(model)
        
        if not filtered_models:
            print("Error: No models match the filter pattern '{}'.".format(filter_pattern))
            print("  (Matched against filename and filename without extension.)")
            print("  Sample of {} candidate basenames:".format(min(10, len(model_files))))
            for m in model_files[:10]:
                print("    {}".format(os.path.basename(m)))
            if len(model_files) > 10:
                print("    ... and {} more".format(len(model_files) - 10))
            sys.exit(1)
        
        model_files = filtered_models

    # Convenience printer for structure lists (for user feedback)
    def _print_structures(label, files):
        print(f"{label} ({len(files)} structures):")
        for f in files:
            print(f"  - {os.path.basename(f)}")

    if use_two_sets:
        # AxB mode: one Coot script for all ref x model pairs (batch exits Coot unless --interactive)
        _print_structures("Reference set (A)", ref_files)
        _print_structures("Model set (B)", model_files_B)

        ref_stems = [os.path.splitext(os.path.basename(p))[0] for p in ref_files]
        stem_counts = {}
        for s in ref_stems:
            stem_counts[s] = stem_counts.get(s, 0) + 1

        output_dirs_per_ref = []
        for ref_path, ref_stem in zip(ref_files, ref_stems):
            if stem_counts[ref_stem] > 1:
                ref_name = "{}__{}".format(
                    ref_stem,
                    hashlib.sha256(os.path.abspath(ref_path).encode("utf-8")).hexdigest()[:8],
                )
            else:
                ref_name = ref_stem
            if output_dir_arg:
                od = expand_output_dir_pattern(
                    output_dir_arg, ref_name=ref_name, filter_pattern=None
                )
            else:
                od = f"LSQaligned2_{ref_name}"
            output_dirs_per_ref.append(od)

        if any(stem_counts[s] > 1 for s in stem_counts):
            print(
                "Note: Duplicate reference basenames: each output dir uses the basename plus "
                "__<8 hex chars> from the reference file path so directories stay distinct."
            )

        existing = [d for d in output_dirs_per_ref if os.path.exists(d)]
        if existing:
            response = input(
                f"{len(existing)} output director(y/ies) already exist. "
                "Files may be overwritten. Continue? (y/n): "
            )
            if response.lower() != 'y':
                print("Operation cancelled.")
                sys.exit(0)

        for od in output_dirs_per_ref:
            if not os.path.exists(od):
                os.makedirs(od)

        tag_parts = ["AxB"]
        if ref_filter:
            tag_parts.append(f"ref_{sanitize_pattern_for_filename(ref_filter)}")
        if model_filter:
            tag_parts.append(f"model_{sanitize_pattern_for_filename(model_filter)}")
        log_suffix = "_".join(tag_parts)
        log_basename = f"coot_log_{log_suffix}.txt"
        log_file = log_basename
        print(
            "Mode: LSQ AxB (non-interactive batch, Coot exits when done)"
            if not axb_keep_coot_open
            else "Mode: LSQ AxB (interactive: Coot stays open after reload)"
        )
        print(f"Creating log file at: {log_file}")

        script_file = "temp_coot_script.py"
        if os.path.exists(script_file):
            print("Warning: Temporary file '{}' exists and will be overwritten".format(script_file))

        script_body = create_axb_lsq_script(
            ref_files,
            model_files_B,
            output_dirs_per_ref,
            keep_coot_open=axb_keep_coot_open,
            ref_chain=ref_chain if explicit_chains else None,
            model_chain=model_chain if explicit_chains else None,
            lsq_match_type=lsq_match_type,
        )

        with open(script_file, "w") as f:
            f.write(script_body)

        try:
            print("Starting Coot process (LSQ AxB)...")
            with open(log_file, "w") as log:
                log.write("# LSQ AxB alignment log\n")
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if ref_filter:
                    log.write(f"# Ref filter: {ref_filter}\n")
                if model_filter:
                    log.write(f"# Model filter: {model_filter}\n")
                log.write(
                    f"# AxB keep_coot_open ( --interactive ): {axb_keep_coot_open}\n"
                )
                log.write(f"# LSQ match type: {lsq_match_type}\n")
                if explicit_chains:
                    log.write(f"# Ref chain: {ref_chain}  Model chain: {model_chain}\n")
                log.write(f"# References (A): {len(ref_files)}\n")
                log.write(f"# Models (B): {len(model_files_B)}\n\n")

            process = Popen(
                ["coot", "--script", script_file],
                stdout=PIPE,
                stderr=STDOUT,
                universal_newlines=True,
                bufsize=1,
            )
            _drain_coot_stdout_and_announce(
                process,
                log_file,
                axb_keep_coot_open,
                "LSQ AxB: {} reference(s), {} model(s) in set B.".format(
                    len(ref_files), len(model_files_B)
                ),
                rmsd_format="lsq",
            )
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)

        return
    
    print("Found {} structure files matching criteria:".format(len(model_files)))
    for f in model_files:
        print("  - {}".format(os.path.basename(f)))
    
    if all_vs_all:
        # All-vs-all superposition
        print("\nPerforming all-vs-all superposition...")
        
        # Output directory: user-specified (patterns [reference_name], [pattern]/[filter]) or default
        if output_dir_arg:
            output_dir = expand_output_dir_pattern(output_dir_arg, ref_name=None, filter_pattern=filter_pattern)
        elif filter_pattern:
            output_dir = "LSQaligned_all_vs_all_{}".format(sanitize_pattern_for_filename(filter_pattern) or filter_pattern)
        else:
            output_dir = "LSQaligned_all_vs_all"
        
        if os.path.exists(output_dir):
            response = input("Directory '{}' already exists. Files may be overwritten. Continue? (y/n): ".format(output_dir))
            if response.lower() != 'y':
                print("Operation cancelled.")
                sys.exit(0)
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Match SSM naming logic: coot_log.txt or coot_log_<sanitized_filter>.txt
        log_suffix = sanitize_pattern_for_filename(filter_pattern) if filter_pattern else ""
        log_basename = f"coot_log_{log_suffix}.txt" if log_suffix else "coot_log.txt"
        log_file = os.path.join(output_dir, log_basename)
        print("Creating log file at: {}".format(log_file))
        
        # Create temporary script file
        script_file = "temp_coot_script.py"
        if os.path.exists(script_file):
            print("Warning: Temporary file '{}' exists and will be overwritten".format(script_file))
        
        with open(script_file, "w") as f:
            f.write(
                create_all_vs_all_lsq_script(
                    model_files,
                    output_dir,
                    keep_coot_open=legacy_keep_coot_open,
                    ref_chain=ref_chain if explicit_chains else None,
                    model_chain=model_chain if explicit_chains else None,
                    lsq_match_type=lsq_match_type,
                )
            )
        
        # Run Coot with the script and capture output
        try:
            print("Starting Coot process...")
            
            # Initialise main log file with header
            with open(log_file, 'w') as log:
                log.write("# LSQ all-vs-all alignment log\n")
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
                log.write(
                    "# Keep Coot open after run: {}\n".format(legacy_keep_coot_open)
                )
                log.write("# LSQ match type: {}\n".format(lsq_match_type))
                if explicit_chains:
                    log.write(
                        "# Ref chain: {}  Model chain: {}\n".format(ref_chain, model_chain)
                    )
                log.write("# Number of models: {}\n\n".format(len(model_files)))
            
            # Run process and capture output
            process = Popen(["coot", "--script", script_file], 
                           stdout=PIPE, stderr=STDOUT,
                           universal_newlines=True, bufsize=1)
            
            num_alignments = len(model_files) * (len(model_files) - 1)
            _drain_coot_stdout_and_announce(
                process,
                log_file,
                legacy_keep_coot_open,
                "LSQ all-vs-all: {} structures, {} directed pairwise alignments. Output directory: {}.".format(
                    len(model_files),
                    num_alignments,
                    output_dir,
                ),
                rmsd_format="lsq",
            )
                
        finally:
            # Clean up temporary script
            if os.path.exists(script_file):
                os.remove(script_file)
        
    else:
        # One-to-many superposition
        if reference_file:
            # Reference file specified from command line
            print("Using reference file: {}".format(os.path.basename(reference_file)))
            
            # Remove the reference file from the list of models to align when present
            model_files = [f for f in model_files if os.path.abspath(f) != reference_file]
            
            if not model_files:
                print("Error: No model files remain after removing reference file.")
                sys.exit(1)
        else:
            # Original behaviour - ask for reference file
            print("\nPlease select a reference file from the list above (enter the number):")
            for i, f in enumerate(model_files):
                print("  {}. {}".format(i+1, os.path.basename(f)))
            
            try:
                choice = int(input("Enter number: "))
                if choice < 1 or choice > len(model_files):
                    print("Invalid selection.")
                    sys.exit(1)
                
                reference_file = model_files[choice-1]
                print("Selected reference: {}".format(os.path.basename(reference_file)))
                
                # Remove reference file from the list of models to align
                model_files = [f for f in model_files if f != reference_file]
            except (ValueError, KeyboardInterrupt):
                print("Invalid input or operation cancelled.")
                sys.exit(1)
        
        # Output directory: user-specified (patterns [reference_name], [pattern]/[filter]) or default
        ref_name = os.path.splitext(os.path.basename(reference_file))[0]
        if output_dir_arg:
            output_dir = expand_output_dir_pattern(output_dir_arg, ref_name=ref_name, filter_pattern=filter_pattern)
        else:
            output_dir = "LSQaligned2_{}".format(ref_name)
        
        if os.path.exists(output_dir):
            response = input("Directory '{}' already exists. Files may be overwritten. Continue? (y/n): ".format(output_dir))
            if response.lower() != 'y':
                print("Operation cancelled.")
                sys.exit(0)
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Match SSM naming logic: coot_log.txt or coot_log_<sanitized_filter>.txt
        log_suffix = sanitize_pattern_for_filename(filter_pattern) if filter_pattern else ""
        log_basename = f"coot_log_{log_suffix}.txt" if log_suffix else "coot_log.txt"
        log_file = os.path.join(output_dir, log_basename)
        print("Creating log file at: {}".format(log_file))
        
        # Create temporary script file
        script_file = "temp_coot_script.py"
        if os.path.exists(script_file):
            print("Warning: Temporary file '{}' exists and will be overwritten".format(script_file))
        
        with open(script_file, "w") as f:
            f.write(
                create_lsq_script(
                    reference_file,
                    model_files,
                    output_dir,
                    keep_coot_open=legacy_keep_coot_open,
                    ref_chain=ref_chain if explicit_chains else None,
                    model_chain=model_chain if explicit_chains else None,
                    lsq_match_type=lsq_match_type,
                )
            )
        
        # Run Coot with the script and capture output
        try:
            print("Starting Coot process...")
            
            # Initialise main log file with header
            with open(log_file, 'w') as log:
                log.write("# LSQ alignment log\n")
                log.write("# Reference: {}\n".format(reference_file))
                log.write("# Directories: {}\n".format(", ".join(directories)))
                log.write("# Number of models: {}\n".format(len(model_files)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
                log.write(
                    "# Keep Coot open after run: {}\n".format(legacy_keep_coot_open)
                )
                log.write("# LSQ match type: {}\n".format(lsq_match_type))
                if explicit_chains:
                    log.write(
                        "# Ref chain: {}  Model chain: {}\n".format(ref_chain, model_chain)
                    )

            # Run process and capture output
            process = Popen(["coot", "--script", script_file], 
                           stdout=PIPE, stderr=STDOUT,
                           universal_newlines=True, bufsize=1)
            
            _drain_coot_stdout_and_announce(
                process,
                log_file,
                legacy_keep_coot_open,
                "LSQ one-to-many: output directory {}.".format(output_dir),
                rmsd_format="lsq",
                echo_all_stdout=True,
            )

        except Exception as e:
            print("Error running Coot: {}".format(e))
            
        except KeyboardInterrupt:
            print("\nProcess interrupted.")
            
        finally:
            # Clean up temporary script
            if os.path.exists(script_file):
                os.remove(script_file)

if __name__ == "__main__":
    main() 