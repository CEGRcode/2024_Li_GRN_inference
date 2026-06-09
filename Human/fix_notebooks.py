#!/usr/bin/env python3
"""
De-hardcode and clean SETIA human notebooks.

For every *.ipynb in the target directory (default: current dir), this script:
  1. clears all cell outputs and execution counts (clean diffs for git);
  2. rewrites every "/Users/rl884/Downloads/REST" string literal into the
     f-string f"{HUMAN_BASE}/REST", preserving REST exactly (no guessing about
     what each path means, no merging of distinct directories);
  3. inserts a config-import cell at the top if any rewrite was made;
  4. prints a report and flags any other absolute paths it did not touch.

Usage:
    python fix_notebooks.py [DIR]          # edits notebooks in place (writes .bak first)
    python fix_notebooks.py [DIR] --dry    # report only, no writes
"""
import json, re, sys, glob, os, shutil

PERSONAL_PREFIX = "/Users/rl884/Downloads/"
# matches  "/Users/rl884/Downloads/REST"  or  '/Users/rl884/Downloads/REST'
# (also strips a leading f/r so we don't end up with double prefixes)
LIT = re.compile(r"""[fr]?(['"])""" + re.escape(PERSONAL_PREFIX) + r"""([^'"]*)\1""")
OTHER_ABS = re.compile(r"""['"](/(?:Users|home|mnt|scratch)/[^'"]*)['"]""")

IMPORT_CELL = (
    "# --- paths come from human/config.py (auto-inserted by fix_notebooks.py) ---\n"
    "import sys; sys.path.append('..')\n"
    "from config import HUMAN_BASE\n"
)

def repl(m):
    rest = m.group(2)
    return 'f"{HUMAN_BASE}/' + rest + '"'

def has_import(nb):
    for c in nb.get("cells", []):
        if c.get("cell_type") == "code" and "from config import HUMAN_BASE" in "".join(c.get("source", [])):
            return True
    return False

def process(path, dry):
    nb = json.load(open(path))
    n_paths = 0
    leftovers = []
    for c in nb.get("cells", []):
        if c.get("cell_type") != "code":
            continue
        # clear outputs / counts
        c["outputs"] = []
        c["execution_count"] = None
        new_src = []
        for line in c.get("source", []):
            line2, k = LIT.subn(repl, line)
            n_paths += k
            for o in OTHER_ABS.findall(line2):
                leftovers.append(o)
            new_src.append(line2)
        c["source"] = new_src
    if n_paths and not has_import(nb):
        nb["cells"].insert(0, {"cell_type": "code", "metadata": {},
                               "execution_count": None, "outputs": [],
                               "source": [IMPORT_CELL]})
    if not dry:
        shutil.copy(path, path + ".bak")
        json.dump(nb, open(path, "w"), indent=1)
    return n_paths, sorted(set(leftovers))

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry" in sys.argv
    d = args[0] if args else "."
    nbs = sorted(glob.glob(os.path.join(d, "*.ipynb")))
    if not nbs:
        print("no notebooks found in", d); return
    total = 0
    for p in nbs:
        k, left = process(p, dry)
        total += k
        print(f"{os.path.basename(p):44s} rewrote {k:3d} path(s)" + (f"  | clear outputs" ))
        if left:
            print("    !! other absolute paths still present (fix by hand):")
            for o in left: print("       ", o)
    print(f"\n{'DRY RUN, nothing written' if dry else 'done (.bak backups written)'}; total paths rewritten: {total}")

if __name__ == "__main__":
    main()
