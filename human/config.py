"""
Centralized paths for the SETIA human application (BMMC + GTEx).

Replace the hardcoded "/Users/rl884/Downloads/..." literals in each notebook
with imports from this file, e.g. at the top of a notebook:

    import sys; sys.path.append("..")        # so the notebook can see human/config.py
    from config import GENE_LIST_PATH, H5AD_PATH, CHIP_OUT, PPI_OUT, OUT_ROOT

Point DATA_ROOT / OUT_ROOT at your machine via environment variables, or edit
the two defaults below. Everything else is derived relative to them, so no
absolute personal paths remain in the notebooks.
"""
import os
from pathlib import Path

# --- set these two (or export SETIA_HUMAN_DATA / SETIA_HUMAN_OUT) ---
DATA_ROOT = Path(os.environ.get("SETIA_HUMAN_DATA", "./data"))     # raw + downloaded inputs
OUT_ROOT  = Path(os.environ.get("SETIA_HUMAN_OUT",  "./output"))   # everything written by the notebooks

# --- gene / TF lists (small, version-controlled, live in human/data/) ---
GENE_LIST_PATH = DATA_ROOT / "bmmc_gene_list.tsv"     # 56-TF panel + target genes (build_* use this)
SETIA56_PATH   = DATA_ROOT / "BMMC_gene_list.txt"     # Figure3_BMMC master/TF list

# --- large external inputs (NOT in git; see README data table for download) ---
H5AD_PATH      = DATA_ROOT / "Granja2019_annotated.h5ad"   # Granja 2019 BMMC scRNA-seq
GTEX_DIR       = DATA_ROOT / "GTEx_v11"                     # GTEx v11 TPM matrix + sample attributes
LAMBERT_TF     = DATA_ROOT / "Lambert2018_TF_list.xlsx"    # Lambert et al. 2018 Cell, mmc supplement

# --- per-step output dirs (derived; created on demand) ---
SETIA_INPUT_OUT = OUT_ROOT / "setia_input"     # build_setia_input: pseudocells, states, CPM matrix
CPM_MATRIX      = SETIA_INPUT_OUT / "setia_input_linear_CPM.tsv"   # consumed by build_promoter_strength
GENE_LENGTH_OUT = OUT_ROOT / "gene_length"
PROMOTER_OUT    = OUT_ROOT / "promoter_strength"
CHIP_OUT        = OUT_ROOT / "chip_prior"      # build_chip_prior: BMMC TF-DNA prior json
PPI_OUT         = OUT_ROOT / "ppi_prior"       # build_ppi_prior_v2: BMMC TF-TF prior json
LG_OUT          = OUT_ROOT / "lg"              # build_lg (needs CHIP_OUT + PPI_OUT)
GRN_OUT         = OUT_ROOT / "grn"             # SETIA core run output (Steady_state_count.txt, etc.)
CACHE           = OUT_ROOT / "cache"           # download caches (ChIP-Atlas, motifs, genome, ...)

for _d in (SETIA_INPUT_OUT, GENE_LENGTH_OUT, PROMOTER_OUT, CHIP_OUT, PPI_OUT, LG_OUT, GRN_OUT, CACHE):
    _d.mkdir(parents=True, exist_ok=True)
