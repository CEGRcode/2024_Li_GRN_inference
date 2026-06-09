"""
Centralized paths for the SETIA human application (BMMC + GTEx).

The notebooks read everything relative to HUMAN_BASE. By default HUMAN_BASE is
the human/data/ directory next to this file; override it with the
SETIA_HUMAN_BASE environment variable to point at a data directory elsewhere.

fix_notebooks.py inserts this import at the top of each notebook:
    import sys; sys.path.append("..")
    from config import HUMAN_BASE

Lay out HUMAN_BASE like this (small lists are committed; large files are
downloaded and git-ignored, see ../.gitignore and the table in README.md):

    <HUMAN_BASE>/
        bmmc_gene_list.tsv                                   (committed)
        BMMC_gene_list.txt                                   (committed)
        GTEx_v11/
            Granja2019_annotated.h5ad                        (download)
            1-s2.0-S0092867418301065-mmc2.xlsx               (Lambert 2018 TF list)
            GTEx_Analysis_2025-08-22_v11_..._gene_tpm.gct    (download)
            GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
"""
import os
from pathlib import Path

HUMAN_BASE = Path(os.environ.get("SETIA_HUMAN_BASE",
                                 Path(__file__).resolve().parent / "data"))

# --- named convenience paths (optional; mirror the layout above) ---
GENE_LIST_PATH   = HUMAN_BASE / "bmmc_gene_list.tsv"
SETIA56_PATH     = HUMAN_BASE / "BMMC_gene_list.txt"
H5AD_PATH        = HUMAN_BASE / "GTEx_v11" / "Granja2019_annotated.h5ad"
LAMBERT_TF       = HUMAN_BASE / "GTEx_v11" / "1-s2.0-S0092867418301065-mmc2.xlsx"
GTEX_TPM         = HUMAN_BASE / "GTEx_v11" / "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct"
GTEX_SAMPLE_ATTR = HUMAN_BASE / "GTEx_v11" / "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
