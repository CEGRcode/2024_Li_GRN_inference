# SETIA human application (BMMC + GTEx)

This directory contains the human arm of SETIA: building the SETIA inputs for
human bone-marrow hematopoiesis (Granja 2019 BMMC) and reproducing the
master-regulator identity-coherence analysis (Figure 3) for both GTEx tissues
and BMMC. The GRN inference and ODE simulation themselves use the SETIA engine
in `../SETIA/`; the notebooks here prepare its inputs and analyze its outputs.

State calling reuses the yeast discrete-state-identification procedure (the
notebooks carry it verbatim from the yeast pipeline, so the two arms call states
the same way).

## Layout

```
human/
├── config.py               # all paths in one place (replaces personal absolute paths)
├── environment.yml         # conda env (setia-human); independent of the yeast env
├── 0_build_inputs/         # build the SETIA inputs for BMMC
├── 1_identity_coherence/   # Figure 3: master-regulator within-tissue consistency
├── 2_perturbation/         # (placeholder) KO/FE perturbation + regulator taxonomy (Fig 5, 6)
└── data/                   # small version-controlled inputs; large inputs downloaded per the table below
```

## Run order

Set `SETIA_HUMAN_DATA` and `SETIA_HUMAN_OUT` (or edit the two defaults in
`config.py`), then run the notebooks in `0_build_inputs/` in this order:

1. `build_setia_input.ipynb` -> pseudocells, discrete states, and the linear-CPM matrix
2. `build_gene_length.ipynb` -> per-gene length file (GENCODE)
3. `build_promoter_strength.ipynb` -> CAGE promoter strength (needs the CPM matrix from step 1)
4. `build_chip_prior.ipynb` -> TF–DNA prior (HOCOMOCO v11 motif scan ∪ ChIP-Atlas Blood)
5. `build_ppi_prior_v2.ipynb` -> TF–TF prior (STRING v12)
6. `build_lg.ipynb` -> linked-group (complex) file (needs the TF–DNA and TF–TF priors from steps 4, 5)

Feed the resulting inputs to the SETIA engine in `../SETIA/`. Then run the two
notebooks in `1_identity_coherence/` (`Figure3_GTEx_from_scratch.ipynb`,
`Figure3_BMMC_from_scratch.ipynb`), which read the steady-state output and
compute master-regulator consistency.

## External data (not in git; download into `human/data/`)

| Dataset | Version | Used by | Source |
| --- | --- | --- | --- |
| Granja BMMC scRNA-seq (annotated) | Granja et al. 2019 | build_setia_input, Figure3_BMMC | h5ad provided by the authors / GEO GSE139369 |
| GTEx bulk RNA-seq (TPM + sample attributes) | v11 | Figure3_GTEx | https://www.gtexportal.org/home/downloads |
| Human transcription-factor list | Lambert et al. 2018, Cell | Figure3_GTEx, Figure3_BMMC | Cell 172:650–665, supplement (mmc) |
| GENCODE annotation | release 44 (GRCh38) | build_gene_length, build_chip_prior | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/ |
| hg38 genome FASTA | GRCh38 (UCSC) | build_chip_prior | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/ |
| HOCOMOCO motif models | v11 (HUMAN, mono) | build_chip_prior | https://hocomoco11.autosome.org/downloads_v11 |
| ChIP-Atlas (Blood-class ChIP-seq) | hg38 | build_chip_prior | https://chip-atlas.dbcls.jp/ (experimentList.tab; target tables) |
| STRING protein–protein interactions | v12.0 | build_ppi_prior_v2 | https://string-db.org/ (REST API) |
| FANTOM5 CAGE peaks/TPM | hg38_v8 (reprocessed) | build_promoter_strength | https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/ |

Small derived inputs that DO live in git: `bmmc_gene_list.tsv`, `BMMC_gene_list.txt`.

## Dependencies

`conda env create -f environment.yml`. Beyond the conda/pip packages, the motif
scan in `build_chip_prior` uses HOMER as an external command-line tool, and the
hg38 FASTA is fetched at run time. MOODS is installed via pip (`MOODS-python`).
