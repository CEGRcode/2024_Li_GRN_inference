# SETIA human application (BMMC + GTEx)

This directory holds the human arm of SETIA. It contains two independent pieces:

1. **`0_build_inputs/`** builds the SETIA inputs for human bone-marrow
   hematopoiesis (Granja 2019 BMMC), which are then fed to the SETIA engine in
   `../SETIA/` to infer and simulate the GRN.
2. **`1_identity_coherence/`** reproduces the master-regulator identity-coherence
   analysis (Figure 3) for GTEx tissues and for BMMC. These two notebooks are
   self-contained: they call discrete expression states directly from the
   GTEx TPM matrix and the Granja BMMC object and compute within-tissue
   consistency. They do not require the GRN inference run.

Discrete state calling reuses the yeast procedure (the notebooks carry the
state-caller verbatim from the yeast pipeline, so both arms call states the
same way).

## Layout

```
human/
├── config.py               # all paths via HUMAN_BASE (set once, see below)
├── environment.yml          # conda env (setia-human); independent of the yeast env
├── fix_notebooks.py         # utility: re-root absolute paths to HUMAN_BASE and clear outputs
├── 0_build_inputs/          # build the SETIA inputs for BMMC
├── 1_identity_coherence/    # Figure 3: master-regulator within-tissue consistency (standalone)
├── 2_perturbation/          # (placeholder) KO/FE perturbation + regulator taxonomy (Fig 5, 6)
└── data/                    # HUMAN_BASE by default; small lists committed, large files downloaded here
```

## Paths

All notebooks read through `HUMAN_BASE`, defined in `config.py`. By default
`HUMAN_BASE` is `human/data/`; override it by exporting `SETIA_HUMAN_BASE`:

```bash
export SETIA_HUMAN_BASE=/path/to/your/data
```

Lay the directory out as shown in `config.py` and the table below. The small
gene/TF lists (`bmmc_gene_list.tsv`, `BMMC_gene_list.txt`) are committed; the
large external files are downloaded there and are git-ignored.

`fix_notebooks.py` is the helper that produced these notebooks: it rewrites any
`/Users/...`-style absolute path to `f"{HUMAN_BASE}/..."` and clears cell
outputs. Re-run it on any notebook that still carries personal paths
(`python fix_notebooks.py <dir>`; it writes a `.bak` first).

## Run order (0_build_inputs)

1. `build_setia_input.ipynb`  -> pseudocells, discrete states, linear-CPM matrix
2. `build_gene_length.ipynb`  -> per-gene length file (GENCODE)
3. `build_promoter_strength.ipynb` -> CAGE promoter strength (needs the CPM matrix from step 1)
4. `build_chip_prior.ipynb`   -> TF–DNA prior (HOCOMOCO v11 motif scan ∪ ChIP-Atlas Blood)
5. `build_ppi_prior_v2.ipynb` -> TF–TF prior (STRING v12)
6. `build_lg.ipynb`           -> linked-group (complex) file (needs the priors from steps 4, 5)

Then feed the inputs to the SETIA engine in `../SETIA/`. The two
`1_identity_coherence/` notebooks can be run independently at any time.

## External data (download into HUMAN_BASE; not committed)

| Dataset | Version | Used by | Source |
| --- | --- | --- | --- |
| Granja BMMC scRNA-seq (annotated) | Granja et al. 2019 | build_setia_input, Figure3_BMMC | h5ad from the authors / GEO GSE139369 |
| GTEx bulk RNA-seq (TPM + sample attributes) | v11 | Figure3_GTEx | https://www.gtexportal.org/home/downloads |
| Human transcription-factor list | Lambert et al. 2018, Cell | Figure3_GTEx, Figure3_BMMC | Cell 172:650–665, supplement (mmc2.xlsx) |
| GENCODE annotation | release 44 (GRCh38) | build_gene_length, build_chip_prior | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/ |
| hg38 genome FASTA | GRCh38 (UCSC) | build_chip_prior | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/ |
| HOCOMOCO motif models | v11 (HUMAN, mono) | build_chip_prior | https://hocomoco11.autosome.org/downloads_v11 |
| ChIP-Atlas (Blood-class ChIP-seq) | hg38 | build_chip_prior | https://chip-atlas.dbcls.jp/ (experimentList.tab; target tables) |
| STRING protein–protein interactions | v12.0 | build_ppi_prior_v2 | https://string-db.org/ (REST API) |
| FANTOM5 CAGE peaks/TPM | hg38_v8 (reprocessed) | build_promoter_strength | https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/ |

Committed small inputs (in `human/data/`): `bmmc_gene_list.tsv`, `BMMC_gene_list.txt`.

## Dependencies

`conda env create -f environment.yml`. The motif scan in `build_chip_prior` also
uses HOMER as an external command-line tool, and the hg38 FASTA is fetched at
run time. MOODS is installed via pip (`MOODS-python`).
