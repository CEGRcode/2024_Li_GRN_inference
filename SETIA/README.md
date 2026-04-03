# Gene Regulatory Network Inference Using SETIA

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA

## Overview
SETIA integrates multiple layers of molecular evidence to infer mechanistically interpretable and dynamically executable gene regulatory networks (GRNs).

Supported evidence layers include:

- RNA-seq transcriptional profiles across genotypes and perturbations
- ChIP-exo TF–DNA binding
- PRO-seq / CAGE-seq promoter activity
- Protein–protein interactions (PPI), including protein–protein colocalization (PPC)

These data are integrated to construct and constrain GRNs whose dynamics are modeled using ordinary differential equations (ODEs).

---

## Key Features
- Infers dynamical GRNs that reproduce stable transcriptional states
- Integrates multi-omics data as mechanistic structural constraints
- Supports combinatorial regulation (AND / OR logic)
- Implements parallelized optimization for scalable GRN inference
- Provides an interactive web platform for simulation and visualization

---

## Pipeline to Reproduce GRNs from the Paper
1. `0_Download_and_touch.sh`  
   Downloads and initializes required data

2. `1_Preprocessing_dicrete_gene_expression_state_identification.sh`  
   Identifies discrete gene expression states

3. `2_GRN_inference.sh`  
   Performs GRN inference

4. `3_Final_GRN_construction.sh`  
   Filters and constructs the final GRN

5. `4_Scale_GRN_to_nonTF_targets.sh`  
   Extends the GRN to downstream non-TF targets

---

## Dependencies
We recommend initializing the following dedicated Anaconda environment for reproducibility:

```bash
conda create -n EvoAlg -c conda-forge -y \
    python=3.10 \
    numpy=1.26 \
    pandas=2.2 \
    scipy=1.12 \
    sympy=1.12 \
    numexpr=2.10 \
    h5py=3.11 \
    matplotlib=3.8 \
    scikit-learn=1.4 \
    pip
conda activate EvoAlg
pip install flask tqdm seaborn joblib
```

---

## Outputs
The inferred GRN is written to the `./result` directory.

Key outputs include:
- `GRN_filtered_Sc_GRN_final_2.json`: final GRN structure
- `GRN_filtered_full_Sc_GRN_final.txt`: fitted model parameters
- `*.pkl`: intermediate and final serialized GRN instances
- `sys_cache_*`: cached searched or simulated GRN configurations

---

## Minimal Example
This example demonstrates how to run SETIA on a small test dataset. Because the code is optimized for HPC environments, local execution may take longer (~10 minutes).

The inferred GRN will be written to:
`./result/GRN_filtered_Sc_GRN_final_raw.json`

This output can be uploaded to the GRN simulator website for interactive exploration:
https://grn.cac.cornell.edu:5000/

```bash
cd SETIA
conda activate EvoAlg

python EGRN_Multi_Genalg_Combinatorial_2025.py \
    -r data/example_RNA_profiles.txt \
    -n 3000 \
    -i 10 \
    -p 0 \
    -t data/example_promoter_strength.txt \
    -l data/example_gene_length.txt \
    -o "Minimal_example" \
    -e 42 \
    -f 0 \
    -k 0

python Export_GRN_from_ResultFlow_Component.py
conda deactivate
```

---

## Usage with Custom Data

### Step 1: Prepare Gene Metadata
Create a gene metadata file describing all genes to be modeled in the GRN.

Required format:
```text
TF_NAME    SGD_ID    GENE_LENGTH
```

Example:
```text
ABF1    YKL112W    2196
AFT1    YGL071W    2073
MET4    YNL103W    2019
```

---

### Step 2: Prepare Normalized RNA-seq Input
RNA-seq data are required.

We recommend TMM normalization, and an example workflow is provided in:
`salmon_tximport_tmm_normalize.R`

Other normalization methods are also supported, provided the final input format remains compatible.

Expression matrix (`-r`) requirements:
- first column: `samples`
- remaining columns: SGD systematic gene IDs
- rows: individual RNA-seq samples
- values: normalized expression values

Example:
```text
samples    YKL112W    YGL071W    YPL202C
WT_batch_1_rep_1    120.5    45.2    18.9
WT_batch_2_rep_1    118.3    44.7    19.1
AZF1_batch_1_rep_1  12.4     39.8    5.2
```

Requirements:
- the `samples` header is required
- sample names must use underscore-delimited formatting
- the text before the first underscore is treated as the genotype label

Recommended naming convention:
```text
GENOTYPE_batch_batch#_rep_replicate#
```

---

### Step 3: Run Input Preprocessing
Use:

```bash
python GRN_input_acquisition.py

Required:
-g: gene metadata file
-r: normalized RNA-seq expression matrix

Optional:
-p: p-value cutoff for discrete state identification (default: 0.01)
-t: TF–DNA prior (json format, see ./data/Rossi_Ruihao_TF_DNA_union_motif_based.json)
-c: protein–protein colocalization prior (json format, see ./data/PPI_network_Cutoff_0_STRING_overlapping_motif_sites_0_2025_union.json)
-a: genome annotation (./data/Sc_genome_annotations.txt)
-b: YEP replicate mapping (./data/YEP_best_rep.txt)
```

This generates SETIA-compatible preprocessing outputs in `./data`.

Optional prior inputs include:
- TF–DNA prior
- protein colocalization prior
- genome annotation
- YEP replicate mapping

If omitted, the corresponding preprocessing steps are automatically skipped.

### Step 3: Run Input Preprocessing
Use `GRN_input_acquisition_v2.py` to convert the gene metadata file, normalized RNA-seq matrix, and optional prior information into SETIA-compatible preprocessing files. The script writes its outputs primarily to `./data`, and also generates summary files in `./result`. The main generated files are:  

**Core outputs used by SETIA**
- `./data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt`  
  Discretized steady-state expression matrix used as the transcriptional-state input for GRN inference.
- `./data/GRN_ssTFs_Salmon_SteadyStates_2025_std.txt`  
  Per-gene standard deviations associated with the discretized steady states.
- `./data/GRN_ssTFs_Sc_gene_length.txt`  
  Gene lengths for the modeled genes.
- `./data/GRN_ssTFs_row_names.txt`  
  Row labels corresponding to genotype / condition names.
- `./data/GRN_ssTFs_column_names.txt`  
  Column labels corresponding to the modeled genes.

**Optional prior-derived outputs**
- `./data/GRN_ssTFs_Sc_TF_DNA.txt`  
  TF–DNA prior adjacency matrix generated when `-t` is provided.
- `./data/GRN_ssTFs_Sc_LG.txt`  
  Local-group / complex-structure file derived from promoter binding and colocalization information when the corresponding optional prior files are provided.
- `./data/GRN_ssTFs_Sc_promoter_strength.txt`  
  Promoter strength matrix generated when the required promoter-annotation resources are available.
- `./data/GRN_ssTFs_Sc_TF_DNA_TPM_union.txt`  
  Union of TF–DNA prior and TPM-derived relationships, generated when `-t` is provided.
- `./data/GRN_ssTFs_Sc_initial_condition.txt`  
  Initial GRN edge-state string derived from TF–DNA and expression similarity, generated when `-t` is provided.

**Additional preprocessing summaries**
- `./data/GRN_ssTFs_Salmon_SteadyStates_2025.txt`  
  Average / grouped steady-state expression values before discretization.
- `./result/Steady_state_count.txt`  
  Summary of the number of inferred discrete states per gene.
- `./result/GMM_figures/AIC/`  
  Diagnostic plots for discrete-state identification.

If optional inputs such as `-t`, `-c`, `-a`, or `-b` are not provided, the corresponding prior-dependent preprocessing steps are skipped automatically.

---

### Step 4: Run SETIA Inference
Run the main evolutionary GRN inference:

```bash
python EGRN_Multi_Genalg_Combinatorial_2025.py
```

using the processed inputs from Step 3.

---

## File Descriptions

### `salmon_tximport_tmm_normalize.R`
Performs TMM normalization on Salmon-aligned RNA-seq samples.

### `EGRN_Multi_Genalg_Combinatorial_2025.py`
Main SETIA executable implementing evolutionary optimization of ODE-based GRNs.

### `GRN_input_acquisition_v2.py`
Preprocesses RNA-seq and optional prior data into SETIA-compatible inputs.

### `GRN_Dynamic_Simulator_Combinatorial_remove_dispensible_edges_2025.py`
Removes dispensable GRN edges that do not affect stable-state recovery.

### `Export_final_GRN_TF_DNA.py`
Exports the final GRN and compares it against the TF–DNA prior.

### `EGRN_Multi_non_TF_2025.py`
Extends the core ssTF GRN to downstream non-TF targets.

### `Export_GRN_from_ResultFlow_Component.py`
Recovers and exports GRNs from intermediate saved optimization states.
