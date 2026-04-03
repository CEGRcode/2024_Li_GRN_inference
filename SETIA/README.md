# Gene regulatory network inference using SETIA

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA

## Overview
SETIA integrates multiple layers of molecular evidence:

- RNA-seq (transcriptional states across genotypes)
- ChIP-exo (TF–DNA binding)
- PRO-seq / CAGE-seq (promoter activity)
- Protein–protein interaction (PPI) including protein-protein colocalization (PPC)

These data are used to construct and constrain GRNs whose dynamics are modeled using ordinary differential equations (ODEs).

---

## Key Features

- Infers dynamical GRNs that reproduce cell-type-specific stable transcriptional states  
- Integrates multi-omics data for mechanistic constraints  
- Supports combinatorial regulation (AND / OR logic)  
- Implements parallelized optimization for scalability  
- Provides an interactive web platform for simulation and visualization  

---

## Pipeline to reproduce GRNs from the paper

1. `0_Download_and_touch.sh`  
   → downloads and initializes required data  

2. `1_Preprocessing_dicrete_gene_expression_state_identification.sh`  
   → identifies discrete expression states  

3. `2_GRN_inference.sh`  
   → performs GRN inference  

4. `3_Final_GRN_construction.sh`  
   → filters and constructs final GRN  

5. `4_Scale_GRN_to_nonTF_targets.sh`  
   → extends network to non-TF genes

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

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

## Outputs
The inferred gene regulatory network (GRN) is written to the `./result` directory. The final GRN structure is saved as `GRN_filtered_Sc_GRN_final_2.json`, and the corresponding model parameters are provided in `GRN_filtered_full_Sc_GRN_final.txt`. Intermediate and final GRN instances are serialized and stored as pickle files, while all searched or simulated GRN configurations during inference are cached in the `sys_cache_*` files.

## Minimal Example

This example demonstrates how to run SETIA on a small test dataset. As the code is optimized for HPC platforms, running it on a local PC may take longer (approximately 10 minutes). The example inferred GRN will be written to `./result/GRN_filtered_Sc_GRN_final_raw.json`, and can be uploaded to our GRN simulator website (https://grn.cac.cornell.edu:5000/) for interactive exploration.

```bash
cd SETIA
conda activate EvoAlg
# -r: input transcriptional profiles
# -n: time span used in solve_ivp
# -i: number of training iterations
# -p: perturbation setting to the initial states
# -t: input promoter strength settings
# -l: input gene lengths
# -o: prefix for output files
# -e: random seed
# -f: allow edges not supported by ChIP
# -k: whether to pickle GRN instances
# The following command runs SETIA on a minimal example dataset and writes output files with the prefix Minimal_example_.
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



## Table of Contents

### salmon_tximport_tmm_normalize.R
Performs the TMM normalization on the RNA-seq samples aligned by Salmon. It produces the normalized read counts in GRN_Sc_TMM_normalized_CPM.txt for downstream analysis.

### EGRN_Multi_Genalg_Combinatorial_2025.py
The main executable for SETIA, implementing evolutionary optimization of ODE-based GRNs to recapitulate observed transcriptional stable states based on given GRN structural information.

### GRN_input_acquisition.py
Preprocess the RNA-seq, ChIP-exo, and other data for GRN inference in yeast.
