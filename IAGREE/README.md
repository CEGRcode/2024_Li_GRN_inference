# Gene regulatory network inference using IAGREE

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

## Directions
To recreate the gene regulatory network for this manuscript, please execute the shell scripts in order: `0_Download_and_touch.sh`, `1_Preprocessing_dicrete_gene_expression_state_identification.sh`, `2_GRN_inference.sh`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n EvoAlg -c conda-forge -y python=3.10 numpy pandas scipy matplotlib scikit-learn pip
```

## Table of Contents

### salmon_tximport_tmm_normalize
Performs the TMM normalization on the RNA-seq samples aligned by Salmon. It produces the normalized read counts in GRN_Sc_TMM_normalized_CPM.txt for downstream analysis.

### GRN_Expanded_Combinatorial_2025.py
Contains classes and functions for representing gene regulatory network (GRN) instances, including their kinetic parameters and associated attributes.

### dynamics.py
Contains helper functions for the dynamics module.

### distance_functions.py
Contains helper functions for computing Hamming distance and mean normalized L1 distance at the per-gene and per-profile levels.

### mutation_functions.py
Contains functions for mutating gene regulatory network configurations under different strategies, including TF–DNA–constrained mutations and unconstrained (free-search) mutations.

### population_functions.py
Contains functions for population-level evolutionary dynamics of gene regulatory networks, including fitness-based proportion updates, natural selection, long-jump mutation selection, and initialization of global mutation direction parameters.

### utility_functions.py
Miscellaneous helper functions.

### Export_GRN_from_ResultFlow_Component.py
Exports gene regulatory networks from intermediate saved states when inference stops before the pre-set number of iterations. Transient GRN states are saved at each iteration and can be recovered using this script.

