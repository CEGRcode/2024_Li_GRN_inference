# Gene regulatory network inference using IAGREE

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

## Directions
To recreate the gene regulatory network for this manuscript, please execute the shell scripts in order: `0_Download_and_touch.sh`, `1_Preprocessing_dicrete_gene_expression_state_identification.sh`, `2_GRN_inference.sh`, and `3_Final_GRN_construction.sh`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n EvoAlg -c conda-forge -y python=3.10 numpy pandas scipy matplotlib scikit-learn pip
```

## Outputs
The inferred gene regulatory network (GRN) is written to the `./result` directory. The final GRN structure is saved as `GRN_filtered_Sc_GRN_final_2.json`, and the corresponding model parameters are provided in `GRN_filtered_full_Sc_GRN_final.txt`. Intermediate and final GRN instances are serialized and stored as pickle files, while all searched or simulated GRN configurations during inference are cached in the `sys_cache_*` files.

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

### merge_GRN.py
Merged two GRNs.

### GRN_Dynamic_Simulator_Combinatorial_remove_dispensible_edges_2025.py
Removes dispensable edges from a gene regulatory network that do not affect the reproduction of transcriptional profiles as stable states.

### Export_final_GRN_TF_DNA.py
Exports the final gene regulatory network and compares it against the TF–DNA binding prior, highlighting shared edges in green.
