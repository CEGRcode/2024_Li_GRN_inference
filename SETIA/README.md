# Gene regulatory network inference using SETIA

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA

## Overview
SETIA integrates multiple layers of molecular evidence:

- RNA-seq (transcriptional states across genotypes)
- ChIP-exo (TF–DNA binding)
- PRO-seq / CAGE-seq (promoter activity)
- Protein–protein colocalization (PPC)

These data are used to construct and constrain GRNs whose dynamics are modeled using ordinary differential equations (ODEs).

---

## Key Features

- Infers dynamical GRNs that reproduce cell-type-specific stable transcriptional states  
- Integrates multi-omics data for mechanistic constraints  
- Supports combinatorial regulation (AND / OR logic)  
- Implements parallelized optimization for scalability  
- Provides an interactive web platform for simulation and visualization  

---

## Repository Structure

```text
.
├── SETIA/                                #Core GRN inference and simulation framework
├── TF_DNA_Binding_Network/               #TF–DNA regulatory networks
├── Protein_protein_interaction_network/  #PPC-based regulation
├── GRN_simulator_website/                #Web interface
└── README.md
```

---

## Pipeline

1. `0_Download_and_touch.sh`  
   → downloads and initializes required data  

2. `1_Preprocessing...`  
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
conda create -n EvoAlg -c conda-forge -y python=3.10 numpy pandas scipy matplotlib scikit-learn pip
conda activate EvoAlg
pip install flask tqdm seaborn joblib
```


## Minimal Example

To test the pipeline without full datasets:

```bash
cd SETIA
conda activate EvoAlg
# -r: input transcriptional profiles
# -n: time scale in solveIVP
# -i: iteration times for model training
# -p: perturbation setting to the initial states
# -t: input promoter strength settings
# -l: input gene lengths
# -o: prefix for output files
# -e: random seed
# -f: allow edges not supported by ChIP
# -k: whether to pickle GRN instances
python EGRN_Multi_Genalg_Combinatorial_2025.py \
    -r data/example_RNA_profiles.txt \
    -n 3000 \
    -i 800 \
    -p 0 \
    -t data/example_promoter_strength.txt \
    -l data/example_gene_length.txt \
    -o "Minimal_example_" \
    -e 42 \
    -f 0 \
    -k 0
```

## Outputs
The inferred gene regulatory network (GRN) is written to the `./result` directory. The final GRN structure is saved as `GRN_filtered_Sc_GRN_final_2.json`, and the corresponding model parameters are provided in `GRN_filtered_full_Sc_GRN_final.txt`. Intermediate and final GRN instances are serialized and stored as pickle files, while all searched or simulated GRN configurations during inference are cached in the `sys_cache_*` files.

## Table of Contents

### GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py
Searches for stable states of a given gene regulatory network (GRN) using a local exploration strategy. The simulator initializes the system from a specified set of initial transcriptional states and evolves the GRN dynamics forward in time to determine whether each trajectory converges to a stable or unstable final state. This mode is intended to assess attractor stability and basin structure near biologically relevant initial conditions, such as the transcriptional profiles of stable cells.

### GRN_Dynamic_Simulator_Combinatorial_Global_multistate_2025.py
Searches for stable states of a given GRN using a global exploration strategy. The simulator samples a large, evenly distributed set of initial states across the full expression space, where each dimension corresponds to a gene’s expression level, and evolves the GRN dynamics from each starting point. This mode is designed to systematically map the global attractor landscape of the GRN and identify all accessible stable states, and therefore will be computationally expensive.

### salmon_tximport_tmm_normalize
Performs the TMM normalization on the RNA-seq samples aligned by Salmon. It produces the normalized read counts in GRN_Sc_TMM_normalized_CPM.txt for downstream analysis.

### GRN_Expanded_Combinatorial_2025.py
Contains classes and functions for representing GRN instances, including their kinetic parameters and associated attributes.

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

### EGRN_Multi_non_TF_2025.py
Extends the core ssTF GRN to downstream non-TF target genes by fitting their stable expression states using ssTF-to-target edges supported by TF–DNA binding evidence.

