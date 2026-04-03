## Table of Contents

### GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py
Searches for stable states of a given gene regulatory network (GRN) using a local exploration strategy. The simulator initializes the system from a specified set of initial transcriptional states and evolves the GRN dynamics forward in time to determine whether each trajectory converges to a stable or unstable final state. This mode is intended to assess attractor stability and basin structure near biologically relevant initial conditions, such as the transcriptional profiles of stable cells.

### GRN_Dynamic_Simulator_Combinatorial_Global_multistate_2025.py
Searches for stable states of a given GRN using a global exploration strategy. The simulator samples a large, evenly distributed set of initial states across the full expression space, where each dimension corresponds to a gene’s expression level, and evolves the GRN dynamics from each starting point. This mode is designed to systematically map the global attractor landscape of the GRN and identify all accessible stable states, and therefore will be computationally expensive.

### salmon_tximport_tmm_normalize.R
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
Merges two GRNs.

### GRN_Dynamic_Simulator_Combinatorial_remove_dispensible_edges_2025.py
Removes dispensable edges from a gene regulatory network that do not affect the reproduction of transcriptional profiles as stable states.

### Export_final_GRN_TF_DNA.py
Exports the final gene regulatory network and compares it against the TF–DNA binding prior, highlighting shared edges in green.

### EGRN_Multi_non_TF_2025.py
Extends the core ssTF GRN to downstream non-TF target genes by fitting their stable expression states using ssTF-to-target edges supported by TF–DNA binding evidence.

### EGRN_Multi_Genalg_Combinatorial_2025.py
The main executable for SETIA
