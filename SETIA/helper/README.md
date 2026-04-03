## Table of Contents

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
