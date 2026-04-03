## Table of Contents

### GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py
Searches for stable states of a given gene regulatory network (GRN) using a local exploration strategy. The simulator initializes the system from a specified set of initial transcriptional states and evolves the GRN dynamics forward in time to determine whether each trajectory converges to a stable or unstable final state. This mode is intended to assess attractor stability and basin structure near biologically relevant initial conditions, such as the transcriptional profiles of stable cells.

### GRN_Dynamic_Simulator_Combinatorial_Global_multistate_2025.py
Searches for stable states of a given GRN using a global exploration strategy. The simulator samples a large, evenly distributed set of initial states across the full expression space, where each dimension corresponds to a gene’s expression level, and evolves the GRN dynamics from each starting point. This mode is designed to systematically map the global attractor landscape of the GRN and identify all accessible stable states, and therefore will be computationally expensive.

### GRN_Expanded_Combinatorial_2025.py
Contains classes and functions for representing GRN instances, including their kinetic parameters and associated attributes.

### dynamics.py
Contains helper functions for the dynamics module.
