#!/bin/bash
# ---------------------- Build protein-protein interaction network --------------------------------------
# ---------------------- rl884@cornell.edu, 25/04/29 ----------------------------------------------------

# To build motif-based PPI
#python PPI_network.py "True"

# To build PPI for all binding sites
python PPI_network.py "False"