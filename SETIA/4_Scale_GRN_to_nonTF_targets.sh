#!/bin/bash
# ---------------------- Connect the ssTF core GRN to non-TF targets by TF-DNA binding evidence and optimal fit to RNA-seq stable states ------------------------
# -------------------------------------------------- rl884@cornell.edu, 24/04/16 --------------------------------------------------------------------------------

python EGRN_Multi_non_TF_2025.py -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt -n 3000 -c data/GRN_ssTFs_Sc_TF_DNA.txt -x ./data/GRN_ssTFs_TF_set.txt -y ./data/GRN_ssTFs_column_names.txt -g ./data/GRN_ssTFs_Sc_LG.txt