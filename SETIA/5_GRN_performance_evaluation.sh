#!/bin/bash
# ------------------------------------------------------------------
# Evaluate GRN performance in reproducing transcriptional profiles
# as stable states
#
# Ruihao Li (rl884@cornell.edu)
# 2024-04-16
# ------------------------------------------------------------------

python GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py \
  -j ./result/GRN_filtered_full_Sc_GRN_final.txt \
  -r ./data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt \
  -l ./data/GRN_ssTFs_Sc_gene_length.txt \
  -t ./data/GRN_ssTFs_Sc_promoter_strength.txt \
  -o GRN_prediction \
  -m 0

# Arguments:
# -j : inferred GRN file
# -r : transcriptional profiles to test
# -l : gene length file
# -t : promoter strength file
# -o : output prefix
# -m : if not 0, use a random GRN
