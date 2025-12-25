#!/bin/bash
# ---------------------- Preprocessing data and discrete gene expression state identification -----------
# ---------------------- rl884@cornell.edu, 25/12/16 ----------------------------------------------------


python GRN_input_acquisition.py \
  -p 0.01 \
  -g ./data/All_gene_ID.txt \
  -r ./data/GRN_Sc_TMM_normalized_CPM.txt \
  -t ./data/Rossi_Ruihao_TF_DNA_union_motif_based.json \
  -c ./data/PPI_network_Cutoff_0_STRING_overlapping_motif_sites_0_2025_union.json \
  -a ./data/Sc_genome_annotations.txt \
  -b ./data/YEP_best_rep.txt