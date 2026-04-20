#!/bin/bash
# ---------------------- Connect the ssTF core GRN to non-TF targets by TF-DNA binding evidence and optimal fit to RNA-seq stable states ------------------------
# -------------------------------------------------- rl884@cornell.edu, 24/04/16 --------------------------------------------------------------------------------
# Important note: the -r -c -y -g files here should contain data for ALL genes, not just the ssTFs. Rerun 1_Preprocessing_dicrete_gene_expression_state_identification.sh for ALL genes (-g ./data/All_gene_ID.txt).
# Otherwise an empty network will be generated.
for j in $(seq 0 49); do
    python EGRN_Multi_non_TF_2025.py \
        -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt \
        -n 3000 \
        -c data/GRN_ssTFs_Sc_TF_DNA.txt \
        -x ./data/GRN_ssTFs_TF_set.txt \
        -y ./data/GRN_ssTFs_column_names.txt \
        -g ./data/GRN_ssTFs_Sc_LG.txt \
        -j $j &
done
wait
python ./helper/merge_GRN.py
rm ./result/GRN_nonTF_{0..49}.json
