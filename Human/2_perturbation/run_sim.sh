#!/bin/bash
NAME=$1
KO=$2
OE=$3

# sim.sub joins multi-gene indices with '-' (e.g. 5-48) so the comma inside a
# value does not collide with HTCondor's comma field delimiter. Translate back
# to the comma form the simulator's -K/-O expect (e.g. 5,48).
[ "$KO" != "NONE" ] && KO=$(echo "$KO" | tr '-' ',')
[ "$OE" != "NONE" ] && OE=$(echo "$OE" | tr '-' ',')

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

tar -xzf SETIA_sim_code.tar.gz
mkdir -p SETIA/result
cp grn_output/GRN_filtered_full_*.txt SETIA/result/

cd SETIA

PERT=""
[ "$KO" != "NONE" ] && PERT="$PERT -K $KO"
[ "$OE" != "NONE" ] && PERT="$PERT -O $OE"

python3 -u GRN_Dynamic_Simulator_parallel.py \
    -r data/setia_input_linear_CPM_pseudocount.tsv \
    -t data/BMMC_promoter_strength.txt \
    -l data/BMMC_gene_length.txt \
    -j result/GRN_filtered_full_BMMC_f1_GRN_final_raw.txt \
    -o "${NAME}" -n 3000 $PERT

cd ..
mkdir -p sim_result
cp SETIA/result/GRN_dynamic_local_search_result_${NAME}.txt sim_result/
