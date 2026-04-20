conda activate EvoAlg

############################## Optional ##############################
# Remove dispensible edges from the GRN.
#python GRN_Dynamic_Simulator_Combinatorial_remove_dispensible_edges_2025.py \
#    -c data/GRN_ssTFs_Sc_TF_DNA.txt \
#    -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt \
#    -g data/GRN_ssTFs_column_names.txt \
#    -n 3000 \
#    -i result/ \
#    -p 0 \
#    -o Sc_beta_complete
############################## Optional ##############################

# Finally export the GRN in comparison against TF-DNA binding prior.
python Export_final_GRN_TF_DNA.py \
    -c data/GRN_ssTFs_Sc_TF_DNA.txt \
    -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt \
    -g data/GRN_ssTFs_column_names.txt \
    -n 3000 \
    -i result/ \
    -p 0 \
    -o Sc_GRN_final

conda deactivate
