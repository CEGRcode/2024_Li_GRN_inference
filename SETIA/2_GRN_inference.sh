#!/bin/bash
# ---------------------- GRN inference was performed on HPC environment ------------------------

#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --mem=32GB
#SBATCH --time=48:00:00
#SBATCH --partition=open
#SBATCH --array=0-9
#SBATCH --exclude=p-sc-2521

seeds=(42 123 456 789 101112 131415 161718 192021 222324 252627)
outputs=("Sc_0" "Sc_1" "Sc_2" "Sc_3" "Sc_4" "Sc_5" "Sc_6" "Sc_7" "Sc_8" "Sc_9")

# Get the corresponding values for this task
seed=${seeds[$SLURM_ARRAY_TASK_ID]}
output=${outputs[$SLURM_ARRAY_TASK_ID]}

# Activate the conda environment
conda activate EvoAlg

# Run the Python script with the specified parameters
python EGRN_Multi_Genalg_Combinatorial_2025.py \
    -c data/GRN_ssTFs_Sc_TF_DNA.txt \
    -g data/GRN_ssTFs_Sc_LG.txt \
    -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt \
    -b data/GRN_ssTFs_Salmon_SteadyStates_2025_std.txt \
    -n 3000 \
    -i 200 \
    -p 0 \
    -t data/GRN_ssTFs_Sc_promoter_strength.txt \
    -l data/GRN_ssTFs_Sc_gene_length.txt \
    -o "${output}" \
    -e "${seed}" \
    -d data/GRN_ssTFs_Sc_initial_condition.txt \
    -f 0 \
    -k 0

# Deactivate the conda environment
conda deactivate
