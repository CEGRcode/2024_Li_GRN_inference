#!/bin/bash
set -e

module load PPI

INPUT_DIR="../Protein_protein_interaction_network"

bash ./static/data/Convert_bam_to_bedgraph.sh "$INPUT_DIR"
python ./static/data/Combine_cdt.py "$INPUT_DIR/SCALE_cdt" ./static/data
python ./static/data/Combine_h5.py ./static/data

# clean intermediate files
rm -f ./static/data/*_BED.h5
