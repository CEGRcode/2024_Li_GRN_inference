#!/bin/bash

# ---------------------- Building TF-DNA binding network from YEP data-----------------------------------
# ---------------------- rl884@cornell.edu, 25/04/16 ----------------------------------------------------

# Exit script if any command fails
set -e

echo "Step 0: Downloading processed data from Rossi Nature 2019..."
bash Download_and_touch.sh

echo "Step 1: Processing the supplementary table 3 from Rossi Nature 2019..."
python ./Dissect_full_spreadsheet.py

echo "Step 2: Identifying bidirectional transcription regions and TFs in between..."
python ./Identify_bidirectional_transcription_regions_and_TF_inbetween.py

echo "Step 3: Building the TF-DNA binding network..."
python ./TF_DNA_binding_network.py

echo "Step 4: Cleaning up intermediate results..."
rm ./result/*.txt

echo "TF-DNA binding network built successfully."
