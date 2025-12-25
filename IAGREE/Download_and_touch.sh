#!/bin/bash
# ---------------------- Download the YEP data from archive and update timestamp ------------------------
# ---------------------- rl884@cornell.edu, 24/04/16 ----------------------------------------------------

cd ./data

# Downlad YEP data
wget -r -np -nH --cut-dirs=5 --timestamping --reject "index.html*" "https://www.datacommons.psu.edu/download/eberly/pughlab/yeast-epigenome-project/"

mkdir ./masterNoTag_20180928
for file in ./*_YEP.zip ./masterNoTag_20180928.zip; do
    # Unzip the file
    unzip "$file" -d ./
    
    # Find all files in the unzipped directory and touch them to update timestamp
    find "${file%.zip}" -type f -exec touch {} \;
done

# Put masterNoTag in folder.
mv masterNoTag_20180928.bam ./masterNoTag_20180928
find ./masterNoTag_20180928/ -exec touch {} +
rm -- ./*.zip

# Copy all ChExMix filtered peaks to a new folder.
mkdir -p ./data/YEP_bed

for d in ./data/*_YEP/; do
  for bed in "$d"/*_chexmix_filtered_peaks.bed; do
    [ -e "$bed" ] || continue
    cp "$bed" ./data/YEP_bed/
  done
done

# Remove all YEP data.
find . -mindepth 1 -maxdepth 1 \
  ! -name "YEP_bed" \
  ! -name "All_gene_ID.txt" \
  ! -name "Sc_genome_annotations.txt" \
  ! -name "YEP_best_rep.txt" \
  ! -name "ssTFs_MATa_Spots_76.txt" \
  ! -name "Rossi_Ruihao_TF_DNA_union_motif_based.json" \
  ! -name "Rossi_Ruihao_TF_DNA_union_all.json" \
  ! -name "PPI_network_Cutoff_0_STRING_overlapping_motif_sites_0_2025_union.json" \
  -exec rm -rf {} +

# Download Supplementary data from Rossi 2021 and Lu 2020.
wget -O Supplementary_Data_2.xlsx \
"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03314-8/MediaObjects/41586_2021_3314_MOESM4_ESM.xlsx"

wget -O 41586_2021_3314_MOESM3_ESM.xlsx \
"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03314-8/MediaObjects/41586_2021_3314_MOESM3_ESM.xlsx"

curl -L \
  -A "Mozilla/5.0" \
  -o Supplemental_Data_S5_S8.xlsx \
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6633255/bin/supp_gr.245456.118_Supplemental_Data_S5_S8.xlsx

mv ./Supplementary_Data_2.xlsx ./41586_2021_3314_MOESM3_ESM.xlsx ./Supplemental_Data_S5_S8.xlsx ./data/

cd ..