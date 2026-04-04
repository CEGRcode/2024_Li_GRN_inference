#!/bin/bash
# ---------------------- Download the YEP data from archive and update timestamp ------------------------
# ---------------------- rl884@cornell.edu, 24/04/16 ----------------------------------------------------

cd ./data

# Download RNAseq data (if not public yet, use the reviewer token to download the GSE317148_GRN_Sc_TMM_normalized_CPM.txt.gz)
#wget -O GSE317148_GRN_Sc_TMM_normalized_CPM.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317148/suppl/GSE317148_GRN_Sc_TMM_normalized_CPM.txt.gz
gunzip -k GSE317148_GRN_Sc_TMM_normalized_CPM.txt.gz
mv GSE317148_GRN_Sc_TMM_normalized_CPM.txt GRN_Sc_TMM_normalized_CPM.txt

# Download YEP data
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
mkdir -p ./YEP_bed

for d in ./*_YEP/; do
  for bed in "$d"/*_chexmix_filtered_peaks.bed; do
    [ -e "$bed" ] || continue
    cp "$bed" ./YEP_bed/
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

# Download Supplementary data from Rossi 2021.
wget -O Supplementary_Data_2.xlsx \
"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03314-8/MediaObjects/41586_2021_3314_MOESM4_ESM.xlsx"

wget -O 41586_2021_3314_MOESM3_ESM.xlsx \
"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03314-8/MediaObjects/41586_2021_3314_MOESM3_ESM.xlsx"

# Download the CAGE-seq data from Lu 2020 for the promoter strength estimation.
# 1) Download the PMC OA package tarball
curl -L -o PMC6633255.tar.gz \
  "https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/94/02/PMC6633255.tar.gz"

# 2) Extract
mkdir -p PMC6633255_pkg
tar -xzf PMC6633255.tar.gz -C PMC6633255_pkg

# 3) Find the supplemental S5–S8 file inside
find PMC6633255_pkg -type f -iname "*Supplemental*Data*S5*S8*.xlsx" -o -iname "*S5*S8*.xlsx"

# 4) Copy it to your working dir (adjust the pattern if needed)
cp "$(find PMC6633255_pkg -type f -iname "*Supplemental*Data*S5*S8*.xlsx" | head -n 1)" \
  supp_gr.245456.118_Supplemental_Data_S5_S8.xlsx

rm PMC6633255.tar.gz
rm -rf ./PMC6633255_pkg

# 5) Verify it's a real OOXML .xlsx (ZIP container)
mv supp_gr.245456.118_Supplemental_Data_S5_S8.xlsx Supplemental_Data_S5_S8.xlsx


cd ..
