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

# Download Supplementary data from Rossi Nature 2021.
wget -O Supplementary_Data_2.xlsx \
"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03314-8/MediaObjects/41586_2021_3314_MOESM4_ESM.xlsx"

wget -O 41586_2021_3314_MOESM3_ESM.xlsx \
"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03314-8/MediaObjects/41586_2021_3314_MOESM3_ESM.xlsx"

mv ./Supplementary_Data_2.xlsx ./41586_2021_3314_MOESM3_ESM.xlsx ./data/



cd ..
