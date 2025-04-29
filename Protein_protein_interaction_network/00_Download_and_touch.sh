#!/bin/bash
# ---------------------- Download the YEP data from archive and update timestamp ------------------------
# ---------------------- rl884@cornell.edu, 24/04/16 ----------------------------------------------------

# Downlad YEP data
wget -r -np -nH --cut-dirs=5 --timestamping --reject "index.html*" "https://www.datacommons.psu.edu/download/eberly/pughlab/yeast-epigenome-project/"

for file in ./*.zip; do
    # Unzip the file
    unzip "$file" -d ./
    
    # Find all files in the unzipped directory and touch them to update timestamp
    find "${file%.zip}" -type f -exec touch {} \;
done

# Put masterNoTag in folder.
mkdir ./masterNoTag_20180928
mv masterNoTag_20180928.bam ./masterNoTag_20180928