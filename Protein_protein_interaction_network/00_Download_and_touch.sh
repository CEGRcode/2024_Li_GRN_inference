#!/bin/bash
# ---------------------- Download the YEP data from archive and update timestamp ------------------------
# ---------------------- rl884@cornell.edu, 24/04/16 ----------------------------------------------------

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

# Downlad protein-protein interaction data in yeast from STRING
wget https://stringdb-downloads.org/download/protein.links.v12.0/4932.protein.links.v12.0.txt.gz \
     https://stringdb-downloads.org/download/protein.info.v12.0/4932.protein.info.v12.0.txt.gz

gunzip 4932.protein.links.v12.0.txt.gz
gunzip 4932.protein.info.v12.0.txt.gz
