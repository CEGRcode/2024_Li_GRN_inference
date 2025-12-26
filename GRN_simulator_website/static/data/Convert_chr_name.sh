#!/bin/bash

# Define the dictionary
declare -A change_dic=(
    ["chr1"]="chrI" ["chr2"]="chrII" ["chr3"]="chrIII" ["chr4"]="chrIV"
    ["chr5"]="chrV" ["chr6"]="chrVI" ["chr7"]="chrVII" ["chr8"]="chrVIII"
    ["chr9"]="chrIX" ["chr10"]="chrX" ["chr11"]="chrXI" ["chr12"]="chrXII"
    ["chr13"]="chrXIII" ["chr14"]="chrXIV" ["chr15"]="chrXV" ["chr16"]="chrXVI"
    ["2-micron"]="chrM"
)

# Define the folder path
folder_path="/home/ubuntu/GRN_UI_app_test/static/data/"

# Process each .bedGraph file in the directory
for file_path in "$folder_path"*.bedGraph; do
    # Create an output file with the new name
    output_file="${file_path%.bedGraph}_processed.bedGraph"
    
    # Use awk to process the file
    awk -v OFS='\t' -v output_file="$output_file" '
    BEGIN {
        # Populate the change_dic in awk
        change_dic["chr1"]="chrI"; change_dic["chr2"]="chrII"; change_dic["chr3"]="chrIII"; 
        change_dic["chr4"]="chrIV"; change_dic["chr5"]="chrV"; change_dic["chr6"]="chrVI";
        change_dic["chr7"]="chrVII"; change_dic["chr8"]="chrVIII"; change_dic["chr9"]="chrIX"; 
        change_dic["chr10"]="chrX"; change_dic["chr11"]="chrXI"; change_dic["chr12"]="chrXII"; 
        change_dic["chr13"]="chrXIII"; change_dic["chr14"]="chrXIV"; change_dic["chr15"]="chrXV";
        change_dic["chr16"]="chrXVI"; change_dic["2-micron"]="chrM";
    }
    {
        # Check if the first column matches a key in change_dic
        if ($1 in change_dic) {
            $1 = change_dic[$1]
        }
        print $0 > output_file
    }' "$file_path"
done
