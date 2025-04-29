#!/bin/bash
# ---------------------- Get the motifs and binding sites for all factors -------------
# ---------------------- rl884@cornell.edu, 25/04/28 ----------------------------------

############################################# Get MEME motifs covering most peaks #########################################
# Input file
file="QC_count.txt"

# Create output directory if it doesn't exist
mkdir -p ./YEP_Motif_BED

# Process each line in the file
while IFS= read -r line; do
    # Extract the filename (last field in the line)
    filename=$(echo "$line" | awk '{print $NF}')

    # Extract the sample ID (everything before "_Motif")
    sample_id=$(echo "$filename" | sed -n 's/_Motif.*//p')

    # Extract the motif number (number between "Motif_" and "_peaks")
    motif_number=$(echo "$filename" | grep -oP 'Motif_\K[0-9]+' || echo "None")

    # Check if motif number exists
    if [[ $motif_number == "None" ]]; then
        echo "No motif number found for sample ID $sample_id. Skipping..."
        continue
    fi

    # Construct the folder and file paths
    folder_path="./${sample_id}_YEP_Done"
    file_to_copy="${sample_id}_Motif_${motif_number}_bound.bed"

    # Check if the folder exists
    if [[ -d $folder_path ]]; then
        # Copy the file
        src_file="$folder_path/$file_to_copy"
        if [[ -f $src_file ]]; then
            cp "$src_file" ./YEP_Motif_BED
            echo "Copied $src_file to current directory."
        else
            echo "File $file_to_copy not found in $folder_path. Skipping..."
        fi
    else
        echo "Folder $folder_path not found. Skipping..."
    fi
done < "$file"
###########################################################################################################################


########################################### Get all binding sites called by CHExMix #######################################
# Create output directory if it doesn't exist
mkdir -p ./YEP_ALL_BED

# Loop through each folder under ./ that matches *YEP
for dir in ./*YEP/; do
    if [ -d "$dir" ]; then
        echo "Processing folder: $dir"

        # Find the file matching *chexmix_peaks.bed inside the folder
        bed_file=$(find "$dir" -maxdepth 1 -name "*chexmix_peaks.bed")

        if [ -f "$bed_file" ]; then
            echo "  Found BED file: $bed_file"

            # Copy it into ./YEP_ALL_BED
            cp "$bed_file" ./YEP_ALL_BED/
        else
            echo "  No BED file found in $dir"
        fi
    fi
done

echo "All done. BED files collected in ./YEP_ALL_BED/"
###########################################################################################################################