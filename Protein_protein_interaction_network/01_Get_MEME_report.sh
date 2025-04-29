#!/bin/bash
# ---------------------- Gather MEME Motif reports from all factors -------------------
# ---------------------- rl884@cornell.edu, 25/04/28 ----------------------------------

# Make sure the destination directory exists
mkdir -p ./Motifs_all

# Loop over each folder that matches *_YEP under ./
for dir in ./*_YEP/; do
    if [ -d "$dir" ]; then
        # Find the MEME_Motifs.txt file inside
        file=$(find "$dir" -maxdepth 1 -name "*_MEME_Motifs.txt")

        if [ -f "$file" ]; then
            echo "Copying $file to ./Motifs_all/"
            cp "$file" ./Motifs_all/
        else
            echo "Warning: No *_MEME_Motifs.txt found in $dir"
        fi
    fi
done

echo "All motif files copied!"
