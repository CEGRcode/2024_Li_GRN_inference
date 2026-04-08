#!/bin/bash
# ---------------------- MAKE BEDGRAPH FOR EACH BAM ----------------------
# ---------------------- rl884@cornell.edu, 24/08/89 ----------------------------------------------------

set -e

module load PPI

OUTDIR="../Protein_protein_interaction_network/static/data"

for folder in ./*_YEP/; do
    [ -d "$folder" ] || continue

    dirname=$(basename "$folder")
    ID="${dirname%%_*}"

    bam="${folder}${ID}_filtered_sorted.bam"
    [ -e "$bam" ] || continue

    java -jar ScriptManager-v0.14.jar \
        bam-format-converter bam-to-bedgraph \
        -o "${folder}${ID}" \
        "$bam"

    mv "${folder}"*.bedGraph "$OUTDIR"/
done
