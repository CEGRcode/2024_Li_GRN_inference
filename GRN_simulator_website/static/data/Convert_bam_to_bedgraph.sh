#!/bin/bash
# ---------------------- MAKE BEDGRAPH FOR EACH BAM ----------------------
# ---------------------- rl884@cornell.edu, 24/08/29 ---------------------

set -e

module load PPI

INPUT_ROOT="../Protein_protein_interaction_network"
OUTDIR="./static/data"

mkdir -p "$OUTDIR"

for folder in "$INPUT_ROOT"/*_YEP/; do
    [ -d "$folder" ] || continue

    dirname=$(basename "$folder")
    ID="${dirname%%_*}"

    bam="${folder}${ID}_filtered_sorted.bam"
    [ -e "$bam" ] || continue

    java -jar ./static/data/ScriptManager-v0.14.jar \
        bam-format-converter bam-to-bedgraph \
        -o "${folder}${ID}" \
        "$bam"

    mv "${folder}"*.bedGraph "$OUTDIR"/
done
