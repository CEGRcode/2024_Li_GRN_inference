#!/bin/bash
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

    java -jar "$INPUT_ROOT/ScriptManager-v0.14.jar" \
        bam-format-converter bam-to-bedgraph \
        -o "${folder}${ID}" \
        "$bam"

    mv "${folder}${ID}.bedGraph" "$OUTDIR"/
done
