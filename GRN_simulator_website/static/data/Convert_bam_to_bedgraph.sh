# ---------------------- MAKE BEDGRAPH FOR EACH BAM ----------------------
# ---------------------- rl884@cornell.edu, 24/08/89 ----------------------------------------------------

for file in *.bam; do
	# get sample ID.
	ID="${file%%_*}"
	
	# samtools sort BAM files.
	samtools sort -o "${ID}_filtered_sorted.bam" "${ID}_filtered.bam"

	# samtools index BAM files.
	samtools index "${ID}_filtered_sorted.bam"

	# expand the bed file.
	java -jar ScriptManager-v0.14.jar bam-format-converter bam-to-bedgraph -o "${ID}" "${ID}_filtered_sorted.bam"

	# clean-up.
	rm "${ID}_filtered.bam"
done
