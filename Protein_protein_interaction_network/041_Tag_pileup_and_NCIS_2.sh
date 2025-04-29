# ---------------------- MAKE COMPOSITE PROFILE FOR EACH PROTEIN TO ALL OTHER PROTEINS ALL BINDING SITES PPI-
# ---------------------- rl884@cornell.edu, 25/04/21 --------------------------------------------------------

# Loop through all directories (excluding symbolic links)
for folder in ./*_YEP/; do
    # Check if the item is a directory
    if [ -d "$folder" ]; then
        BED_ID="${folder%%_*}"
	BED_ID="${BED_ID##*/}"
	folder="./${BED_ID}_YEP/"

	for folder in ./*_YEP/ ./*_YEP_Done/; do
		# get sample ID.
		ID="${folder%%_*}"
		ID="${ID#./}"	
	
		# pile up tags and output in cdt format.
		java -jar ./ScriptManager-v0.14.jar read-analysis tag-pileup "${BED_ID}_YEP/${BED_ID}_chexmix_filtered_peaks_500bp_remove_convert.bed" "${folder}${ID}_filtered_sorted.bam" --cdt -M "${BED_ID}_BED_${ID}_BAM"

		# obtain the scaling factor.
		Factor=$(tail -n 1 "${folder}${ID}_filtered_sorted_ScalingFactors.out" | awk '{print $NF}')

		# apply the scaling factor to the cdt.
		java -jar ./ScriptManager-v0.14.jar read-analysis scale-matrix "${BED_ID}_BED_${ID}_BAM_sense.cdt" -s $Factor -r 2 -l 3 -o "./SCALE_cdt/${BED_ID}_BED_${ID}_BAM_sense_SCALE.cdt"
		java -jar ./ScriptManager-v0.14.jar read-analysis scale-matrix "${BED_ID}_BED_${ID}_BAM_anti.cdt" -s $Factor -r 2 -l 3 -o "./SCALE_cdt/${BED_ID}_BED_${ID}_BAM_anti_SCALE.cdt"
	
		# clean-up.
		rm "${BED_ID}_BED_${ID}_BAM_sense.cdt" "${BED_ID}_BED_${ID}_BAM_anti.cdt"
	done
	zip "./SCALE_cdt/${BED_ID}_zipped_BED_BAM_sense.zip" "./SCALE_cdt/${BED_ID}_BED_"*
	count=$(unzip -l "./SCALE_cdt/${BED_ID}_zipped_BED_BAM_sense.zip" | grep -c '^')
	echo "$BED_ID: $count" >> QC_count.txt
	rm "./SCALE_cdt/${BED_ID}_BED_"*
	mv "./${BED_ID}_YEP/" "./${BED_ID}_YEP_Done/"
    fi
done