#!/bin/bash
# ---------------------- MAKE COMPOSITE PROFILE FOR EACH PROTEIN ITSELF ALL BINDING SITES PPI -----------
# ---------------------- rl884@cornell.edu, 25/04/21 ----------------------------------------------------

module load PPI

samtools sort -o "./masterNoTag_20180928/masterNoTag_20180928_sorted.bam" "./masterNoTag_20180928/masterNoTag_20180928.bam"
samtools index "./masterNoTag_20180928/masterNoTag_20180928_sorted.bam"

for folder in ./*_YEP/; do
	# get sample ID.
	ID="${folder%%_*}"
	ID="${ID#./}"
	
	# samtools sort BAM files.
	samtools sort -o "${folder}${ID}_filtered_sorted.bam" "${folder}${ID}_filtered.bam"

	# samtools index BAM files.
	samtools index "${folder}${ID}_filtered_sorted.bam"

	# expand the bed file.
	java -jar ./ScriptManager-v0.14.jar coordinate-manipulation expand-bed -c 500 "${folder}${ID}_chexmix_filtered_peaks.bed" -o "${folder}${ID}_chexmix_filtered_peaks_500bp.bed"

	# remove the last two rows in the chexmix bed file.
	awk '{NF=NF-2}1' "${folder}${ID}_chexmix_filtered_peaks_500bp.bed" > "${folder}${ID}_chexmix_filtered_peaks_500bp_remove.bed"
	sed 's/$/\r/' "${folder}${ID}_chexmix_filtered_peaks_500bp_remove.bed" > "${folder}${ID}_chexmix_filtered_peaks_500bp_remove_temp.bed"
	sed 's/ /\t/g' "${folder}${ID}_chexmix_filtered_peaks_500bp_remove_temp.bed" > "${folder}${ID}_chexmix_filtered_peaks_500bp_remove_convert.bed"

	# pile up tags and output in cdt format.
	java -jar ./ScriptManager-v0.14.jar read-analysis tag-pileup "${folder}${ID}_chexmix_filtered_peaks_500bp_remove_convert.bed" "${folder}${ID}_filtered_sorted.bam" --cdt -M "${ID}_BED_${ID}_BAM"

	# calculate the scaling factor using NCIS.
	java -jar ./ScriptManager-v0.14.jar read-analysis scaling-factor -n -c ./masterNoTag_20180928/masterNoTag_20180928_sorted.bam "${folder}${ID}_filtered_sorted.bam"

	# obtain the scaling factor.
	Factor=$(tail -n 1 "./${ID}_filtered_sorted_ScalingFactors.out" | awk '{print $NF}')

	# apply the scaling factor to the cdt.
	java -jar ./ScriptManager-v0.14.jar read-analysis scale-matrix "${ID}_BED_${ID}_BAM_sense.cdt" -s $Factor -r 2 -l 3 -o "./SCALE_cdt/${ID}_BED_${ID}_BAM_sense_SCALE.cdt"
	java -jar ./ScriptManager-v0.14.jar read-analysis scale-matrix "${ID}_BED_${ID}_BAM_anti.cdt" -s $Factor -r 2 -l 3 -o "./SCALE_cdt/${ID}_BED_${ID}_BAM_anti_SCALE.cdt"
	
	# clean-up.
	mv "./${ID}_filtered_sorted_ScalingFactors.out" "./${ID}_filtered_sorted.NCIS_scaling-ccr.count" "./${ID}_filtered_sorted.NCIS_scaling-marginal.count" $folder
	rm "${folder}${ID}_chexmix_filtered_peaks_500bp_remove_temp.bed" "${folder}${ID}_chexmix_filtered_peaks_500bp_remove.bed" "${ID}_BED_${ID}_BAM_sense.cdt" "${ID}_BED_${ID}_BAM_anti.cdt" "${folder}${ID}_filtered.bam"
done