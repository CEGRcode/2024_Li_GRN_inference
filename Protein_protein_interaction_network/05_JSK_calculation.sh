#!/bin/bash
# ---------------------- Calculate the Jensen-Shannon divergence and K ratio for each factor ------------
# ---------------------- rl884@cornell.edu, 25/04/25 ----------------------------------------------------

# Define the target directory
target_dir="./SCALE_cdt"

# Loop through all .zip files in the directory
for zip_file in "$target_dir"/*_zipped_BED_BAM_sense.zip; do
  # Check if there are any .zip files
  if [[ -f "$zip_file" ]]; then
    BED_ID=${zip_file##*/}
    BED_ID=${BED_ID%%_*}

    # Check if BED_ID is already in JS_calculation.log
    if grep -q "^$BED_ID$" ./JS_calculation.log; then
      echo "$BED_ID already processed. Skipping..."
      continue
    fi

    # Unzip the file
    unzip -o "$zip_file"

    # Update the timestamp of extracted files
    touch ./SCALE_cdt/"$BED_ID"*

    # Run JS and K ratio calculation
    python JSK_calculation.py "$BED_ID"

    # Clean up
    rm ./SCALE_cdt/"$BED_ID"*.cdt

    # Log the processed BED_ID
    echo "$BED_ID" >> ./JS_calculation.log
  else
    echo "No zip files found in $target_dir."
    break
  fi
done
