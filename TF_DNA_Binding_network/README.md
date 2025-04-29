# TF-DNA binding network inference from ChIP-exo data on YEP

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

## Directions
To recreate the TF-DNA binding network for this manuscript, please execute the shell script `Build_TF_DNA_network`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n TF_DNA -c bioconda -c conda-forge numpy pandas os math datetime
```

## Table of Contents

### Dissect_full_spreadsheet
Preprocesses the supplementary material spreadsheet, 41586_2021_3314_MOESM3_ESM, from (Rossi et al. Nature, 2021).

### ssTFs_common_names
Contains the sequence-specific transcription factors to be included in the TF-DNA binding network.

### 2021-Rossi_Nature-master
Contains the processed data downloaded from (Rossi et al. Nature, 2021).
