# Gene regulatory network inference using IAGREE

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

## Directions
To recreate the TF-DNA binding network for this manuscript, please execute the shell script `Build_TF_DNA_network`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n TF_DNA -c bioconda -c conda-forge numpy pandas os json math datetime
```

## Table of Contents

### salmon_tximport_tmm_normalize
Performs the TMM normalization on the RNA-seq samples aligned by Salmon. It produces the normalized read counts in GRN_Sc_TMM_normalized_CPM.txt for downstream analysis.
