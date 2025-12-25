# Gene regulatory network inference using IAGREE

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

## Directions
To recreate the gene regulatory network for this manuscript, please execute the shell scripts in order: `0_Download_and_touch.sh`, `1_Preprocessing_dicrete_gene_expression_state_identification.sh`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n EvoAlg -c conda-forge -y python=3.10 numpy pandas scipy matplotlib scikit-learn pip
```

## Table of Contents

### salmon_tximport_tmm_normalize
Performs the TMM normalization on the RNA-seq samples aligned by Salmon. It produces the normalized read counts in GRN_Sc_TMM_normalized_CPM.txt for downstream analysis.
