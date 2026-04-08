# Genome-wide colocalization based protein-protein interaction network inference from ChIP-exo data on YEP

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA

### Pipeline 1: Motif-based PPC network  
1. `00_Download_and_touch.sh` → downloads and initializes required data  
2. `01_Get_MEME_report.sh` → generates motif discovery results  
3. `02_Get_YEP_motifs_and_binding_sites.sh` → identifies motif-supported binding sites  
4. `030_Tag_pileup_and_NCIS_1.sh` → computes tag pileup and normalization (part 1)  
5. `040_Tag_pileup_and_NCIS_2.sh` → computes tag pileup and normalization (part 2)  
6. `05_JSK_calculation.sh` → calculates Jensen–Shannon-based similarity  
7. `06_Secondary_sort.sh` → refines interaction ranking  
8. `070_Build_PPI_network.sh` → constructs the final motif-based PPI network  

### Pipeline 2: PPC network using all bound locations  
1. `00_Download_and_touch.sh` → downloads and initializes required data  
2. `01_Get_MEME_report.sh` → generates motif discovery results  
3. `02_Get_YEP_motifs_and_binding_sites.sh` → identifies all binding sites  
4. `031_Tag_pileup_and_NCIS_1.sh` → computes tag pileup and normalization (part 1)  
5. `041_Tag_pileup_and_NCIS_2.sh` → computes tag pileup and normalization (part 2)  
6. `05_JSK_calculation.sh` → calculates Jensen–Shannon-based similarity  
7. `06_Secondary_sort.sh` → refines interaction ranking  
8. `071_Build_PPI_network.sh` → constructs the final PPI network  

### Notes
- Scripts must be executed in the specified order.  
- The two pipelines differ in how binding sites are defined (motif-based vs. all bound locations).  
- This workflow is computationally intensive and may require parallel computing (e.g., HPC clusters).

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n PPI -c bioconda -c conda-forge \
    python=3.10 \
    bedtools=2.31 \
    samtools=1.19 \
    bwa=0.7 \
    sra-tools=3.0 \
    pybigwig=0.3 \
    numpy=1.26 \
    pandas=2.2 \
    scipy=1.12 \
    scikit-learn=1.4 \
    matplotlib=3.8 \
    opencv=4.9 \
    wget=1.21
```

## Table of Contents

### 00_Download_and_touch
Downloads the ChIP-exo data from the yeast epigenome project (Rossi et al. Nature, 2021) and updates the timestamps.

### 01_Get_MEME_report
Collects the MEME reports for all factors to a directory.

### 02_Get_YEP_motifs_and_binding_sites
Aggregates all binding sites called by ChExMix and collects motifs identified by MEME for each factor into two separate directories, YEP_Motif_BED and YEP_ALL_BED.

### 030_Tag_pileup_and_NCIS_1
Selects the primary motif that covers the most peaks, if motifs exist, for each factor, does NCIS normalization against masterNoTag, and performs bulk read pileups on the motif instances.

### 031_Tag_pileup_and_NCIS_1
Does NCIS normalization against masterNoTag, and performs bulk read pileups on all binding sites identified by ChExMix.

### 040_Tag_pileup_and_NCIS_2
Performs bulk read pileups for a factor (BAM) on the motif instances of all other factors.

### 041_Tag_pileup_and_NCIS_2
Performs bulk read pileups for a factor (BAM) on the ChExMix binding sites of all other factors.

### 05_JSK_calculation
Calculates the Jensen-Shannon divergence and K ratio distributions for each factor.

### 06_Secondary_sort
Performs the secondary sorting to the colocalization factors and filters out the low-read artifacts.

### 4932.protein.info.v12.0 and 4932.protein.links.v12.0
Protein-protein interaction network downloaded from the STRING database.

### FeatureClass_Genes and SupplementaryData-Table4_Sample-Key_tabular
Common names, systematic IDs, classes of all genes in yeast and the ChIP-exo sample codes as obtained from (Rossi et al. Nature, 2021) supplementary material.

### ssTFs_common_names
The sequence-specific transcription factors selected to build the PPI in this study.

### YEP_GO
The gene ontology terms for genes in yeast.

### ScriptManager-v0.14
ScriptManager used in this study to make NCIS normalization and bulk read pileups.
