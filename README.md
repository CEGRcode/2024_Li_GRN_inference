# A unified framework for causal gene regulatory network inference grounded in orthogonal molecular evidence

## Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA
### Abstract
Gene regulatory networks (GRNs) govern gene expression, cellular differentiation, and stable transcriptional states. Yet inferring GRNs that integrate molecular regulatory mechanisms and reproduce transcriptional states as stable outcomes remains a central challenge. Here we present SETIA, a framework that infers GRNs whose explicit dynamical models reproduce transcriptional profiles as one or more stable states across conditions. Applied to RNA–seq data from wild–type and transcription factor knockout strains in *Saccharomyces cerevisiae*, SETIA infers GRNs that accurately reproduce held–out transcriptional states in cross–validation experiments. Incorporating TF–promoter binding and protein–protein interaction priors, SETIA yields GRNs ranging from mechanistically grounded architectures to flexible models that capture indirect regulatory influences. SETIA reveals that gene expression organizes into discrete stable states that represent distinct transcriptional programs, all emerging as stable attractors of a single underlying GRN whose dynamics are predominantly explained by TF–DNA binding and protein–protein interactions from orthogonal molecular evidence.


### Repository Structure

```text
.
├── SETIA/                                # Core GRN inference and simulation framework
├── TF_DNA_Binding_Network/               # TF–DNA regulatory networks
├── Protein_protein_interaction_network/  # Colocalization-based protein-protein interaction networks
├── human/                                # Human application (BMMC + GTEx): SETIA input building and Figure 3
├── GRN_simulator_website/                # Web interface
└── README.md
```
The repository grew yeast-first, so the yeast analysis lives in the component directories above (SETIA/, TF_DNA_Binding_Network/, Protein_protein_interaction_network/). The human work was added later as a self-contained application under human/, which reuses the SETIA/ engine.
