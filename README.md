# SETIA: causal gene regulatory network inference that discovers transcription-factor determinants of cell identity


## Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA
### Abstract
Cell types correspond to stable attractors of a gene regulatory network (GRN), yet this picture has lacked an executable dynamical model to make it concrete.  Data-driven methods infer static regulatory structure without runnable dynamics that converge to cell-type attractors, and existing dynamical models have largely been small, hand-built, and binary. Here we develop SETIA to infer GRNs whose ordinary-differential-equation dynamics reproduce cell-type expression profiles as stable attractors, integrating transcription-factor (TF)–DNA and TF–TF interaction evidence as priors. SETIA captures the multiple stable gene expression states in yeast and predicts those of TF perturbations held out from training, revealing regulatory logic that governs the unseen profiles. In human hematopoiesis, SETIA simulates TF knockout and forced expression to distinguish the sufficiency of a TF to induce a cell type from its necessity to maintain one, classifying TFs into drivers, stabilizers, master regulators, and bystanders. These classifications recover annotated hematopoietic regulator roles at high precision. Overall, SETIA turns the attractor picture of cell types into a predictive, molecularly grounded GRN model, applicable from yeast to human.



### Repository Structure

```text
.
├── SETIA/                                # Core GRN inference and simulation framework
├── TF_DNA_Binding_Network/               # TF–DNA regulatory networks
├── Protein_protein_interaction_network/  # Colocalization-based protein-protein interaction networks
├── Human/                                # Human application (BMMC + GTEx): SETIA input building and Figure 3
├── GRN_simulator_website/                # Web interface
└── README.md
```
The repository grew yeast-first, so the yeast analysis lives in the component directories above (SETIA/, TF_DNA_Binding_Network/, Protein_protein_interaction_network/). The human work was added later as a self-contained application under human/, which reuses the SETIA/ engine.
