# A causal gene regulatory network grounded in orthogonal molecular evidence
## Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA
### Abstract
The coordinated activity of genes within an organism's genome forms a gene regulatory network (GRN) that governs gene expression, cellular differentiation, and maintenance of stable transcriptional states. Transcriptional profiles across cell types, conditions, and genetic perturbations reflect underlying dynamics of GRNs, yet inferring GRNs that both integrate molecular regulatory mechanisms and reproduce transcriptional states as stable outcomes remains a central challenge. We present SETIA, a framework that identifies multiple discrete gene expression states from RNA–seq data and infers GRNs whose explicit dynamical models reproduce transcriptional profiles as stable states across conditions. Applied to RNA–seq data from wild–type and transcription factor knockout strains in *Saccharomyces cerevisiae*, SETIA infers GRNs that accurately reproduce observed transcriptional profiles and generalize to held–out transcriptional states in cross–validation. Incorporating TF–promoter binding and protein–protein interaction priors at differing strengths, SETIA yields GRNs ranging from mechanistically grounded architectures to flexible models that capture indirect regulatory influences. SETIA provides a systems–level framework for GRN inference that explains and predicts transcriptional states while maintaining fidelity to its molecular regulatory structure.

### Repository Structure

```text
.
├── SETIA/                                # Core GRN inference and simulation framework
├── TF_DNA_Binding_Network/               # TF–DNA regulatory networks
├── Protein_protein_interaction_network/  # PPC-based regulation
├── GRN_simulator_website/                # Web interface
└── README.md
```
