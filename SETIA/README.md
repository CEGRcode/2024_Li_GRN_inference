# Gene Regulatory Network Inference Using SETIA

### Ruihao Li<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1</sup>
<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA  
<sup>2</sup>Department of Computational Biology, Cornell University, Ithaca, New York, 14853, USA

## Overview
SETIA integrates multiple layers of molecular evidence to infer mechanistically interpretable and dynamically executable gene regulatory networks (GRNs).

Supported evidence layers include:

- RNA-seq transcriptional profiles across genotypes and perturbations
- ChIP-exo TF–DNA binding
- PRO-seq / CAGE-seq promoter activity
- Protein–protein interactions (PPI), including protein–protein colocalization (PPC)

These data are integrated to construct and constrain GRNs whose dynamics are modeled using ordinary differential equations (ODEs).

---

## Key Features
- Infers dynamical GRNs that reproduce stable transcriptional states
- Integrates multi-omics data as mechanistic structural constraints
- Supports combinatorial regulation
- Implements parallelized optimization for scalable GRN inference
- Provides an interactive web platform for simulation and visualization

---

## Pipeline to Reproduce GRNs from the Paper
1. `0_Download_and_touch.sh`  
   Downloads and initializes required data

2. `1_Preprocessing_dicrete_gene_expression_state_identification.sh`  
   Identifies discrete gene expression states

3. `2_GRN_inference.sh`  
   Performs GRN inference

4. `3_Final_GRN_construction.sh`  
   Filters and constructs the final GRN

5. `4_Scale_GRN_to_nonTF_targets.sh`  
   Extends the GRN to downstream non-TF targets

---

## Dependencies
We recommend initializing the following dedicated Anaconda environment for reproducibility:

```bash
conda create -n EvoAlg -c conda-forge -y     python=3.10     numpy=1.26     pandas=2.2     scipy=1.12     sympy=1.12     numexpr=2.10     h5py=3.11     matplotlib=3.8     scikit-learn=1.4     pip
conda activate EvoAlg
pip install flask tqdm seaborn joblib
```

---
## Repository Structure
A recommended repository organization for reproducing the SETIA workflow is:

```text
SETIA/
├── data/
│   ├── example_RNA_profiles.txt
│   ├── example_promoter_strength.txt
│   ├── example_gene_length.txt
│   ├── GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt
│   ├── GRN_ssTFs_Salmon_SteadyStates_2025_std.txt
│   ├── GRN_ssTFs_Sc_gene_length.txt
│   ├── GRN_ssTFs_Sc_promoter_strength.txt
│   ├── GRN_ssTFs_Sc_TF_DNA.txt
│   ├── GRN_ssTFs_Sc_LG.txt
│   └── ...
├── result/
│   ├── GRN_filtered_Sc_GRN_final_2.json
│   ├── GRN_filtered_full_Sc_GRN_final.txt
│   ├── *.pkl
│   └── GMM_figures/
├── helper/
│   ├── path_setup.py
│   └── utility scripts
├── dynamics_module/
│   ├── dynamics.py
│   └── ODE simulation modules
├── EGRN_Multi_Genalg_Combinatorial_2025.py
├── GRN_input_acquisition_v2.py
├── GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py
├── GRN_Dynamic_Simulator_Combinatorial_Global_multistate_2025.py
├── Export_GRN_from_ResultFlow_Component.py
├── salmon_tximport_tmm_normalize.R
└── README.md
```
---
## Outputs
The inferred GRN is written to the `./result` directory.

Key outputs include:
- `GRN_filtered_Sc_GRN_final_2.json`: final GRN structure in json format
- `GRN_filtered_full_Sc_GRN_final.txt`: fitted model parameters in text
- `*.pkl`: intermediate and final serialized GRN instances
- `sys_cache_*`: cached searched or simulated GRN configurations

---

## Minimal Example
This example demonstrates how to run SETIA on a small test dataset. Because the code is optimized for HPC environments, local execution may take longer (~10 minutes).

The inferred GRN will be written to:
`./result/GRN_filtered_Sc_GRN_final_raw.json`

This output can be uploaded to the GRN simulator website for interactive exploration:
https://grn.cac.cornell.edu:5000/

```bash
cd SETIA
conda activate EvoAlg

python EGRN_Multi_Genalg_Combinatorial_2025.py     -r data/example_RNA_profiles.txt     -n 3000     -i 10     -p 0     -t data/example_promoter_strength.txt     -l data/example_gene_length.txt     -o "Minimal_example"     -e 42     -f 0     -k 0

python Export_GRN_from_ResultFlow_Component.py
conda deactivate
```

---

## Usage with Custom Data

### Step 1: Prepare Gene Metadata
Create a gene metadata file describing all genes to be modeled in the GRN.

Required format:
```text
TF_NAME    SGD_ID    GENE_LENGTH
```

Example:
```text
ABF1    YKL112W    2196
AFT1    YGL071W    2073
MET4    YNL103W    2019
```

---

### Step 2: Prepare Normalized RNA-seq Input
RNA-seq data are required.

We recommend TMM normalization, and an example workflow is provided in:
`salmon_tximport_tmm_normalize.R`

Other normalization methods are also supported, provided the final input format remains compatible.

Expression matrix (`-r`) requirements:
- first column: `samples`
- remaining columns: SGD systematic gene IDs
- rows: individual RNA-seq samples
- values: normalized expression values

Example:
```text
samples    YKL112W    YGL071W    YPL202C
WT_batch_1_rep_1    120.5    45.2    18.9
WT_batch_2_rep_1    118.3    44.7    19.1
AZF1_batch_1_rep_1  12.4     39.8    5.2
```

Requirements:
- the `samples` header is required
- sample names must use underscore-delimited formatting
- the text before the first underscore is treated as the genotype label

Recommended naming convention:
```text
GENOTYPE_batch_batch#_rep_replicate#
```

---

### Step 3: Run Input Preprocessing
Use `GRN_input_acquisition_v2.py` to convert the gene metadata file, normalized RNA-seq matrix, and optional prior information into SETIA-compatible preprocessing files.

Minimal required usage:

```bash
python GRN_input_acquisition_v2.py   -g ./data/ssTFs_MATa_Spots_76.txt   -r ./data/GRN_Sc_TMM_normalized_CPM.txt
```

Full usage with optional prior information:

```bash
python GRN_input_acquisition_v2.py   -p 0.01   -g ./data/ssTFs_MATa_Spots_76.txt   -r ./data/GRN_Sc_TMM_normalized_CPM.txt   -t ./data/Rossi_Ruihao_TF_DNA_union_motif_based.json   -c ./data/PPI_network_Cutoff_0_STRING_overlapping_motif_sites_0_2025_union.json   -a ./data/Sc_genome_annotations.txt   -b ./data/YEP_best_rep.txt
```

Command-line arguments:

**Required**
- `-g`: gene metadata file
- `-r`: normalized RNA-seq expression matrix

**Optional**
- `-p`: p-value cutoff for discrete state identification (default: `0.01`)
- `-t`: TF–DNA prior in JSON format
- `-c`: protein–protein colocalization prior in JSON format
- `-a`: genome annotation file
- `-b`: YEP replicate mapping file

This script writes preprocessing outputs primarily to `./data` and diagnostic summaries to `./result`.

**Core outputs used by SETIA**
- `./data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt`  
  Discretized steady-state expression matrix used as the transcriptional-state input for GRN inference
- `./data/GRN_ssTFs_Salmon_SteadyStates_2025_std.txt`  
  Per-gene standard deviations associated with the discretized steady states
- `./data/GRN_ssTFs_Sc_gene_length.txt`  
  Gene lengths for the modeled genes
- `./data/GRN_ssTFs_row_names.txt`  
  Row labels corresponding to genotype or condition names
- `./data/GRN_ssTFs_column_names.txt`  
  Column labels corresponding to the modeled genes

**Optional prior-derived outputs**
- `./data/GRN_ssTFs_Sc_TF_DNA.txt`  
  TF–DNA prior adjacency matrix generated when `-t` is provided
- `./data/GRN_ssTFs_Sc_LG.txt`  
  Local-group / complex-structure file derived from promoter binding and colocalization information when the corresponding optional prior files are provided
- `./data/GRN_ssTFs_Sc_promoter_strength.txt`  
  Promoter strength matrix generated when the required promoter-annotation resources are available
- `./data/GRN_ssTFs_Sc_TF_DNA_TPM_union.txt`  
  Union of TF–DNA prior and TPM-derived relationships, generated when `-t` is provided
- `./data/GRN_ssTFs_Sc_initial_condition.txt`  
  Initial GRN edge-state string derived from TF–DNA and expression similarity, generated when `-t` is provided

**Additional preprocessing summaries**
- `./data/GRN_ssTFs_Salmon_SteadyStates_2025.txt`  
  Average / grouped steady-state expression values before discretization
- `./result/Steady_state_count.txt`  
  Summary of the number of inferred discrete states per gene
- `./result/GMM_figures/AIC/`  
  Diagnostic plots for discrete-state identification

If optional inputs such as `-t`, `-c`, `-a`, or `-b` are not provided, the corresponding prior-dependent preprocessing steps are skipped automatically.

---

### Step 4: Run SETIA Inference
Run the main evolutionary GRN inference using the processed inputs from Step 3.

General example command:

```bash
python EGRN_Multi_Genalg_Combinatorial_2025.py     -r ./data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt     -n 3000     -i 100     -p 0     -t ./data/GRN_ssTFs_Sc_promoter_strength.txt     -l ./data/GRN_ssTFs_Sc_gene_length.txt     -o "SETIA_run"     -e 42     -f 0     -k 1
```

Command-line arguments:

**Required**
- `-r`: discretized transcriptional steady-state matrix
- `-t`: promoter strength matrix
- `-l`: gene length file

**Optional**
- `-c`: structural prior adjacency matrix
- `-g`: logic-gate / local-group file
- `-b`: per-gene steady-state standard deviation file
- `-d`: initial GRN edge-state string
- `-n`: simulation time span used in `solve_ivp`
- `-i`: number of evolutionary optimization iterations
- `-p`: perturbation magnitude applied to initial states
- `-o`: output file prefix
- `-e`: random seed
- `-f`: constraint mode
- `-k`: whether to serialize intermediate GRN instances as pickle files

The files `GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt`, `GRN_ssTFs_Sc_promoter_strength.txt`, and `GRN_ssTFs_Sc_gene_length.txt` generated in Step 3 can be used directly here.

#### GRN configurations used in this study
We evaluated multiple GRN configurations in the manuscript by changing the structural prior supplied to `-c` and the constraint mode set by `-f`.

**GRN A: unconstrained GRN**
```bash
python EGRN_Multi_Genalg_Combinatorial_2025.py     -c data/GRN_ssTFs_Sc_TF_DNA.txt     -g data/GRN_ssTFs_Sc_LG.txt     -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt     -b data/GRN_ssTFs_Salmon_SteadyStates_2025_std.txt     -n 3000     -i 200     -p 0     -t data/GRN_ssTFs_Sc_promoter_strength.txt     -l data/GRN_ssTFs_Sc_gene_length.txt     -o "${output}"     -e "${seed}"     -d data/GRN_ssTFs_Sc_initial_condition.txt     -f 0     -k 0
```

**GRN B**
- `-f 1`
- `-c data/GRN_ssTFs_Sc_TF_DNA_TPM_union.txt`

**GRN C**
- `-f 1`
- `-c data/GRN_ssTFs_Sc_TF_DNA.txt`

**GRN D**
- `-f 1`
- `-c data/GRN_ssTFs_Sc_TF_DNA_Motif.txt`

Unless otherwise specified, the remaining command-line arguments are the same as those shown for GRN A.

---

### Step 5: Evaluate GRN Dynamical Performance
Use `GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py` to evaluate an inferred GRN for its ability to reproduce the target transcriptional profiles as dynamical stable states.

Example command:

```bash
python GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py     -j ./result/GRN_filtered_full_Sc_GRN_final_raw.txt     -r ./data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt     -l ./data/GRN_ssTFs_Sc_gene_length.txt     -t ./data/GRN_ssTFs_Sc_promoter_strength.txt     -o GRN_dynamical_performance     -m 0
```

Command-line arguments:

**Required**
- `-j`: inferred GRN parameter file to be evaluated
- `-r`: discretized transcriptional steady-state matrix used as the target profiles
- `-l`: gene length file
- `-t`: promoter strength matrix
- `-o`: output name prefix

**Optional**
- `-m`: random-mode switch for random-network evaluation (`0` for the inferred GRN; nonzero / true-like values enable random evaluation mode)
- `-p`: perturbation magnitude applied to the initial transcriptional states before simulation
- `-n`: simulation time span used during numerical integration
- `-i`: number of simulation iterations
- `-e`: random seed
- `-c`: gene index to focus on during detailed output reporting
- `-f`: whether to apply refinement results from `Refinement_for_gene_*.txt`

Inputs:
- The GRN file supplied with `-j` is typically the fitted GRN parameter file generated by the SETIA inference workflow.
- The files supplied with `-r`, `-l`, and `-t` are compatible with the preprocessing outputs generated in Step 3.

Main output:
- `./GRN_dynamic_local_search_result_<output_name>.txt`

This output file reports, for each transcriptional profile:
- the simulated initial state
- the simulated final state
- whether the final state is classified as a point attractor
- the attractor distance between the simulated final state and the target transcriptional profile

At the end of the file, the script also reports:
- the mean attractor distance across all evaluated profiles (`overall performance`)
- the sum of attractor distances across all evaluated profiles

Interpretation:
Lower attractor distances indicate better dynamical agreement between the inferred GRN and the target transcriptional profiles. A strong-performing GRN should drive perturbed initial states back toward the expected transcriptional stable states while maintaining point-attractor behavior for the evaluated profiles. The mean per-gene normalized L1 distances shown in Figures 3b and 6c of the paper are computed from these initial and final point-attractor profiles.

---

## File Descriptions

### `salmon_tximport_tmm_normalize.R`
Performs TMM normalization on Salmon-aligned RNA-seq samples.

### `EGRN_Multi_Genalg_Combinatorial_2025.py`
Main SETIA executable implementing evolutionary optimization of ODE-based GRNs.

### `GRN_input_acquisition_v2.py`
Preprocesses RNA-seq and optional prior data into SETIA-compatible inputs.

### `GRN_Dynamic_Simulator_Combinatorial_remove_dispensible_edges_2025.py`
Removes dispensable GRN edges that do not affect stable-state recovery.

### `Export_final_GRN_TF_DNA.py`
Exports the final GRN and compares it against the TF–DNA prior.

### `EGRN_Multi_non_TF_2025.py`
Extends the core ssTF GRN to downstream non-TF targets.

### `Export_GRN_from_ResultFlow_Component.py`
Recovers and exports GRNs from intermediate saved optimization states.

### `GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py`
Searches for stable states of a given gene regulatory network (GRN) using a local exploration strategy. The simulator initializes the system from a specified set of initial transcriptional states and evolves the GRN dynamics forward in time to determine whether each trajectory converges to a stable or unstable final state. This mode is intended to assess attractor stability and basin structure near biologically relevant initial conditions, such as the transcriptional profiles of stable cells.

### `GRN_Dynamic_Simulator_Combinatorial_Global_multistate_2025.py`
Searches for stable states of a given GRN using a global exploration strategy. The simulator samples a large, evenly distributed set of initial states across the full expression space, where each dimension corresponds to a gene’s expression level, and evolves the GRN dynamics from each starting point. This mode is designed to systematically map the global attractor landscape of the GRN and identify all accessible stable states, and therefore will be computationally expensive.
