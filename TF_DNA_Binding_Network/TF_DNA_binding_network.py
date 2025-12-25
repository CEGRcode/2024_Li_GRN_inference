import numpy as np
import pandas as pd
import math
import os
import json
from datetime import datetime

# obtain gene index from genome annotation
inputfile = open('./result/Sc_genome_annotations.txt', 'r')
GeneNames = []
for line in inputfile:
    if line.split('\t')[5] == '':
        GeneNames.append(line.split('\t')[4])
    else:
        GeneNames.append(line.split('\t')[5])
inputfile.close()
GeneNames = GeneNames[1:]

# Get the TFs in the NDR/NFR for each Tandem gene
inputfile = open('./result/Tandem_Genes_and_Stuff_inbetween.txt', 'r')

TF_Gene_Matrix = {}
for line in inputfile:
    if line.split('\t')[5] == '':
        genename = line.split('\t')[4]
    else:
        genename = line.split('\t')[5]
    temp_TFs = []
    for each in line.split('\t')[8:-1]:
        temp_TFs.append(each.split('_')[0])
    temp_TFs = set(temp_TFs)
    temp_vector = np.zeros(shape=np.array(GeneNames).shape)
    for each in temp_TFs:
        if each.upper() == 'NUC':
            continue
        else:
            temp_vector[GeneNames.index(each.upper())] = 1
    if genename in TF_Gene_Matrix:
        raise Exception('Existed gene!')
    else:
        TF_Gene_Matrix[genename] = temp_vector
inputfile.close()

#Get the TFs in the NDR/NFR for each pair of H-H gene
inputfile = open('./result/Divergent_Genes_and_Stuff_inbetween.txt', 'r')
for line in inputfile:
    temp_vector_left = np.zeros(shape=np.array(GeneNames).shape)
    temp_vector_right = np.zeros(shape=np.array(GeneNames).shape)
    left_ = int(float(line.split('\t')[3]))
    right_ = int(float(line.split('\t')[9]))
    insulator_ = []
    
    if abs(right_-left_) <= 300:
        # the NDR in between is short, consider the genes as co-regulated.
        pass
    else:
        for each in line.split('\t')[10:-1]:
            # judge the rap1, reb1, and abf1: insulator in the middle or repressor close to a gene.
            if (each.split('_')[0] in ['Rap1', 'Reb1', 'Abf1']) and abs(int(each.split('_')[-1])-0.5*(left_+right_)) <= 0.15*abs(right_-left_):
                insulator_.append(int(each.split('_')[-1]))
            elif each.split('_')[0] in ['Rap1', 'Reb1', 'Abf1']:
                if abs(int(each.split('_')[-1])-left_) < abs(int(each.split('_')[-1])-right_):
                    temp_vector_left[GeneNames.index(each.split('_')[0].upper())] = 1
                else:
                    temp_vector_right[GeneNames.index(each.split('_')[0].upper())] = 1
            else:
                pass
            
    for each in line.split('\t')[10:-1]:
        if each.split('_')[0] not in ['Rap1', 'Reb1', 'Abf1', 'nuc']:
            if len(insulator_) == 0:
                temp_vector_left[GeneNames.index(each.split('_')[0].upper())] = 1
                temp_vector_right[GeneNames.index(each.split('_')[0].upper())] = 1
            else:
                if abs(int(each.split('_')[-1])-left_) < abs(min(insulator_)-left_):
                    temp_vector_left[GeneNames.index(each.split('_')[0].upper())] = 1
                elif abs(int(each.split('_')[-1])-right_) < abs(max(insulator_)-right_):
                    temp_vector_right[GeneNames.index(each.split('_')[0].upper())] = 1
                else:
                    pass
                    #print('buried in insulator.')
        else:
            pass
    # for the left gene:
    genename = line.split('\t')[0]
    if genename in TF_Gene_Matrix:
        TF_Gene_Matrix[genename] = TF_Gene_Matrix[genename] + temp_vector_left
    else:
        TF_Gene_Matrix[genename] = temp_vector_left
    #for the right gene:
    genename = line.split('\t')[5]
    if genename in TF_Gene_Matrix:
        TF_Gene_Matrix[genename] = TF_Gene_Matrix[genename] + temp_vector_right
    else:
        TF_Gene_Matrix[genename] = temp_vector_right
    
inputfile.close()

for each_key in TF_Gene_Matrix:
    for each_i in range(0, len(TF_Gene_Matrix[each_key])):
        if TF_Gene_Matrix[each_key][each_i] == 0:
            continue
        else:
            TF_Gene_Matrix[each_key][each_i] = 1

# get the ssTFs to build the TF-DNA binding network
# ---- read TF names safely ----
ssTF_names = []
with open('./data/All_gene_ID.txt', 'r') as infile:
    for line in infile:
        token = line.split()
        if token:
            ssTF_names.append(token[0])

# quick helpers: detect TF_Gene_Matrix type
is_df = isinstance(TF_Gene_Matrix, pd.DataFrame)

# build a fast lookup for gene name -> index (for sequence-like rows)
gene_index = {g: i for i, g in enumerate(GeneNames)}

# validate TF names: keep only those present in the matrix and also present in GeneNames (for columns)
valid_tfs = []
missing_tfs = []
missing_genes = set()
for name in ssTF_names:
    ok = False
    if is_df:
        # DataFrame: check both index and that we can find at least one column for safety
        if name in TF_Gene_Matrix.index:
            ok = True
    else:
        # assume dict-like
        try:
            _ = TF_Gene_Matrix[name]
            ok = True
        except Exception:
            ok = False

    if ok:
        valid_tfs.append(name)
    else:
        missing_tfs.append(name)

if missing_tfs:
    print(f"Warning: the following TFs from input file were not found in TF_Gene_Matrix and will be skipped: {missing_tfs}")

# If no valid TFs, exit gracefully
if not valid_tfs:
    raise RuntimeError("No valid TF names found in TF_Gene_Matrix. Aborting.")

# ---- build OutMatrix (rows = valid_tfs, cols = valid_tfs) ----
OutMatrix_list = []
for row_tf in valid_tfs:
    row_vals = []
    for col_tf in valid_tfs:
        val = 0
        try:
            if is_df:
                # prefer label access by column (column = gene name)
                # if TF_Gene_Matrix is TF x gene (index: TF, columns: gene names)
                # then use .at[row_tf, col_tf] only if column name equals gene name
                if col_tf in TF_Gene_Matrix.columns:
                    val = TF_Gene_Matrix.at[row_tf, col_tf]
                else:
                    # try using gene index if DataFrame columns are indices/positions
                    if col_tf in gene_index:
                        # attempt iloc by gene index
                        idx = gene_index[col_tf]
                        val = TF_Gene_Matrix.iloc[TF_Gene_Matrix.index.get_loc(row_tf), idx]
                    else:
                        # missing gene column
                        val = 0
                        missing_genes.add(col_tf)
            else:
                # dict-like case
                row = TF_Gene_Matrix.get(row_tf)
                if row is None:
                    val = 0
                else:
                    # row could be a dict (gene->value) or sequence (indexed by gene index)
                    if isinstance(row, dict):
                        val = row.get(col_tf, 0)
                        if col_tf not in row:
                            missing_genes.add(col_tf)
                    else:
                        # sequence-like: use gene_index
                        if col_tf in gene_index:
                            idx = gene_index[col_tf]
                            try:
                                val = row[idx]
                            except Exception:
                                val = 0
                                missing_genes.add(col_tf)
                        else:
                            val = 0
                            missing_genes.add(col_tf)
        except Exception:
            # protect against any unexpected error per cell
            val = 0
            missing_genes.add(col_tf)

        # coerce to int/0/1 if necessary
        try:
            # If val is array-like, try to extract scalar
            if isinstance(val, (np.ndarray, list)) and len(val) == 1:
                val = val[0]
            val = int(val)
        except Exception:
            # fallback: treat as zero
            val = 0

        row_vals.append(val)
    OutMatrix_list.append(row_vals)

# create DataFrame (use valid_tfs order)
OutMatrix = pd.DataFrame(np.array(OutMatrix_list).T, index=valid_tfs, columns=valid_tfs)

if missing_genes:
    print(f"Note: some column gene names were missing when reading rows: {sorted(missing_genes)} (treated as 0)")

# ---- build JSON structure ----
nodes = []
for name in valid_tfs:
    nodes.append({
        "id": name,
        "label": name,
    })

edges = []
# find entries equal to 1 (or nonzero)
rows, cols = np.where(OutMatrix.values == 1)
for r, c in zip(rows, cols):
    src = OutMatrix.index[r]
    tgt = OutMatrix.columns[c]
    edges.append({
        "source": src,
        "target": tgt,
        "label": "",
        "style": ["solid", "triangle"]
    })

json_obj = {"nodes": nodes, "edges": edges}

# ---- write to file (create dir if needed) ----
out_dir = './result'
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, f"TF_DNA_binding_network_{datetime.today().strftime('%d%m%y')}.json")

with open(out_path, 'w', encoding='utf-8') as outfile:
    json.dump(json_obj, outfile, indent=2)

print(f"Wrote {len(nodes)} nodes and {len(edges)} edges to {out_path}")
