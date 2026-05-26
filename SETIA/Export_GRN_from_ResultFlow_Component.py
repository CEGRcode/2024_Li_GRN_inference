import numpy as np
import os
import sys
import json
import path_setup
from utility_functions import *
from pathlib import Path

directory = './result/'
sys_output_name = 'yeast_f1_GRN_final'

# 读GRN_components
GRN_Components_Dict = {}
for filename in sorted(os.listdir(directory)):
    if not filename.endswith('_GRN_components.txt'):
        continue
    seed_key = filename.split('_')[0]  # seed0, seed1, ...
    with open(os.path.join(directory, filename), 'r') as f:
        for line in f.readlines():
            ls = line.replace('\n', '').split('\t')
            if ls[0] == 'GeneIndex':
                continue
            if seed_key not in GRN_Components_Dict:
                GRN_Components_Dict[seed_key] = {}
            GRN_Components_Dict[seed_key][ls[0]] = ls[1:]

print(f"Read seeds: {sorted(GRN_Components_Dict.keys())}")

# 读ResultFlow
Result_Flow_Dict = {}
for filename in sorted(os.listdir(directory)):
    if not filename.endswith('_ResultFlow.txt'):
        continue
    seed_key = filename.split('_')[0]
    with open(os.path.join(directory, filename), 'r') as f:
        lines = f.readlines()
        if lines:
            Result_Flow_Dict[seed_key] = list(map(float, lines[-1].strip().split()))

print(f"ResultFlow seeds: {sorted(Result_Flow_Dict.keys())}")

# gene list
indexes_of_diff_genes = sorted(list(map(int, GRN_Components_Dict[next(iter(GRN_Components_Dict))].keys())))
print(f"diff_genes数: {len(indexes_of_diff_genes)}")

# 写indexes_of_diff_gene.txt
with open('./result/indexes_of_diff_gene.txt', 'w') as f:
    f.write('\t'.join(map(str, indexes_of_diff_genes)))

# gene名称
file_path = Path('./data/GRN_ssTFs_column_names.txt')
Order_of_genes = []
if file_path.exists():
    with file_path.open() as f:
        for line in f:
            Order_of_genes = line.split()
if not Order_of_genes:
    Order_of_genes = [chr(i) for i in range(ord('A'), ord('Z') + 1)]
print(f"Order_of_genes数: {len(Order_of_genes)}")

# TF_DNA_net（用于JSON可视化）
def get_edges_from_json(json_file):
    edge_dict = {}
    if not os.path.exists(json_file):
        return edge_dict
    with open(json_file, 'r') as file:
        data = json.load(file)
    edges = data.get('edges', [])
    for edge in edges:
        source = edge['source']
        target = edge['target']
        if source in edge_dict:
            edge_dict[source].append(target)
        else:
            edge_dict[source] = [target]
    return edge_dict

TF_DNA_net = get_edges_from_json('./data/Rossi_Ruihao_TF_DNA_union_all.json')

# 互补：每个基因取最优seed
seeds = sorted(Result_Flow_Dict.keys())
n_genes = len(indexes_of_diff_genes)

min_keys = []
min_values = []
for j in range(n_genes):
    best_seed = min(seeds, key=lambda s: Result_Flow_Dict[s][j])
    min_keys.append(best_seed)
    min_values.append(Result_Flow_Dict[best_seed][j])

# fitness>=0.5的基因去掉
TG_to_remove = [indexes_of_diff_genes[i] for i, v in enumerate(min_values) if v >= 0.5]
print(f"TG_to_remove: {TG_to_remove}")

# 组装各列表
AM_List = []; LG_List = []; f0_c_List = []; Sigmoid_k_List = []
Leakage_List = []; Degradation_mRNA_List = []; Transcription_Rate_List = []
Transcriptional_Threshold_List = []
kept_diff_genes = []

with open('./result/Preliminary_training_result.txt', 'w') as outfile:
    for i, (key, value) in enumerate(zip(min_keys, min_values)):
        gi = indexes_of_diff_genes[i]
        outfile.write(f"{gi}\t{value}\n")
        print(f"Target Gene {gi} ({Order_of_genes[gi] if gi < len(Order_of_genes) else '?'}): seed={key}, fitness={value:.4f}")
        if gi not in TG_to_remove:
            am = GRN_Components_Dict[str(key)][str(gi)]
            AM_List.append(am[0])
            LG_List.append(am[1])
            f0_c_List.append(list(map(float, am[2].split(','))))
            Sigmoid_k_List.append(float(am[3]))
            Leakage_List.append(float(am[4]))
            Degradation_mRNA_List.append(float(am[5]))
            Transcription_Rate_List.append(float(am[6]))
            Transcriptional_Threshold_List.append(eval(am[7]))
            kept_diff_genes.append(gi)

indexes_of_diff_genes = kept_diff_genes
print(f"diff_genes: {len(indexes_of_diff_genes)}个")

TotalNumberOfGenes = len(AM_List[0])
Order_of_genes = Order_of_genes[:TotalNumberOfGenes]

# 组装Final_AM
Final_AM = ''
for j in range(TotalNumberOfGenes):
    for k in range(TotalNumberOfGenes):
        if k in indexes_of_diff_genes:
            Final_AM += str(AM_List[indexes_of_diff_genes.index(k)][j])
        else:
            Final_AM += '0'

# 组装Final_LG等
Final_LG = []; Final_Sigmoid_k = []; Final_mRNA_DegradationRate = []
Final_Leakage = []; Final_f0_c = []; Final_Transcription_Rate = []
Final_Transcriptional_Threshold = []

for j in range(TotalNumberOfGenes):
    if j in indexes_of_diff_genes:
        p = indexes_of_diff_genes.index(j)
        Final_LG.append(list(map(int, LG_List[p].split(','))))
        Final_Sigmoid_k.append(Sigmoid_k_List[p])
        Final_mRNA_DegradationRate.append(Degradation_mRNA_List[p])
        Final_Leakage.append(Leakage_List[p])
        Final_f0_c.append(f0_c_List[p])
        Final_Transcription_Rate.append(Transcription_Rate_List[p])
        Final_Transcriptional_Threshold.append(Transcriptional_Threshold_List[p])
    else:
        Final_LG.append([unit for unit in range(TotalNumberOfGenes)])
        Final_Sigmoid_k.append(0)
        Final_mRNA_DegradationRate.append(0)
        Final_Leakage.append(0)
        Final_f0_c.append([0,0,0,0,0,0])
        Final_Transcription_Rate.append(0)
        Final_Transcriptional_Threshold.append([[0,0,0] for _ in range(TotalNumberOfGenes)])

# 保存brief和full
with open(f'./result/GRN_filtered_brief_{sys_output_name}_raw.txt', 'w') as f:
    f.write('{}\t{}\n'.format(Final_AM, Matrix2String01_LG_Expanded(Final_LG)))

with open(f'./result/GRN_filtered_full_{sys_output_name}_raw.txt', 'w') as f:
    f.write('Adjacency Matrix:\t{}\nLogic Gate:\t{}\nf0:\t{}\nHill Coefficient:\t{}\nTranscription Rate:\t{}\nmRNA Degradation Rate:\t{}\nLeakage Rate:\t{}\nTF Effective Threshold:\t{}\n\n'.format(
        Final_AM, Matrix2String01_LG_Expanded(Final_LG), Final_f0_c,
        Final_Sigmoid_k, Final_Transcription_Rate, Final_mRNA_DegradationRate,
        Final_Leakage, Final_Transcriptional_Threshold))

print("GRN Done")
print(f"Saved：result/GRN_filtered_brief_{sys_output_name}_raw.txt")
print(f"Saved：result/GRN_filtered_full_{sys_output_name}_raw.txt")
