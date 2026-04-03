import numpy as np
import os
from utility_functions import *
from pathlib import Path

def get_edges_from_json(json_file):
    edge_dict = {}
    if not os.path.exists(json_file):
        return edge_dict  # return empty dict
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

directory = './result/'
sys_output_name = 'Sc_GRN_final'

GRN_Components_Dict = {}
for filename in os.listdir(directory):
    if filename.endswith('_GRN_components.txt'):
        filepath = os.path.join(directory, filename)
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    line_split = (line.replace('\n', '')).split('\t')
                    if line_split[0] == 'GeneIndex':
                        pass
                    else:
                        if filename.split('_')[1] not in GRN_Components_Dict:
                            GRN_Components_Dict[filename.split('_')[1]] = {line_split[0]: line_split[1:]}
                        else:
                            GRN_Components_Dict[filename.split('_')[1]][line_split[0]] = line_split[1:]
    else:
        continue
print(GRN_Components_Dict.keys())

indexes_of_diff_genes = sorted(list(map(int, GRN_Components_Dict[next(iter(GRN_Components_Dict))].keys())))
outfile = open('./result/indexes_of_diff_gene.txt', 'w')
outfile.write('\t'.join(map(str, indexes_of_diff_genes)))
outfile.close()

file_path = Path('./data/GRN_ssTFs_column_names.txt')
Order_of_genes = []
if file_path.exists():
    with file_path.open() as f:
        for line in f:
            Order_of_genes = line.split()
if not Order_of_genes:
    Order_of_genes = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

print(Order_of_genes)
Order_of_genes = Order_of_genes[:len(indexes_of_diff_genes)]
print(Order_of_genes)

Result_Flow_Dict = {}
for filename in os.listdir(directory):
    if filename.endswith('ResultFlow.txt'):
        filepath = os.path.join(directory, filename)
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                lines = f.readlines()
                if lines:  # if file is not empty
                    last_row = lines[-1].strip()  # Remove trailing newline/whitespace
                    #print(f"File: {filename} -> Last row: {last_row}")
                    Result_Flow_Dict[filename.split('_')[1]] = list(map(float, last_row.split()))
    else:
        continue

print(Result_Flow_Dict)

TF_DNA_net = get_edges_from_json('./data/Rossi_Ruihao_TF_DNA_union_all.json')

# Transpose lists to iterate over positions
min_keys = []
min_values = []

for i, values in enumerate(zip(*Result_Flow_Dict.values())):
    min_key = min(Result_Flow_Dict.keys(), key=lambda k: Result_Flow_Dict[k][i])
    min_value = Result_Flow_Dict[min_key][i]
    min_keys.append(min_key)
    min_values.append(min_value)

AM_List = []
LG_List = []
f0_c_List = []
Sigmoid_k_List = []
Degradation_mRNA_List = []
Leakage_List = []
Transcription_Rate_List = []
Transcriptional_Threshold_List = []
TG_to_remove = []

for i, (key, value) in enumerate(zip(min_keys, min_values)):
    if value >= 0.5:
        TG_to_remove.append(indexes_of_diff_genes[i])
    else:
        continue
#print('TG_to_remove: ', TG_to_remove)
outfile = open('./result/Preliminary_training_result.txt', 'a')
for i, (key, value) in enumerate(zip(min_keys, min_values)):
    outfile.write(f"{indexes_of_diff_genes[i]}\t{value}\n")
    print(f"Target Gene {indexes_of_diff_genes[i]}: GRN = {key}, Min Distance = {value}")
    if i not in TG_to_remove:
        AM_List.append(GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][0])
        LG_List.append(GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][1])
        f0_c_List.append(list(map(float, GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][2].split(','))))
        Sigmoid_k_List.append(float(GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][3]))
        Leakage_List.append(float(GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][4]))
        Degradation_mRNA_List.append(float(GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][5]))
        Transcription_Rate_List.append(float(GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][6]))
        Transcriptional_Threshold_List.append(eval(GRN_Components_Dict[str(key)][str(indexes_of_diff_genes[i])][7]))
    else:
        continue
outfile.close()
#print(indexes_of_diff_genes)
indexes_of_diff_genes = [x for x in indexes_of_diff_genes if x not in TG_to_remove]
print(indexes_of_diff_genes)

TotalNumberOfGenes = len(AM_List[0])
Final_AM = ''
for j in range(0, TotalNumberOfGenes):
    for k in range(0, TotalNumberOfGenes):
        if k in indexes_of_diff_genes:
            Final_AM = Final_AM + str(AM_List[indexes_of_diff_genes.index(k)][j])
        else:
            Final_AM = Final_AM + '0'
Final_LG = []
Final_Sigmoid_k = []
Final_mRNA_DegradationRate = []
Final_Leakage = []
Final_f0_c = []
Final_Transcription_Rate = []
Final_Transcriptional_Threshold = []
for j in range(0, TotalNumberOfGenes):
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
        Final_LG.append([unit for unit in range(0, TotalNumberOfGenes)])
        Final_Sigmoid_k.append(0)
        Final_mRNA_DegradationRate.append(0)
        Final_Leakage.append(0)
        Final_f0_c.append([0,0,0,0,0,0])
        Final_Transcription_Rate.append(0)
        Final_Transcriptional_Threshold.append([[0,0,0] for q in range(0, TotalNumberOfGenes)])
#print(np.sqrt(len(Final_AM)))

#######################################
### SAVE ADDITIONAL OUTPUT TO FILES ###
#######################################
OutputFile = open(directory + 'GRN_filtered_brief_{}_raw.txt'.format(sys_output_name), 'a')
OutputFile.write('{}\t{}\n'.format(Final_AM, Matrix2String01_LG_Expanded(Final_LG)))
OutputFile.close()

OutputFile = open(directory + 'GRN_filtered_full_{}_raw.txt'.format(sys_output_name), 'a')
OutputFile.write(
    'Adjacency Matrix:\t{}\nLogic Gate:\t{}\nf0:\t{}\nHill Coefficient:\t{}\nTranscription Rate:\t{}\nmRNA Degradation Rate:\t{}\nLeakage Rate:\t{}\nTF Effective Threshold:\t{}\n\n'.format(
        Final_AM,
        Matrix2String01_LG_Expanded(Final_LG),
        Final_f0_c,
        Final_Sigmoid_k,
        Final_Transcription_Rate,
        Final_mRNA_DegradationRate,
        Final_Leakage,
        Final_Transcriptional_Threshold))
OutputFile.close()

###########################
### SAVE IN JSON FORMAT ###
###########################
AM = Final_AM
LG = Final_LG
f0_c = eval(str(Final_f0_c))
DegradationRatemRNA = eval(str(Final_mRNA_DegradationRate))
Leakage = eval(str(Final_Leakage))
TranscriptionRate = eval(str(Final_Transcription_Rate))
Sigmoid_k_init = eval(str(Final_Sigmoid_k))
TranscriptionThreshold = eval(str(Final_Transcriptional_Threshold))
for i in range(0, len(TranscriptionThreshold)):
    for j in range(0, len(TranscriptionThreshold[i])):
        if TranscriptionThreshold[i][j][0] == TranscriptionThreshold[i][j][1] and TranscriptionThreshold[i][j][1]  == TranscriptionThreshold[i][j][2]:
            TranscriptionThreshold[i][j] = [0.5, 0.5, 0.5]
        elif TranscriptionThreshold[i][j][0] == TranscriptionThreshold[i][j][1] and TranscriptionThreshold[i][j][0] != TranscriptionThreshold[i][j][2]:
            TranscriptionThreshold[i][j] = [0.66, 0.66, 0.33]
        elif TranscriptionThreshold[i][j][0] != TranscriptionThreshold[i][j][1] and TranscriptionThreshold[i][j][0] == TranscriptionThreshold[i][j][2]:
            TranscriptionThreshold[i][j] = [0.66, 0.33, 0.66]
        else:
            raise Exception('Error.')

if type(LG) == str:
    LG = LG.split(',')
else:
    LG = [item for sublist in LG for item in sublist]

json_ = ''
json_ = json_ + '{\n  "nodes": [\n'
for i in range(0, len(Order_of_genes)):
    if i in indexes_of_diff_genes:
        json_ = json_ + '\t{\n\t  ' + '"id": "{}",\n'.format(Order_of_genes[i]) + '\t  "label": "{}",\n'.format(Order_of_genes[i]) + '\t  "TR": {},\n'.format(TranscriptionRate[i]) + '\t  "f0": {},\n'.format(f0_c[i][0]) + '\t  "f0p": {},\n'.format(f0_c[i][1])  + '\t  "t1": {},\n'.format(TranscriptionThreshold[i][0][0]) + '\t  "t2": {},\n'.format(TranscriptionThreshold[i][0][1]) + '\t  "t3": {},\n'.format(TranscriptionThreshold[i][0][2]) + '\t  "c1": {},\n'.format(f0_c[i][2]) + '\t  "c2": {},\n'.format(f0_c[i][3]) + '\t  "c3": {},\n'.format(f0_c[i][4]) + '\t  "c4": {},\n'.format(f0_c[i][5]) + '\t  "k": {},\n'.format(Sigmoid_k_init[i]) + '\t  "Lk": {},\n'.format(Leakage[i]) + '\t  "Deg": {},\n'.format(DegradationRatemRNA[i]) +  '\t  "color": \"black\",\n' +  '\t  "background_color": \"white\"' + '\n  \t},\n'
    else:
        json_ = json_ + '\t{\n\t  ' + '"id": "{}",\n'.format(Order_of_genes[i]) + '\t  "label": "{}",\n'.format(Order_of_genes[i]) + '\t  "TR": {},\n'.format(TranscriptionRate[i]) + '\t  "f0": {},\n'.format(f0_c[i][0]) + '\t  "f0p": {},\n'.format(f0_c[i][1])  + '\t  "t1": {},\n'.format(TranscriptionThreshold[i][0][0]) + '\t  "t2": {},\n'.format(TranscriptionThreshold[i][0][1]) + '\t  "t3": {},\n'.format(TranscriptionThreshold[i][0][2]) + '\t  "c1": {},\n'.format(f0_c[i][2]) + '\t  "c2": {},\n'.format(f0_c[i][3]) + '\t  "c3": {},\n'.format(f0_c[i][4]) + '\t  "c4": {},\n'.format(f0_c[i][5]) + '\t  "k": {},\n'.format(Sigmoid_k_init[i]) + '\t  "Lk": {},\n'.format(Leakage[i]) + '\t  "Deg": {}'.format(DegradationRatemRNA[i]) + '\n  \t},\n'
json_ = json_[:-2]

# Add the edges
json_ = json_ + '\n  ],\n  "edges": [\n'
num_nodes = len(f0_c)
#print('AM: ', AM, len(AM))
for i in range(0, len(AM)):
    if AM[i] == '0':
        continue
    elif AM[i] == '1':
        LG_slide = LG[num_nodes*(i%num_nodes):num_nodes*(1+i%num_nodes)]
        #print(i, LG_slide, i//num_nodes)
        subunit_indexes = [x for x, v in enumerate(LG_slide) if v == LG_slide[i//num_nodes]]
        activator_indexes = [x for x, v in enumerate(AM[(i%num_nodes)::num_nodes]) if v == '1']
        existing_subunits = list(set(subunit_indexes) & set(indexes_of_diff_genes) & set(activator_indexes))
        #print('->', activator_indexes, '\n')
        if LG_slide.count(LG_slide[i//num_nodes]) > 1 and len(existing_subunits) > 1:
            if Order_of_genes[i%num_nodes] in TF_DNA_net.get(Order_of_genes[i//num_nodes], []):
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
            else:
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle"\n\t,\n\t\t"#000000"\n\t  ]'+ '\n  \t},\n'
        else:
            if Order_of_genes[i%num_nodes] in TF_DNA_net.get(Order_of_genes[i//num_nodes], []):
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"triangle"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
            else:
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"triangle"\n\t,\n\t\t"#000000"\n\t  ]'+ '\n  \t},\n'
    else:
        LG_slide = LG[num_nodes*(i%num_nodes):num_nodes*(1+i%num_nodes)]
        subunit_indexes = [x for x, v in enumerate(LG_slide) if v == LG_slide[i//num_nodes]]
        inhibitor_indexes = [x for x, v in enumerate(AM[(i%num_nodes)::num_nodes]) if v == '2']
        existing_subunits = list(set(subunit_indexes) & set(indexes_of_diff_genes) & set(inhibitor_indexes))
        #print('-->', existing_subunits, subunit_indexes, '\n')
        if LG_slide.count(LG_slide[i//num_nodes]) > 1 and len(existing_subunits) > 1:
            if Order_of_genes[i%num_nodes] in TF_DNA_net.get(Order_of_genes[i//num_nodes], []):
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"tee"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
            else:
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"tee"\n\t,\n\t\t"#000000"\n\t  ]'+ '\n  \t},\n'
        else:
            if Order_of_genes[i%num_nodes] in TF_DNA_net.get(Order_of_genes[i//num_nodes], []):
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"tee"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
            else:
                json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"tee"\n\t,\n\t\t"#000000"\n\t  ]'+ '\n  \t},\n'
json_ = json_[:-2]
json_ = json_ + '\n  ]\n}'

outfile = open(directory+'/GRN_filtered_{}_raw.json'.format(sys_output_name), 'a')
outfile.write(json_)
outfile.close()
