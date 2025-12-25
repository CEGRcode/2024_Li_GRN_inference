import sys
import getopt
import os, re, json
import pandas as pd
import numpy as np
import random
import copy
import time
from scipy.integrate import solve_ivp
from dynamics import *
from distance_functions import *
from GRN_Expanded_Combinatorial_2025 import GRN
from mutation_functions import *
from population_functions import *
from utility_functions import *
from matplotlib import pyplot as plt
import itertools

def get_edges_from_json(json_file):
    # Dictionary to store edges
    edge_dict = {}
    # Load the JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)
    # Assuming edges are stored under the 'edges' key
    edges = data.get('edges', [])
    # Iterate over the edges and add them to the dictionary
    for edge in edges:
        source = edge['source']
        target = edge['target']
        if source in edge_dict:
            edge_dict[source].append(target)
        else:
            edge_dict[source] = [target]  # Store the edge as a tuple key
    return edge_dict

TF_DNA_net = get_edges_from_json('./data/Rossi_Ruihao_TF_DNA_union_all.json')

sys_path = os.getcwd()
sys_input_RNAseq = ""
sys_input_ChIP = ""
sys_input_Path = ""
sys_training_count = 3000
sys_PerturbationPower = 0
sys_output_name = ""
Order_of_genes = ''

try:
    opts, args = getopt.getopt(sys.argv[1:],
                                "hr:c::i:n::p::o:g:",
                                ["help", "input_RNAseq", "input_ChIP", "input_path",
                                "training_count", "PerturbationPower",
                                "output_name", "gene_names"])
except getopt.GetoptError:
    print('Usage: remove_dispensable_edges.py [-h] -r <input_RNAseq> [-c <input_ChIP>] -i <input_path> '
            ' [-n <training_count>]'
            ' [-p <PerturbationPower>]'
            ' -o <output_name> -g <gene_names>')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print('Usage: remove_dispensable_edges.py [-h] -r <input_RNAseq> [-c <input_ChIP>] -i <input_path> '
            ' [-n <training_count>]'
            ' [-p <PerturbationPower>]'
            ' -o <output_name> -g <gene_names>')
        sys.exit()
    elif opt in ("-r", "--input_RNAseq"):
        sys_input_RNAseq = sys_path + '/' + arg
        sys_input_RNAseq = np.array(pd.read_csv(
            sys_input_RNAseq, header=None, delimiter='\t', dtype=str))
        sys_WTTP = {}
        for i in range(0, sys_input_RNAseq.shape[0]):
            # print(sys_input_RNAseq[i][0], type(sys_input_RNAseq[i][0]))
            if sys_input_RNAseq[i][0] == '-1':
                sys_WTTP[str(i)] = [[], np.array(sys_input_RNAseq[i][1:], dtype=float)]
            elif len(sys_input_RNAseq[i][0]) == 1:
                sys_WTTP[str(i)] = [[int(sys_input_RNAseq[i][0])], np.array(sys_input_RNAseq[i][1:], dtype=float)]
            else:
                sys_WTTP[str(i)] = [np.array(sys_input_RNAseq[i][0].split(','), dtype=int).tolist(),
                                    np.array(sys_input_RNAseq[i][1:], dtype=float)]
    elif opt in ("-c", "--input_ChIP"):
        sys_input_ChIP = sys_path + '/' + arg
        sys_input_ChIP = np.array(pd.read_csv(sys_input_ChIP, header=None, delimiter='\t'), dtype=int)
    elif opt in ("-i", "--input_path"):
        sys_input_Path = arg
    elif opt in ("-n", "--training_count"):
        sys_training_count = arg
    elif opt in ("-p", "--PerturbationPower"):
        sys_PerturbationPower = arg
    elif opt in ("-o", "--output_name"):
        sys_output_name = arg
    elif opt in ("-g", "--gene_names"):
        with open(arg, "r") as file:
            first_row = file.readline().strip()
            Order_of_genes = [x for x in first_row.split("\t") if x]
    else:
        pass
########################################################################### Setting systematic Parameters #############################################################################
Result_Flow_Dict = {}
for filename in os.listdir(sys_input_Path):
    if filename.endswith('ResultFlow.txt'):
        filepath = os.path.join(sys_input_Path, filename)
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                lines = f.readlines()
                if lines:  # if file is not empty
                    last_row = lines[-1].strip()  # Remove trailing newline/whitespace
                    #print(f"File: {filename} -> Last row: {last_row}")
                    Result_Flow_Dict[filename.split('_')[1]] = list(map(float, last_row.split()))
    else:
        continue

GRN_Components_Dict = {}
for filename in os.listdir(sys_input_Path):
    if filename.endswith('_GRN_components.txt'):
        filepath = os.path.join(sys_input_Path, filename)
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

indexes_of_diff_genes = sorted(list(map(int, GRN_Components_Dict[next(iter(GRN_Components_Dict))].keys())))

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
for i, (key, value) in enumerate(zip(min_keys, min_values)):    
    #print(f"Target Gene {indexes_of_diff_genes[i]}: GRN = {key}, Min Distance = {value}")
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

TotalNumberOfGenes = len(AM_List[0])

if Order_of_genes == '':
    Order_of_genes = ['G{}'.format(i+1) for i in range(0, TotalNumberOfGenes)]
else:
    pass

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

TrainingCount = 3000
sys_PerturbationPower = 0.0

def ParametersInitiations(null_net):
    '''Parameters Initiations'''
    # Name
    # mRNA
    ####################################################
    RandomStringList = np.random.randint(2, size=(1, TotalNumberOfGenes ** 2))[0]
    RandomString = ''
    if null_net:
        for INT in RandomStringList:
            RandomString = RandomString + '0'
    else:
        for INT in RandomStringList:
            RandomString = RandomString + str(INT)

    #Configuration = String012ToMatrix(RandomString)
    Configuration = String012ToMatrix(Final_AM)
    ####################################################

    MutationRate = 0  # must be integer

    Sigmoid_k_init = Final_Sigmoid_k
    LogicGates = Final_LG

    TranscriptionThreshold = Final_Transcriptional_Threshold
    f0_c = Final_f0_c
    DegradationRatemRNA = Final_mRNA_DegradationRate
    Leakage = Final_Leakage
    TranscriptionRate = Final_Transcription_Rate

    return Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_k_init, Leakage, f0_c

#######################################################################################################################################################################################

########################################################################### Target transcriptional profiles #############################################################################
WTTP = sys_WTTP
TotalNumberOfGenes = len(Final_f0_c)
#######################################################################################################################################################################################

StableStatesCollector = {}

###########################
### PREPARING MAIN LOOP ###
###########################
BaseMatrixCollector = []
TranscriptionPofileMax = []
TranscriptionPofileMin = []
TranscriptionPofileAve = []
Overexpression = np.zeros((len(sys_WTTP), TotalNumberOfGenes))
Knockout = np.zeros((len(sys_WTTP), TotalNumberOfGenes))

for keys in sys_WTTP:
    for each in sys_WTTP[keys][0]:
        if each >= 0:
            Knockout[int(keys)][each] = sys_WTTP[keys][1][each]
            sys_WTTP[keys][1][each] = np.nan
        elif each <= -2:
            Overexpression[int(keys)][-each - 2] = sys_WTTP[keys][1][-each - 2]
            sys_WTTP[keys][1][-each - 2] = np.nan
        else:
            pass
    BaseMatrixCollector.append(sys_WTTP[keys][1])

BaseMatrixCollector = np.array(BaseMatrixCollector)

for row in BaseMatrixCollector.T:
    TranscriptionPofileMax.append(np.nanmax(row))
    TranscriptionPofileMin.append(np.nanmin(row))
    TranscriptionPofileAve.append(0.5 * (TranscriptionPofileMax[-1] + TranscriptionPofileMin[-1]))

### Get the genes having same TPM across all samples ###
indexes_of_diff_gene = []
for x in range(0, len(TranscriptionPofileMax)):
    if TranscriptionPofileMax[x] == TranscriptionPofileMin[x]:
        pass
    else:
        indexes_of_diff_gene.append(x)
Diff_gene_cannot_be_inferred = [4, 11]
indexes_of_diff_gene = [x for x in indexes_of_diff_gene if x not in Diff_gene_cannot_be_inferred]
Num_of_genes_diff_TPM = len(indexes_of_diff_gene)
### Done Getting the genes having same TPM across all samples ###

########################################################################### Setting systematic Parameters #############################################################################
(Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0_c) = ParametersInitiations(False)
this_GRN = GRN('this_GRN', [], Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0_c, sys_input_ChIP)
#######################################################################################################################################################################################
print('indexes_of_diff_gene: ', len(indexes_of_diff_gene), indexes_of_diff_gene)

filepath = os.path.join("result", 'temp_remove_dispensible_edge.txt')
with open(filepath, "r") as fh:
    first_line = fh.readline().rstrip("\n")

Final_AM = first_line.split()[1]
Final_LG = LogicGatesString2Matrix_Expanded(updateLGonAM(Final_AM, Final_LG))
#######################################
### SAVE ADDITIONAL OUTPUT TO FILES ###
#######################################
OutputFile = open(sys_path + '/result/' + 'GRN_filtered_brief_{}.txt'.format(sys_output_name), 'a')
OutputFile.write('{}\t{}\n'.format(Final_AM, Matrix2String01_LG_Expanded(Final_LG)))
OutputFile.close()

OutputFile = open(sys_path + '/result/' + 'GRN_filtered_full_{}.txt'.format(sys_output_name), 'a')
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
intersect_ = {}
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
#for each_source in TF_DNA_net:
#    for each_target in TF_DNA_net[each_source]:
#        if each_target in intersect_.get(each_source, []):
#            pass
#        else:
#            json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(each_source) + '\t  "target": "{}",\n'.format(each_target) + '\t  "label": "{}",\n'.format('') + '\t  "style": [\n\t\t"solid",\n\t\t"triangle",\n\t\t"#FFA500"\n\t  ]'+ '\n  \t},\n'

#print('From intersection: ', intersect_)
#print('From TF_DNA_net: ', TF_DNA_net)

outfile = open(sys_path+'/result/GRN_filtered_{}_2.json'.format(sys_output_name), 'a')
outfile.write(json_)
outfile.close()