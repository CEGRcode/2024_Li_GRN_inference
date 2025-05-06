import numpy as np
import pandas as pd
import math
import os
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

# get the ssTFs to build the TF-DNA binding network
inputfile = open('./data/ssTFs_common_names.txt', 'r')
ssTF_names = []
for line in inputfile:
    ssTF_names.append(line.split()[0])
inputfile.close()

OutMatrix = []
GeneNames.index(ssTF_names[0])
for each_row in ssTF_names:
    OutMatrix.append([])
    for each_column in ssTF_names:
        OutMatrix[-1].append(TF_Gene_Matrix[each_row][GeneNames.index(each_column)])
OutMatrix = pd.DataFrame(np.array(OutMatrix).T, index=ssTF_names, columns=ssTF_names)

# output to JSON
# nodes:
json_output = ''
json_output = json_output + '{\n  "nodes": [\n'
for i in range(0, len(ssTF_names)):
    json_output = json_output + '\t{\n\t  ' + '"id": "{}",\n'.format(ssTF_names[i]) + '\t  "label": "{}",\n'.format(ssTF_names[i]) + '\t  "sua7Occupancy": {}\n'.format(1) + '\n  \t},\n'
json_output = json_output[:-2]


# edges:
json_output = json_output + '\n  ],\n  "edges": [\n'
for each in np.argwhere(OutMatrix  == 1):
    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format(ssTF_names[each[0]]) + '\t  "target": "{}",\n'.format(ssTF_names[each[1]]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"triangle"\n\t  ]'+ '\n  \t},\n'
json_output = json_output[:-2]
json_output = json_output + '\n  ]\n}'
outfile = open('./result/TF_DNA_binding_network_{}.json'.format(datetime.today().strftime('%d%m%y')), 'a')
outfile.write(json_output)
outfile.close()