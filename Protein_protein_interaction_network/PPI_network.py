import numpy as np
import math
import pandas as pd
import copy
import os
import re
import sys

def Get_stats_in_MEME(input_file):
    #print(input_file)
    output_ = []
    with open(input_file, "r") as file:
        for line in file:
            # Match lines with the desired format
            match = re.search(r"MOTIF\s+(\d+).*?sites\s+=\s+(\d+).*?E-value\s+=\s+([\d\.eE+-]+)", line)
            if match:
                motif = match.group(1)
                sites = match.group(2)
                e_value = match.group(3)
                output_.append([motif, sites, e_value])
    return output_


############################################# Parameter initialization #############################################
Overall_CutOff = False # 'default'. Otherwise set a number between 0 and 1 as the JS factor-specific cutoff.
Const_CutOff = 0 # JS absolute cutoff
dic_PPI = {}
Number_Of_Edges = 0
Motif_based_PPI = sys.argv[1].strip().lower() == 'true' # True if constructing Motif-based PPI; False if conscructing PPI for all binding sites.
########################################################################################################################
########################################################################################################################


############################### Blacklist the low quality motifs for the Motif-based PPI ###############################
if Motif_based_PPI:
    infile = open('./QC_count.txt', 'r')
    Motif_Quality_Dic = {}

    for line in infile:
        #print("line: ", line)
        Motif_num = re.search(r"Motif_(\d+)_peaks", line)
        if Motif_num:
            #print('Factor_ID: ', line.split(':')[0], 'Motif number: ', Motif_num.group(1))
            if os.path.isfile('./Motifs_all/{}_MEME_Motifs.txt'.format(line.split(':')[0])):
                #print(Get_stats_in_MEME('/Users/rl884/Downloads/241119_motif_based_JS/Motifs_all/{}_MEME_Motifs.txt'.format(line.split(':')[0]))[-1+int(Motif_num.group(1))])
                Motif_Quality_Dic[line.split(':')[0]+'_'+Motif_num.group(1)] = [Get_stats_in_MEME('./Motifs_all/{}_MEME_Motifs.txt'.format(line.split(':')[0]))[-1+int(Motif_num.group(1))][0], Get_stats_in_MEME('./Motifs_all/{}_MEME_Motifs.txt'.format(line.split(':')[0]))[-1+int(Motif_num.group(1))][-1]]
            else:
                Motif_Quality_Dic[line.split(':')[0]+'_'+Motif_num.group(1)] = []
        else:
            pass
        #print('****************************************************************\n\n')
    infile.close()

    Factor_blacklist = []
    for keys in Motif_Quality_Dic:
        #print(keys, Motif_Quality_Dic[keys])
        if Motif_Quality_Dic[keys] == []:
            pass
        elif float(Motif_Quality_Dic[keys][-1]) >= 0.1:
            Factor_blacklist.append(keys.split('_')[0])
        else:
            pass
else:
    Factor_blacklist = []
print(len(Factor_blacklist))
########################################################################################################################
########################################################################################################################


############################################# Construct PPI from ChIP-exo data #########################################
BED_itself = []
The_Best_Replicate = []
The_Best_Replicate_Name = {}
infile = open('./YEP_best_rep.txt', 'r')
for line in infile:
    line_split = line.split()
    The_Best_Replicate.append(line_split[1])
    The_Best_Replicate_Name[line_split[1]] = line_split[0].capitalize()
infile.close()

score_recorder = []
for filename in os.listdir('./Second_sort/'):
    Content_Dic = {}
    if filename.endswith('.txt'):
        infile = open('./Second_sort/'+filename, 'r')
        BED_ID = filename.split('_')[0]
        if BED_ID not in The_Best_Replicate or BED_ID in Factor_blacklist:
            continue
        else:
            for line in infile:
                line_split = line.replace("\n", "").split('\t')
                if line_split[0] in The_Best_Replicate:
                    Content_Dic[line_split[0]] = []
                    for i in range(1, len(line_split)):
                        Content_Dic[line_split[0]].append(line_split[i])
                else:
                    pass
            infile.close()
            if BED_ID not in Content_Dic:
                print('skipped for itself not present', BED_ID)
                continue
            elif Content_Dic[BED_ID][0] == 'Bad Quality' or float(Content_Dic[BED_ID][-1]) <= 3:
                print('skipped for bad quality', BED_ID)
                continue
            else:
                BED_itself.append(float(Content_Dic[BED_ID][-1]))
                pass

        for each_key in list(Content_Dic.keys()):
            if Content_Dic[each_key][-1] == 'Bad Quality':
                del Content_Dic[each_key]
        
        sorted_keys = sorted(Content_Dic, key=lambda x: float(Content_Dic[x][-1]), reverse=True)

        max_score = Content_Dic[sorted_keys[0]][-1]
        if max_score == 'Bad Quality':
            print('skipped for bad quality', BED_ID)
            continue
        else:
            pass
        
        if Overall_CutOff:
            print(np.round((1+sorted_keys.index(BED_ID))/len(sorted_keys), 3), len(sorted_keys), BED_ID)
            BED_ID_Cutoff = Const_CutOff*(float(max_score)*((1+sorted_keys.index(BED_ID))/len(sorted_keys)))
        else:
            BED_ID_Cutoff = Const_CutOff*float(max_score)
        #print(max_score, BED_ID_Cutoff)
        for keys in Content_Dic:
            if keys == BED_ID or Content_Dic[keys][-1] == 'Bad Quality':
                pass
            else:
                #print(keys, The_Best_Replicate_Name[keys])
                if float(Content_Dic[keys][-1]) >= max([BED_ID_Cutoff, 3]):
                    #print(BED_ID, keys, float(Content_Dic[keys][-1]), BED_ID_Cutoff)
                    score_recorder.append(float(Content_Dic[keys][-1]))
                    Number_Of_Edges = Number_Of_Edges + 1
                    if The_Best_Replicate_Name[BED_ID] not in dic_PPI:
                        dic_PPI[The_Best_Replicate_Name[BED_ID]] = [The_Best_Replicate_Name[keys]]
                    else:
                        dic_PPI[The_Best_Replicate_Name[BED_ID]].append(The_Best_Replicate_Name[keys])
                else: 
                    pass
    else:
        continue

print(Number_Of_Edges)
print(len(dic_PPI))
########################################################################################################################
########################################################################################################################


############################################# Construct PPI from ChIP-exo data #########################################
STRING_ID_file = open('./4932.protein.info.v12.0.txt', 'r')
STRING_ID = {}
for line in STRING_ID_file:
    line_split = line.split('\t')
    STRING_ID[line_split[0]] = line_split[1]
    #print(line_split)
STRING_ID_file.close()

YEP_Factors = []
for each in The_Best_Replicate_Name:
    YEP_Factors.append(The_Best_Replicate_Name[each])

STRING_DATA_file = open('./4932.protein.links.v12.0.tsv', 'r')
temp_Dic_very_high = {}
temp_Dic_high = {}
temp_Dic_medium = {}
temp_Dic_exploratory = {}
for line in STRING_DATA_file:
    line_split = line.split()
    if line_split[0] not in STRING_ID or line_split[1] not in STRING_ID or STRING_ID[line_split[0]].capitalize() not in YEP_Factors or STRING_ID[line_split[1]].capitalize() not in YEP_Factors:
        continue
    else:
        pass
    if int(line_split[2]) >= 900:
        if STRING_ID[line_split[0]].capitalize() not in temp_Dic_very_high:
            temp_Dic_very_high[STRING_ID[line_split[0]].capitalize()] = [STRING_ID[line_split[1]].capitalize()]
        else:
            temp_Dic_very_high[STRING_ID[line_split[0]].capitalize()].append(STRING_ID[line_split[1]].capitalize())
    elif int(line_split[2]) >= 700 and int(line_split[2]) < 900:
        if STRING_ID[line_split[0]].capitalize() not in temp_Dic_high:
            temp_Dic_high[STRING_ID[line_split[0]].capitalize()] = [STRING_ID[line_split[1]].capitalize()]
        else:
            temp_Dic_high[STRING_ID[line_split[0]].capitalize()].append(STRING_ID[line_split[1]].capitalize())
    elif int(line_split[2]) >= 400 and int(line_split[2]) < 700:
        if STRING_ID[line_split[0]].capitalize() not in temp_Dic_medium:
            temp_Dic_medium[STRING_ID[line_split[0]].capitalize()] = [STRING_ID[line_split[1]].capitalize()]
        else:
            temp_Dic_medium[STRING_ID[line_split[0]].capitalize()].append(STRING_ID[line_split[1]].capitalize())
    elif int(line_split[2]) >= 0 and int(line_split[2]) < 400:
        if STRING_ID[line_split[0]].capitalize() not in temp_Dic_exploratory:
            temp_Dic_exploratory[STRING_ID[line_split[0]].capitalize()] = [STRING_ID[line_split[1]].capitalize()]
        else:
            temp_Dic_exploratory[STRING_ID[line_split[0]].capitalize()].append(STRING_ID[line_split[1]].capitalize())
    else:
        pass
STRING_DATA_file.close()

num_very_high = 0
num_high = 0
num_medium = 0
num_exploratory = 0
for keys in temp_Dic_very_high:
    temp_Dic_very_high[keys] = set(temp_Dic_very_high[keys])
    num_very_high = num_very_high + len(temp_Dic_very_high[keys])
    
for keys in temp_Dic_high:
    temp_Dic_high[keys] = set(temp_Dic_high[keys])
    num_high = num_high + len(temp_Dic_high[keys])

for keys in temp_Dic_medium:
    temp_Dic_medium[keys] = set(temp_Dic_medium[keys])
    num_medium = num_medium + len(temp_Dic_medium[keys])
    
for keys in temp_Dic_exploratory:
    temp_Dic_exploratory[keys] = set(temp_Dic_exploratory[keys])
    num_exploratory = num_exploratory + len(temp_Dic_exploratory[keys])

print('very high: ', num_very_high, '; high: ', num_high, '; medium: ', num_medium, '; exploratory: ', num_exploratory)

# Analyze bed files for the overlapping binding sites.
All_factors = ['AFT1']
for keys in dic_PPI:
    if keys in All_factors or keys == 'Rcs1':
        continue
    else:
        All_factors.append(keys)
    for each_factor in dic_PPI[keys]:
        if each_factor not in All_factors:
            All_factors.append(each_factor)
        else:
            pass
approximation_threshold = 100

BED_Dic = {}
if Motif_based_PPI:
    dir_string = './YEP_Motif_BED'
else:
    dir_string = './YEP_ALL_BED'
for filename in os.listdir(dir_string):
    if filename.split('_')[0] not in The_Best_Replicate:
        continue
    else:
        pass 
    with open(os.path.join(dir_string, filename), 'r', encoding='utf-8', errors='ignore') as infile:
        BED_Dic[filename.split('_')[0]] = {}
        for line in infile:
            if line.split()[0] not in BED_Dic[filename.split('_')[0]]:
                BED_Dic[filename.split('_')[0]][line.split()[0]] = [line.split()[1]]
            else:
                BED_Dic[filename.split('_')[0]][line.split()[0]].append(line.split()[1])
                
total_number_of_peaks = {}
positive_hits = {}
for factor1 in BED_Dic:
    total_number_of_peaks[factor1] = sum([len(BED_Dic[factor1][x]) for x in BED_Dic[factor1]])
    for factor2 in BED_Dic:
        if factor1 == factor2:
            continue
        else:
            for chrom in BED_Dic[factor1]:
                if chrom not in BED_Dic[factor2]:
                    pass
                else:
                    for binding_site1 in BED_Dic[factor1][chrom]:
                        distance_list = []
                        for binding_site2 in BED_Dic[factor2][chrom]:
                            distance_list.append(abs(int(binding_site1)-int(binding_site2)))
                        if min(distance_list) <= approximation_threshold:
                            if factor1+'_'+factor2 in positive_hits:
                                positive_hits[factor1+'_'+factor2] = positive_hits[factor1+'_'+factor2] + 1
                            else:
                                positive_hits[factor1+'_'+factor2] = 1
                        else:
                            pass

Proportion_binding_sites = 0 # threshold for the proportion of overlapping binding sites

dic_PPI_archive = copy.deepcopy(dic_PPI)
dic_PPI = {}
for keys in dic_PPI_archive:
    if keys == 'Rcs1':
        pass
    else:
        key_to_use = keys
    if key_to_use not in All_factors:
        continue
    else:
        pass
    
    for each in dic_PPI_archive[keys]:
        if each == keys or each not in All_factors:
            pass
        else:
            if [key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0] not in total_number_of_peaks:
                pass
            else:
                if total_number_of_peaks[[key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0]] < 0:
                    pass
                if [key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0]+'_'+[key for key, value in The_Best_Replicate_Name.items() if value == each][0] not in positive_hits or [key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0] not in total_number_of_peaks:
                    pass
                elif positive_hits[[key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0]+'_'+[key for key, value in The_Best_Replicate_Name.items() if value == each][0]]/total_number_of_peaks[[key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0]] < Proportion_binding_sites:
                    pass
                else:
                    if keys in dic_PPI:
                        dic_PPI[keys].append(each)
                    else:
                        dic_PPI[keys] = [each]
print(len(dic_PPI_archive), '=>', len(dic_PPI))
########################################################################################################################
########################################################################################################################


####################################### Convert dic_PPI to .json just for all factors with colors indicating comparison outcome with STRING and opacity indicateing number of overlapped ChexMix called Peaks #######################################
ssTF_names = []
infile = open('./ssTFs_common_names.txt', 'r')
for line in infile:
    ssTF_names.append(line.split()[0])
infile.close()

dic_PPI_ssTFs = {}
for each_factor in dic_PPI:
    if each_factor.upper() not in ssTF_names:
        continue
    else:
        dic_PPI_ssTFs[each_factor] = []
        for each_target in dic_PPI[each_factor]:
            if each_target.upper() in ssTF_names:
                dic_PPI_ssTFs[each_factor].append(each_target)
            else:
                pass
dic_PPI = copy.deepcopy(dic_PPI_ssTFs)

# nodes:
json_output = ''
json_output = json_output + '{\n  "nodes": [\n'
for i in range(0, len(ssTF_names)):
    json_output = json_output + '\t{\n\t  ' + '"id": "{}",\n'.format(ssTF_names[i].upper()) + '\t  "label": "{}",\n'.format(ssTF_names[i].upper()) + '\t  "sua7Occupancy": {},\n'.format(1) + '\t  "f0": {}'.format(0)+ '\n  \t},\n'
json_output = json_output[:-2]

# edges:
json_output = json_output + '\n  ],\n  "edges": [\n'
for keys in dic_PPI:
    if keys == 'Rcs1':
        key_to_use = 'AFT1'
    else:
        key_to_use = keys
    if key_to_use not in All_factors:
        continue
    else:
        pass
    
    for each in dic_PPI[keys]:
        if each == keys or each not in All_factors:
            pass
        else:
            if [key for key, value in The_Best_Replicate_Name.items() if value == each] == [] or [key for key, value in The_Best_Replicate_Name.items() if value == key_to_use] == []:
                Opacity = 0    
            elif [key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0]+'_'+[key for key, value in The_Best_Replicate_Name.items() if value == each][0] in positive_hits:
                Opacity = positive_hits[[key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0]+'_'+[key for key, value in The_Best_Replicate_Name.items() if value == each][0]]/total_number_of_peaks[[key for key, value in The_Best_Replicate_Name.items() if value == key_to_use][0]]
            else:
                Opacity = 0
            if each == 'RCS1':
                if 'Aft1' in temp_Dic_very_high and key_to_use.capitalize() in temp_Dic_very_high['Aft1']:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format('AFT1') + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#00A86B"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                elif 'Aft1' in temp_Dic_high and key_to_use.capitalize() in temp_Dic_high['Aft1']:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format('AFT1') + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#F7F700"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                elif 'Aft1' in temp_Dic_medium and key_to_use.capitalize() in temp_Dic_medium['Aft1']:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format('AFT1') + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#FFA500"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                elif 'Aft1' in temp_Dic_exploratory and key_to_use.capitalize() in temp_Dic_exploratory['Aft1']:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format('AFT1') + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#DC143C"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                else:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format('AFT1') + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle"\n\t,\n\t\t"#000000"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
            else:
                if key_to_use in temp_Dic_very_high and each.capitalize() in temp_Dic_very_high[key_to_use]:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format(each.upper()) + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#00A86B"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                elif key_to_use in temp_Dic_high and each.capitalize() in temp_Dic_high[key_to_use]:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format(each.upper()) + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#F7F700"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                elif key_to_use in temp_Dic_medium and each.capitalize() in temp_Dic_medium[key_to_use]:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format(each.upper()) + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#FFA500"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                elif key_to_use in temp_Dic_exploratory and each.capitalize() in temp_Dic_exploratory[key_to_use]:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format(each.upper()) + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle",\n\t\t"#DC143C"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
                else:
                    json_output = json_output + '\t{\n\t  ' + '"source": "{}",\n'.format(each.upper()) + '\t  "target": "{}",\n'.format(key_to_use.upper()) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle"\n\t,\n\t\t"#000000"\n\t,\n\t\t"{}"\n\t  ]'.format(Opacity)+ '\n  \t},\n'
json_output = json_output[:-2]
json_output = json_output + '\n  ]\n}'

if Motif_based_PPI:
    output_name_str = 'motif_sites'
else:
    output_name_str = 'all_binding_sites'

outfile = open('./SCALE_cdt/Results/PPI_network_Cutoff_{}_STRING_overlapping_{}_{}_ssTFs_2025.json'.format(Const_CutOff, output_name_str, Proportion_binding_sites), 'a')
outfile.write(json_output)
outfile.close()
#####################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################
