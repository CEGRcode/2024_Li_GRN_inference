import numpy as np
import pandas as pd
import math
import os

def Define_BiDirectional_Region(genename_, chr_, strand_, left_, right_, window_size):
    genome_file = open('./result/Sc_genome_annotations.txt','r')
    output = {}
    for line in genome_file:
        if line.split('\t')[0] != chr_ or line.split('\t')[2] == '' or line.split('\t')[2] == 'NA':
            continue
        else:
            if strand_ == '+' and line.split('\t')[1] == '-':
                if 0< -int(float(line.split('\t')[3]))+left_ <= 500:
                    if line.split('\t')[5] == '':
                        if '{}_{}'.format(genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_), line.split()[4]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3]) in output:
                            pass
                        else:
                            output['{}_{}'.format(genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_), line.split()[4]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3])] = abs(int(float(line.split('\t')[3]))-left_)
                    else:
                        if '{}_{}'.format(genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_), line.split()[5]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3]) in output:
                            pass
                        else:
                            output['{}_{}'.format(genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_), line.split()[5]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3])] = abs(int(float(line.split('\t')[3]))-left_)
            elif strand_ == '-' and line.split('\t')[1] == '+':
                if 0 < int(float(line.split('\t')[2]))-right_ <= 500:
                    if line.split('\t')[5] == '':
                        if '{}_{}'.format(line.split()[4]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3], genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_)) in output:
                            pass
                        else:
                            output['{}_{}'.format(line.split()[4]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3], genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_))] = abs(int(float(line.split('\t')[2]))-right_)
                    else:
                        if '{}_{}'.format(line.split()[5]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3], genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_)) in output:
                            pass
                        else:
                            output['{}_{}'.format(line.split()[5]+'_'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2]+'_'+line.split('\t')[3], genename_+'_'+chr_+'_'+strand_+'_'+str(left_)+'_'+str(right_))] = abs(int(float(line.split('\t')[2]))-right_)
            else:
                continue
    genome_file.close()
    return output

def Find_intermediate_nucleosome(chr_, left_, right_):
    output_ = []
    
    nucleosome_file = open('./result/PlusOneNuc.txt','r')
    for line in nucleosome_file:
        if line[:3] == 'nuc':
            if line.split()[0].split('_')[1].split('-')[0] == chr_:
                if left_ < int(line.split()[0].split('_')[1].split('-')[1]) < right_:
                    output_.append(line.split()[0].split('_')[1])
            else:
                continue
        else:
            continue
    nucleosome_file.close()
    
    nucleosome_file = open('./result/MinusOneNuc.txt','r')
    for line in nucleosome_file:
        if line[:3] == 'nuc':
            if line.split()[0].split('_')[1].split('-')[0] == chr_:
                if left_ < int(line.split()[0].split('_')[1].split('-')[1]) < right_:
                    output_.append(line.split()[0].split('_')[1])
            else:
                continue
        else:
            continue
    nucleosome_file.close()
    return set(output_)

def Find_intermediate_TFs(chr_, left_, right_, YEP_ID):
    output_ = []
    bed_file = open('./data/{}_YEP/{}_Motif_1_bound.bed'.format(YEP_ID, YEP_ID),'r')
    for line in bed_file:
        if chr_ == line.split()[0] and left_ < 0.5*(int(line.split()[1])+int(line.split()[2])) < right_:
            output_.append(line.split()[0]+'_'+line.split()[1]+'_'+line.split()[2])
        else:
            pass
    bed_file.close()
    return output_

# Get the YEP ID for each factor
YEP_ID_dic = {}
inputfile = open('./data/YEP_best_rep.txt', 'r')
for line in inputfile:
    YEP_ID_dic[line.split()[0].upper()] = line.split()[1]
    if line.split()[0].upper() == 'RCS1':
        YEP_ID_dic['AFT1'] = line.split()[1]
    else:
        pass
inputfile.close()

# Define Head-to-head divergent genes
Bidirectional_Regions = {}
inputfile = open('./result/Sc_genome_annotations.txt','r')
for line in inputfile:
    if line.split('\t')[0] == 'Chrom' or line.split('\t')[2] == '' or line.split('\t')[2] == 'NA':
        continue
    elif line.split('\t')[5] == '':
        Bidirectional_Regions.update(Define_BiDirectional_Region(line.split('\t')[4], line.split('\t')[0], line.split('\t')[1], int(float(line.split('\t')[2])), int(float(line.split('\t')[3])), 500))
    else:
        Bidirectional_Regions.update(Define_BiDirectional_Region(line.split('\t')[5], line.split('\t')[0], line.split('\t')[1], int(float(line.split('\t')[2])), int(float(line.split('\t')[3])), 500))

outfile = open('./result/Divergent_Genes.txt','a')
for keys in Bidirectional_Regions:
    for each in keys.split('_'):
        outfile.write(each+'\t')
    outfile.write('\n')
outfile.close()

# Find the nucleosomes inbetween
Divergent_gene_file = open('./result/Divergent_Genes.txt','r')
Divergent_gene_nuc_file = open('./result/Divergent_Genes_nuc.txt','a')
for line in Divergent_gene_file:
    #print(Find_intermediate_nucleosome(line.split()[1], int(line.split()[9]), int(line.split()[3])))
    Divergent_gene_nuc_file.write(line[:-1])
    for each in Find_intermediate_nucleosome(line.split()[1], int(float(line.split()[9])), int(float(line.split()[3]))):
        Divergent_gene_nuc_file.write('nuc_'+each.split('-')[0]+'_'+each.split('-')[1]+'\t')
    Divergent_gene_nuc_file.write('\n')
Divergent_gene_file.close()
Divergent_gene_nuc_file.close()

# Find the TFs inbetween
inputfile = open('./data/ssTFs_common_names.txt', 'r')
ssTF_names = []
for line in inputfile:
    ssTF_names.append(line.split()[0])
inputfile.close()

for each in ssTF_names:
    if not os.path.isfile('./data/{}_YEP/{}_Motif_1_bound.bed'.format(YEP_ID_dic[each.title().upper()], YEP_ID_dic[each.title().upper()])):
        print('No CX.bed for '+each.title())
        continue
    else:
        pass
    Divergent_gene_file = open('./result/Divergent_Genes_nuc.txt','r')
    outputfile = open('./result/Divergent_Genes_nuc_1.txt','a')
    for line in Divergent_gene_file:
        #print(Find_intermediate_nucleosome(line.split()[1], int(line.split()[9]), int(line.split()[3])))
        outputfile.write(line[:-1])
        for each_ in Find_intermediate_TFs(line.split()[1], int(float(line.split()[9])), int(float(line.split()[3])), YEP_ID_dic[each.title().upper()]):
            outputfile.write('{}_'.format(each.title())+each_+'\t')
        outputfile.write('\n')
    Divergent_gene_file.close()
    outputfile.close()
    os.unlink('./result/Divergent_Genes_nuc.txt')
    os.rename('./result/Divergent_Genes_nuc_1.txt','./result/Divergent_Genes_nuc.txt')

# Calculate the % of the divergent genes that have at least one of the three insulators and another ssTF
os.rename('./result/Divergent_Genes_nuc.txt','./result/Divergent_Genes_and_Stuff_inbetween.txt')
inputfile = open('./result/Divergent_Genes_and_Stuff_inbetween.txt', 'r')
counter = 0
for line in inputfile:
    Is_other_TF = False
    if 'Reb1' in line or 'Rap1' in line or 'Abf1' in line:
        for each in line.split()[10:]:
            if each[:4] not in ['nuc_', 'Reb1', 'Rap1', 'Abf1']:
                Is_other_TF = True
            else:
                pass
    else:
        pass
    if Is_other_TF:
        counter = counter + 1
    else:
        pass
inputfile.close()

# Find Tandem Genes and NDR or NFRs
inputfile = open('./result/Divergent_Genes.txt', 'r')
Divergent_genes = []
for line in inputfile:
    Divergent_genes.append(line.split()[0])
    Divergent_genes.append(line.split()[5])
inputfile.close()

inputfile = open('./result/Sc_genome_annotations.txt', 'r')
outputfile = open('./result/Tandem_Genes.txt', 'a')
for line in inputfile:
    if line.split('\t')[5] == '':
        genename_ = line.split('\t')[4]
    else:
        genename_ = line.split('\t')[5]
    if genename_ not in Divergent_genes:
        outputfile.write(line)
    else:
        continue
inputfile.close()
outputfile.close()

Tandem_gene_file = open('./result/Tandem_Genes.txt','r')
Tandem_gene_nuc_file = open('./result/Tandem_Genes_nuc.txt','a')
for line in Tandem_gene_file:
    if line.split()[0] == 'Chrom':
        continue
    else:
        #print(Find_intermediate_nucleosome(line.split()[1], int(line.split()[9]), int(line.split()[3])))
        Tandem_gene_nuc_file.write(line[:-1]+'\t')
        for each in Find_intermediate_nucleosome(line.split()[0], int(float(line.split()[-2])), int(float(line.split()[-1]))):
            Tandem_gene_nuc_file.write('nuc_'+each.split('-')[0]+'-'+each.split('-')[1]+'\t')
        Tandem_gene_nuc_file.write('\n')
Tandem_gene_file.close()
Tandem_gene_nuc_file.close()

inputfile = open('./data/ssTFs_common_names.txt', 'r')
ssTF_names = []
for line in inputfile:
    ssTF_names.append(line.split()[0])
inputfile.close()

for each in ssTF_names:
    if not os.path.isfile('./data/{}_YEP/{}_Motif_1_bound.bed'.format(YEP_ID_dic[each.title().upper()], YEP_ID_dic[each.title().upper()])):
        print('No CX.bed for '+each.title())
        continue
    else:
        pass
    Tandem_gene_file = open('./result/Tandem_Genes_nuc.txt','r')
    outputfile = open('./result/Tandem_Genes_nuc_1.txt','a')
    for line in Tandem_gene_file:
        #print(Find_intermediate_nucleosome(line.split()[1], int(line.split()[9]), int(line.split()[3])))
        outputfile.write(line[:-1])
        for each_ in Find_intermediate_TFs(line.split()[0], int(float(line.split('\t')[6])), int(float(line.split('\t')[7])), YEP_ID_dic[each.title().upper()]):
            outputfile.write('{}_'.format(each.title())+each_+'\t')
        outputfile.write('\n')
    Tandem_gene_file.close()
    outputfile.close()
    os.unlink('./result/Tandem_Genes_nuc.txt')
    os.rename('./result/Tandem_Genes_nuc_1.txt','./result/Tandem_Genes_nuc.txt')
os.rename('./result/Tandem_Genes_nuc.txt','./result/Tandem_Genes_and_Stuff_inbetween.txt')