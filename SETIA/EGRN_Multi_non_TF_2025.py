import sys
import getopt
import os
import path_setup
import pandas as pd
import numpy as np
import random
import copy
import multiprocessing
import time
import pickle
from scipy.integrate import solve_ivp
from dynamics import *
from distance_functions import *
from GRN_Expanded_Combinatorial_2025 import GRN
from mutation_functions import *
from population_functions import *
from utility_functions import *
from matplotlib import pyplot as plt
import itertools

def pick_sublists(lists, x, seed=None):
    """
    If len(lists) > x, randomly pick x sublists.
    Otherwise, return all sublists.
    """
    if seed is not None:
        random.seed(seed)

    if len(lists) <= x:
        return lists.copy()

    return random.sample(lists, x)

def split_into_chunks_get_group(nums, N, i):
    """
    Split nums into N contiguous chunks as evenly as possible (by length),
    and return the i-th chunk (0-based).
    """
    if N <= 0:
        raise ValueError("N must be positive")
    if i < 0 or i >= N:
        raise IndexError("i out of range")

    L = len(nums)
    base = L // N
    extra = L % N  # first 'extra' groups get one more element

    start = i * base + min(i, extra)
    size = base + (1 if i < extra else 0)
    return nums[start:start + size]


def get_possible_AM(lst, choices=(0,1,2)):
    # find positions of 1s
    ones_idx = [i for i, v in enumerate(lst) if v == 1]
    if not ones_idx:
        return [lst.copy()]   # no ones -> just return the original list
    results = []
    for combo in itertools.product(choices, repeat=len(ones_idx)):
        new = lst.copy()
        for idx, val in zip(ones_idx, combo):
            new[idx] = val
        results.append(new)
    return results

if __name__ == '__main__':
    np.seterr(over='ignore')

    sys_path = os.getcwd()
    sys_input_RNAseq = ""
    sys_input_ChIP = ""
    sys_promoter_strengths = ""
    sys_mRNA_elongation_rate = 4.8
    sys_gene_length = []
    sys_training_count = 250
    sys_PerturbationPower = 0
    sys_iteration_num = 800
    sys_output_name = ""
    sys_LG = ''
    sys_checkpoint = 0
    sys_random_seed = 42
    sys_only_TF_DNA = False
    sys_max_attempts = 1000
    sys_inital_AM = ''
    sys_specific_gene_to_train = []
    sys_TF = []
    sys_all_genes = []
    sys_chunck_i = None
    sys_AM_trial_num = 100

    #############################
    ### PARSE INPUT ARGUMENTS ###
    #############################
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "hr:c::g::i::t:n::l:o:a::m::k:f::s::x:y:b:j:",
                                   ["help", "input_RNAseq", "input_ChIP", "input_LG",
                                    "iteration_num",
                                    "training_count",
                                    "gene_length", "output_name", "max_attempts", "mRNA_elongation_rate", "specific_gene_to_train", "TF", "all_genes", "upper_AM", "chunk"])
    except getopt.GetoptError:
        print('Usage: EGRN_Multi_non_TF_2025.py [-h] -r <input_RNAseq> [-c <input_ChIP>] [-g <input_LG>] [-i <iteration_num>]'
              '  [-n <training_count>]'
              ' -l <gene_length> -o <output_name> [-a <max_attempts>] [-m <mRNA_elongation_rate>] [-s <specific_gene_to_train>] -x <TF> -y <all_genes> -b <upper_AM> -j <chunk>'
              'Example: python EGRN_Multi_non_TF_2025.py -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt -n 3000 -c data/GRN_ssTFs_Sc_TF_DNA.txt -x ./data/GRN_ssTFs_TF_set.txt -y ./data/GRN_ssTFs_column_names.txt -g ./data/GRN_ssTFs_Sc_LG.txt')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: EGRN_Multi_non_TF_2025.py [-h] -r <input_RNAseq> [-c <input_ChIP>] [-g <input_LG>] [-i <iteration_num>]'
                '  [-n <training_count>]'
                ' -l <gene_length> -o <output_name> [-a <max_attempts>] [-m <mRNA_elongation_rate>] [-s <specific_gene_to_train>] -x <TF> -y <all_genes> -b <upper_AM> -j <chunk>'
                'Example: python EGRN_Multi_non_TF_2025.py -r data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt -n 3000 -c data/GRN_ssTFs_Sc_TF_DNA.txt -x ./data/GRN_ssTFs_TF_set.txt -y ./data/GRN_ssTFs_column_names.txt -g ./data/GRN_ssTFs_Sc_LG.txt')
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
        elif opt in ("-a", "--max_attempts"):
            sys_max_attempts = int(arg)
        elif opt in ("-m", "--mRNA_elongation_rate"):
            sys_mRNA_elongation_rate = float(arg)
        elif opt in ("-j", "--chunk"):
            sys_chunck_i = int(arg)
        elif opt in ("-b", "--upper_AM"):
            sys_AM_trial_num = int(arg)
        elif opt in ("-i", "--iteration_num"):
            sys_iteration_num = arg
        elif opt in ("-s", "--specific_gene_to_train"):
            sys_specific_gene_to_train = list(map(int, arg.split(',')))
        elif opt in ("-g", "--input_LG"):
            sys_LG = sys_path + '/' + arg
        elif opt in ("-l", "--gene_length"):
            sys_gene_length = sys_path + '/' + arg
            sys_gene_length = np.array(pd.read_csv(
                sys_gene_length, header=None, delimiter='\t'), dtype=int)[0]
        elif opt in ("-n", "--training_count"):
            sys_training_count = arg
        elif opt in ("-o", "--output_name"):
            sys_output_name = arg
        elif opt in ("-x", "--TF"):
            sys_TF = sys_path + '/' + arg
            sys_TF = np.array(pd.read_csv(
                sys_TF, header=None, delimiter='\t'), dtype=str)[0]
            if sys_TF[-1] == 'nan':
                sys_TF = sys_TF[:-1]
            else:
                pass
        elif opt in ("-y", "--all_genes"):
            sys_all_genes = sys_path + '/' + arg
            sys_all_genes = np.array(pd.read_csv(
                sys_all_genes, header=None, delimiter='\t'), dtype=str)[0]
            if sys_all_genes[-1] == 'nan':
                sys_all_genes = sys_all_genes[:-1]
            else:
                pass
        else:
            pass

    print('sys_TF: ', len(sys_TF), '\n')
    TF_pos = {v:i for i,v in enumerate(sys_all_genes)}
    TF_indices_in_all_genes = [TF_pos.get(x, -1) for x in sys_TF]
    print('TF_indices_in_all_genes: ', TF_indices_in_all_genes, '\n')
    print('sys_all_genes: ', len(sys_all_genes), '\n')

    ###########################
    ### PREPARING MAIN LOOP ###
    ###########################
    Loopcounter = 0
    TotalNumberOfGenes = len(sys_all_genes)
    minDistance = []
    meanDistance = []
    TrainingCount = int(sys_training_count)
    Target_Distance = 0
    # Recombination_Frequency = 1
    SelectionPower = 0.25
    Natural_Selection_Memory = 0
    TimeToMakeALongJump = 0
    InitialGlobalMutationRate = 10
    GlobalMutationRate = InitialGlobalMutationRate
    PerturbationPower = float(sys_PerturbationPower)
    np.random.seed(sys_random_seed)
    random.seed(sys_random_seed)


    '''Population proportions and TranscriptionProfile/mRNAList initiation'''
    InitialTranscriptionProfile = []
    WTTP = copy.deepcopy(sys_WTTP)
    #weights_of_WTTP = WTTP_weights(copy.deepcopy(WTTP))
    BaseMatrixCollector = []
    TranscriptionPofileMax = []
    TranscriptionPofileMin = []
    TranscriptionPofileAve = []
    Overexpression = np.zeros((len(sys_WTTP), TotalNumberOfGenes))
    Knockout = np.ones((len(sys_WTTP), TotalNumberOfGenes))
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
    Num_of_genes_diff_TPM = TotalNumberOfGenes - len(indexes_of_diff_gene)
    ### Done Getting the genes having same TPM across all samples ###
    if sys_specific_gene_to_train == []:
        pass
    else:
        indexes_of_diff_gene = sys_specific_gene_to_train
    #print('indexes_of_diff_gene: ', indexes_of_diff_gene)
    diff_TFs_index = sorted(list(set(TF_indices_in_all_genes) & set(indexes_of_diff_gene)))
    print('diff_TFs_index: ', diff_TFs_index)


    def ParametersInitiations(AM_for_the_gene, LG_for_the_gene, the_gene_i, diff_TFs_index):
        '''Parameters Initiations'''
        # Name
        # mRNA
        ####################################################
        Configuration = AM_for_the_gene
        ####################################################

        MutationRate = 1  # must be integer

        TranscriptionRate = 10

        Sigmoid_k_init = 30

        LogicGates = LG_for_the_gene

        #######################################################################################
        TranscriptionThreshold = []
        f0_c = []
        t1_t2_t3 = []
        for i in [gene_i] + diff_TFs_index:
            expression_levels = np.array(list(set(BaseMatrixCollector.T[i])))
            expression_levels = expression_levels[~np.isnan(expression_levels)]
            Result_t, Result_f0_c = get_param_by_TPM(np.array(expression_levels))
            TranscriptionThreshold_i = []
            for j in [gene_i] + diff_TFs_index:
                expression_level_j = np.array(list(set(BaseMatrixCollector.T[j])))
                TranscriptionThreshold_ij = []
                if len(expression_level_j) - sum(1 for x in expression_level_j if isinstance(x, float) and math.isnan(x)) == 1:
                    if sum(1 for x in BaseMatrixCollector.T[j] if isinstance(x, float) and math.isnan(x)) in [0, 2]:
                        TranscriptionThreshold_ij = [expression_level_j[~np.isnan(expression_level_j)][0] for duplicate_i in range(0, 3)]
                    elif sum(1 for x in BaseMatrixCollector.T[j] if isinstance(x, float) and math.isnan(x)) == 1:
                        if 0 in Knockout[:, j]:
                            TranscriptionThreshold_ij = [1e-10, 1e-10, 1e-10]
                        else:
                            TranscriptionThreshold_ij = [expression_level_j[~np.isnan(expression_level_j)][0] for duplicate_i in range(0, 3)]
                    else:
                        raise Exception('Illegal number of np.nan!')
                else:
                    expression_level_j = expression_level_j[~np.isnan(expression_level_j)]                    
                    for each_t in Result_t:
                        if each_t == 0.5:
                            TranscriptionThreshold_ij.append(0.5*(TranscriptionPofileMax[j]+TranscriptionPofileMin[j]))
                        else:
                            average_list = split_and_average(sorted(expression_level_j))
                            if each_t == 0.33:
                                TranscriptionThreshold_ij.append(average_list[0])
                            elif each_t == 0.66:
                                TranscriptionThreshold_ij.append(average_list[1])
                TranscriptionThreshold_i.append(np.array(TranscriptionThreshold_ij))
            TranscriptionThreshold.append(TranscriptionThreshold_i)

            f0_c.append(Result_f0_c)
            t1_t2_t3.append(Result_t)
        #######################################################################################

        DegradationRatemRNA = 10/TranscriptionPofileMax[the_gene_i]
        Leakage = TranscriptionRate*TranscriptionPofileMin[the_gene_i]/TranscriptionPofileMax[the_gene_i]

        return Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_k_init, Leakage, f0_c, t1_t2_t3

    ###########################
    ### MAIN EVOLUTION LOOP ###
    ###########################

    #print('WTTP: ', WTTP)
    output_nodes = []
    output_edges_1_solid = []
    output_edges_1_dashed = []
    output_edges_2_solid = []
    output_edges_2_dashed = []

    indexes_of_diff_gene_to_do = split_into_chunks_get_group(indexes_of_diff_gene, 50, sys_chunck_i)
    print('len(indexes_of_diff_gene): ', len(indexes_of_diff_gene), 'this run', len(indexes_of_diff_gene_to_do), flush=True)

    for gene_i in indexes_of_diff_gene_to_do:
        TF_DNA_gene_i = list(sys_input_ChIP[:, gene_i])
        TF_DNA_gene_i_diff_TFs_index = [TF_DNA_gene_i[x] for x in diff_TFs_index]
        t_ticks = [int(TrainingCount/250)*tick for tick in range(0, 250+1)]
        if gene_i in TF_indices_in_all_genes or all(x == 0 for x in TF_DNA_gene_i_diff_TFs_index):
            continue
        else:
            print('gene_i: ', gene_i, sys_all_genes[gene_i], '\n')
            #print('TF_DNA_gene_i_diff_TFs_index: ', TF_DNA_gene_i_diff_TFs_index)
            TF_index_for_the_TFs_binding_to_this_gene = [diff_TFs_index[TF_index] for TF_index in [TF_i for TF_i, TF_v in enumerate(TF_DNA_gene_i_diff_TFs_index) if TF_v == 1]]
            #print('sys_all_genes index for the TFs binding to this gene: ', TF_index_for_the_TFs_binding_to_this_gene)
            #print('Names for the TFs binding to this gene: ', [sys_all_genes[real_index] for real_index in TF_index_for_the_TFs_binding_to_this_gene])
            #print('Expression: ', list(BaseMatrixCollector.T[gene_i]))
            expression_levels = np.array(list(set(BaseMatrixCollector.T[gene_i])))
            expression_levels = expression_levels[~np.isnan(expression_levels)]
            #print('Expression: ', expression_levels)
            LG_for_gene_i = get_LG_for_gene(sys_LG, gene_i*len(sys_all_genes), count=len(sys_all_genes))
            #print('LG_for_gene_i: ', LG_for_gene_i)
            LG_of_TF_for_gene_i = [LG_for_gene_i[x] for x in diff_TFs_index]
            #print('LG_for_gene_i: ', LG_of_TF_for_gene_i)
            LG_of_TF_for_gene_i = [int({v:i for i,v in enumerate(sorted(set(map(int,LG_of_TF_for_gene_i))))}[int(x)]) for x in LG_of_TF_for_gene_i]
            LG_of_TF_for_gene_i = [0] + [x + 1 for x in LG_of_TF_for_gene_i]
            #print('LG_for_gene_i: ', LG_of_TF_for_gene_i)

            CurrentFitness = {}
            Possible_AM = pick_sublists(get_possible_AM(TF_DNA_gene_i_diff_TFs_index), sys_AM_trial_num)
            for each_AM in Possible_AM:
                
                (Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0_c, t1_t2_t3) = ParametersInitiations(each_AM, LG_of_TF_for_gene_i, gene_i, diff_TFs_index)
                Configuration = [0] + Configuration
                exec('GRN_for_the_gene = GRN(\'GRN_0\', {},{},{},{},{},{},{},{},{},{},{})'.format('[]', 'Configuration', 'MutationRate', 'TranscriptionRate', 'DegradationRatemRNA', 'TranscriptionThreshold', 'LogicGates', 'Sigmoid_ks', 'Leakages', 'f0_c', 'sys_input_ChIP'))
                GRN_for_the_gene.SetAM(''.join(map(str, Configuration)))


                CurrentDistance = []
                for i in range(0, len(WTTP)):

                    idx_map = {val: idx for idx, val in enumerate(diff_TFs_index)}
                    temp_KO_list = [idx_map[x] for x in WTTP[str(i)][0] if x in idx_map]
                    for temp_KO_list_i in range(0, len(temp_KO_list)):
                        if temp_KO_list[temp_KO_list_i] == -1:
                            pass
                        elif temp_KO_list[temp_KO_list_i] >= 0:
                            temp_KO_list[temp_KO_list_i] = temp_KO_list[temp_KO_list_i] + 1
                        elif temp_KO_list[temp_KO_list_i] < -1:
                            temp_KO_list[temp_KO_list_i] = temp_KO_list[temp_KO_list_i] -1
                        else:
                            pass
                    GRN_for_the_gene.SetmRNA(WTTP[str(i)][1][[gene_i]+diff_TFs_index], temp_KO_list, Overexpression[i])
                    #print('diff_TFs_index_expression: ', [gene_i]+diff_TFs_index, list(GRN_for_the_gene.mRNA))
                    
                    scipystring = GRN_for_the_gene.Delta_mRNA(temp_KO_list, Overexpression[i], list(GRN_for_the_gene.mRNA))
                    TrainningCount_func = sys_training_count
                    scipystring = (scipystring
                                + '\nInitialState = [GRN_for_the_gene.mRNA[0]]'
                                + '\nsol = solve_ivp(update_mRNA_protein, [0, TrainningCount_func], InitialState, '
                                +                    'args=([temp_KO_list]), t_eval=t_ticks, method=\'RK45\')')
                    scipystring = scipystring + '\nGRN_for_the_gene.sol = sol.y'
                    #print(GRN_for_the_gene.Configuration)
                    #print(scipystring)
                    exec(scipystring)
                    sol = GRN_for_the_gene.sol
                    GRN_for_the_gene.sol = []
                    sol = [int(x.item()) for x in list(np.round(sol.T, 2))]
                    
                    y = sp.symbols('y')
                    point = {}
                    point[y] = sol[-1]
                    Derivative_values = []
                    for eq_i in range(0, len(scipystring.split('\n'))):
                        eq_string = scipystring.split('\n')[eq_i].replace(" ", "")
                        #print('KO_List: ', WTTP[str(i)][0], 'npmRNA_continuous[-1]: ', list(np.round(npmRNA_continuous[-1])))
                        if eq_string[:10] == 'delta_mRNA':
                            temp_eq_string = eq_string.split('=')[1]
                            #print('temp_eq_string: ', temp_eq_string)
                            exec('eq_of_string = {}'.format(eq_string.split('=')[1]))
                            if sol[-1]+eq_of_string.subs(point) < 0:
                                Derivative_values.append(-sol[-1])
                            elif gene_i in WTTP[str(i)][0]:
                                Derivative_values.append(0)
                            else:
                                Derivative_values.append(eq_of_string.subs(point))
                        else:
                            pass

                    if sp.diff(eq_of_string, y).subs(point) < 0 and (abs(Derivative_values[0]) < max(0.5, 0.01*max(sol)) or sp.diff(eq_of_string, y).subs(point)+abs(Derivative_values[0])<0):
                        IsPointAttractor = True
                        CurrentDistance.append(abs(GRN_for_the_gene.mRNA[0] - sol[-1]))
                    else:
                        IsPointAttractor = False
                        CurrentDistance.append(1e10)
                    #print('IsPointAttractor: ', IsPointAttractor)
                    #print('\n\n')
                
                CurrentFitness[''.join(map(str, each_AM))] = sum(CurrentDistance)

            best_AM = min(CurrentFitness, key=CurrentFitness.get)
            outfile = open('./result/nonTF_genes_fit_{}.txt'.format(sys_chunck_i), 'a')
            outfile.write(str(gene_i)+'\t'+str((CurrentFitness[best_AM]/len(WTTP))/(max(expression_levels)-min(expression_levels)))+'\n')
            outfile.flush()
            outfile.close()
            #print('Gene', gene_i, ' AM: ', CurrentFitness, 'Distance:', CurrentFitness[best_AM], '\n', flush=True)
            output_nodes.append(gene_i)

            #print('best_AM: ', best_AM, '\n')

            for best_AM_i in range(1, len(best_AM)):
                if best_AM[best_AM_i] == '1':
                    if len(list(set([LG_i for LG_i, v in enumerate(LG_of_TF_for_gene_i) if v == LG_of_TF_for_gene_i[best_AM_i]]) & set([AM_activator_i for AM_activator_i, ch in enumerate(best_AM) if ch == '1']))) > 1:
                        output_edges_1_dashed.append([diff_TFs_index[best_AM_i], gene_i, LG_of_TF_for_gene_i[best_AM_i]])
                    else:
                        output_edges_1_solid.append([diff_TFs_index[best_AM_i], gene_i])
                elif best_AM[best_AM_i] == '2':
                    if len(list(set([LG_i for LG_i, v in enumerate(LG_of_TF_for_gene_i) if v == LG_of_TF_for_gene_i[best_AM_i]]) & set([AM_inhibitor_i for AM_inhibitor_i, ch in enumerate(best_AM) if ch == '2']))) > 1:
                        output_edges_2_dashed.append([diff_TFs_index[best_AM_i], gene_i, LG_of_TF_for_gene_i[best_AM_i]])
                    else:
                        output_edges_2_solid.append([diff_TFs_index[best_AM_i], gene_i])
    
    json_ = ''
    json_ = json_ + '{\n  "nodes": [\n'
    for node_index in output_nodes:
        (Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0_c, t1_t2_t3) = ParametersInitiations(GRN_for_the_gene.Configuration, GRN_for_the_gene.LogicGates, node_index, diff_TFs_index)
        json_ = json_ + '\t{\n\t  ' + '"id": "{}",\n'.format(sys_all_genes[node_index]) + '\t  "label": "{}",\n'.format(sys_all_genes[node_index]) + '\t  "TR": {},\n'.format(10) + '\t  "f0": {},\n'.format(f0_c[0][0]) + '\t  "f0p": {},\n'.format(f0_c[0][1])  + '\t  "t1": {},\n'.format(t1_t2_t3[0][0]) + '\t  "t2": {},\n'.format(t1_t2_t3[0][1]) + '\t  "t3": {},\n'.format(t1_t2_t3[0][2]) + '\t  "c1": {},\n'.format(f0_c[0][2]) + '\t  "c2": {},\n'.format(f0_c[0][3]) + '\t  "c3": {},\n'.format(f0_c[0][4]) + '\t  "c4": {},\n'.format(f0_c[0][5]) + '\t  "k": {},\n'.format(Sigmoid_ks) + '\t  "Lk": {},\n'.format(Leakages) + '\t  "Deg": {}'.format(DegradationRatemRNA) + '\n  \t},\n'
    json_ = json_[:-2]
                
    json_ = json_ + '\n  ],\n  "edges": [\n'

    for each_edge in output_edges_1_solid:
        json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(sys_all_genes[each_edge[0]]) + '\t  "target": "{}",\n'.format(sys_all_genes[each_edge[1]]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"triangle"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
    for each_edge in output_edges_2_solid:
        json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(sys_all_genes[each_edge[0]]) + '\t  "target": "{}",\n'.format(sys_all_genes[each_edge[1]]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"tee"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
    for each_edge in output_edges_1_dashed:
        json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(sys_all_genes[each_edge[0]]) + '\t  "target": "{}",\n'.format(sys_all_genes[each_edge[1]]) + '\t  "label": "{}",\n'.format(sys_all_genes[each_edge[2]]) + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
    for each_edge in output_edges_2_dashed:
        json_ = json_ + '\t{\n\t  ' + '"source": "{}",\n'.format(sys_all_genes[each_edge[0]]) + '\t  "target": "{}",\n'.format(sys_all_genes[each_edge[1]]) + '\t  "label": "{}",\n'.format(sys_all_genes[each_edge[2]]) + '\t  "style": [\n\t\t"dashed",\n\t\t"tee"\n\t,\n\t\t"#00A86B"\n\t  ]'+ '\n  \t},\n'
    json_ = json_[:-2]
    json_ = json_ + '\n  ]\n}'

    outfile = open('./result'+'/GRN_nonTF_{}.json'.format(sys_chunck_i), 'a')
    outfile.write(json_)
    outfile.close()
