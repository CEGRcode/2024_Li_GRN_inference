import sys
import getopt
import os
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

def Multi_GenerateMutation(GRN_instances_dict, instance_key, sys_cache, j, indexes_of_diff_gene, sys_LG, sys_only_TF_DNA, sys_max_attempts):
    GRN_instances = GRN_instances_dict[instance_key]

    if GRN_instances.f0_c[indexes_of_diff_gene[j]][:2] in ([0, 0], [1, 1]):
        sys_cache_00 = copy.deepcopy(sys_cache)
        sys_cache_00[j] = {
            key[3:]: value
            for key, value in sys_cache[j].items()
            if key.startswith('00\t')
        }

        sys_cache_11 = copy.deepcopy(sys_cache)
        sys_cache_11[j] = {
            key[3:]: value
            for key, value in sys_cache[j].items()
            if key.startswith('11\t')
        }

        if GRN_instances.f0_c[indexes_of_diff_gene[j]][:2] == [0, 0]:
            sys_cache_list = [sys_cache_00, sys_cache_11]
        else:
            sys_cache_list = [sys_cache_11, sys_cache_00]

        GRN_instances.GenerateMutation(
            sys_cache_list[0],
            indexes_of_diff_gene,
            sys_LG,
            sys_only_TF_DNA,
            sys_max_attempts
        )

        if GRN_instances.MutationRate == 'All tested!':
            GRN_instances.SetMutationRate(1)
            if GRN_instances.f0_c[indexes_of_diff_gene[j]][:2] == [0, 0]:
                GRN_instances.Setf0([1, 1, 1, 1, 1, 1])
            else:
                GRN_instances.Setf0([0, 0, 1, 1, 1, 1])

            GRN_instances.GenerateMutation(
                sys_cache_list[1],
                indexes_of_diff_gene,
                sys_LG,
                sys_only_TF_DNA,
                sys_max_attempts
            )

    else:
        GRN_instances.GenerateMutation(
            sys_cache,
            indexes_of_diff_gene,
            sys_LG,
            sys_only_TF_DNA,
            sys_max_attempts
        )

    GRN_instances_dict[instance_key] = GRN_instances
    return

if __name__ == '__main__':
    np.seterr(over='ignore')
    os.makedirs("./result", exist_ok=True)

    sys_path = os.getcwd()
    sys_input_RNAseq = ""
    sys_input_RNAseq_2 = ""
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

    #############################
    ### PARSE INPUT ARGUMENTS ###
    #############################
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "hr:b::c::g::i::t:n::p::l:o:a::m::k:e:f::s::d::",
                                   ["help", "input_RNAseq", "input_RNAseq_2", "input_ChIP", "input_LG",
                                    "iteration_num", "promoter_strengths",
                                    "training_count", "PerturbationPower",
                                    "gene_length", "output_name", "max_attempts", "mRNA_elongation_rate", "checkpoint", "randomseed", "only_TFDNA", "specific_gene_to_train", "initial_GRN_configuration"])
    except getopt.GetoptError:
        print('Usage: EGRNM.py [-h] -r <input_RNAseq> [-b <input_RNAseq_2>] [-c <input_ChIP>] [-g <input_LG>] [-i <iteration_num>]'
              ' -t <promoter_strengths> [-n <training_count>]'
              ' [-p <PerturbationPower>]'
              ' -l <gene_length> -o <output_name> [-a <max_attempts>] [-m <mRNA_elongation_rate>] -k <checkpoint> [-e <randomseed>] [-f <only_TFDNA>] [-s <specific_gene_to_train>] [-d <initial_GRN_configuration>]')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: EGRNM.py [-h] -r <input_RNAseq> [-b <input_RNAseq_2>] [-c <input_ChIP>] [-g <input_LG>] [-i <iteration_num>]'
                ' -t <promoter_strengths> [-n <training_count>]'
                ' [-p <PerturbationPower>]'
                ' -l <gene_length> -o <output_name> [-a <max_attempts>] [-m <mRNA_elongation_rate>] -k <checkpoint> [-e <randomseed>] [-f <only_TFDNA>] [-s <specific_gene_to_train>] [-d <initial_GRN_configuration>]')
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
        elif opt in ("-b", "--input_RNAseq_2"):
            sys_input_RNAseq_2 = sys_path + '/' + arg
            sys_input_RNAseq_2 = np.array(pd.read_csv(
                sys_input_RNAseq_2, header=None, delimiter='\t', dtype=str))
            sys_WTTP_2 = {}
            for i in range(0, sys_input_RNAseq_2.shape[0]):
                # print(sys_input_RNAseq[i][0], type(sys_input_RNAseq[i][0]))
                if sys_input_RNAseq_2[i][0] == '-1':
                    sys_WTTP_2[str(i)] = [[], np.array(sys_input_RNAseq_2[i][1:], dtype=float)]
                elif len(sys_input_RNAseq_2[i][0]) == 1:
                    sys_WTTP_2[str(i)] = [[int(sys_input_RNAseq_2[i][0])], np.array(sys_input_RNAseq_2[i][1:], dtype=float)]
                else:
                    sys_WTTP_2[str(i)] = [np.array(sys_input_RNAseq_2[i][0].split(','), dtype=int).tolist(),
                                        np.array(sys_input_RNAseq_2[i][1:], dtype=float)]
        elif opt in ("-c", "--input_ChIP"):
            sys_input_ChIP = sys_path + '/' + arg
            sys_input_ChIP = np.array(pd.read_csv(sys_input_ChIP, header=None, delimiter='\t'), dtype=int)
        elif opt in ("-a", "--max_attempts"):
            sys_max_attempts = int(arg)
        elif opt in ("-m", "--mRNA_elongation_rate"):
            sys_mRNA_elongation_rate = float(arg)
        elif opt in ("-i", "--iteration_num"):
            sys_iteration_num = arg
        elif opt in ("-s", "--specific_gene_to_train"):
            sys_specific_gene_to_train = list(map(int, arg.split(',')))
        elif opt in ("-g", "--input_LG"):
            sys_LG = sys_path + '/' + arg
            sys_LG = open(sys_LG).readline().strip()
        elif opt in ("-d", "--initial_GRN_configuration"):
            sys_inital_AM = sys_path + '/' + arg
            sys_inital_AM = open(sys_inital_AM).readline().strip()
        elif opt in ("-e", "--randomseed"):
            sys_random_seed = int(arg)
        elif opt in ("-l", "--gene_length"):
            sys_gene_length = sys_path + '/' + arg
            sys_gene_length = np.array(pd.read_csv(
                sys_gene_length, header=None, delimiter='\t'), dtype=int)[0]
        elif opt in ("-t", "--promoter_strengths"):
            sys_promoter_strengths = sys_path + '/' + arg
            sys_promoter_strengths = np.array(pd.read_csv(
                sys_promoter_strengths, header=None, delimiter='\t'), dtype=float)
        elif opt in ("-k", "--checkpoint"):
            sys_checkpoint = float(arg)
        elif opt in ("-n", "--training_count"):
            sys_training_count = arg
        elif opt in ("-p", "--PerturbationPower"):
            sys_PerturbationPower = arg
        elif opt in ("-o", "--output_name"):
            sys_output_name = arg
        elif opt in ("-f", "--only_TFDNA"):
            if arg in ['1', 'TRUE', 'True']:
                sys_only_TF_DNA = True
            elif arg in ['0', 'FALSE', 'False']:
                sys_only_TF_DNA = False
        else:
            pass

    print('sys_only_TF_DNA: ', sys_only_TF_DNA)
    ###########################
    ### PREPARING MAIN LOOP ###
    ###########################
    Loopcounter = 0
    TotalNumberOfGenes = len(sys_gene_length)
    print(TotalNumberOfGenes, "TotalNumberOfGenes")
    if len(sys_LG) != 0:
        sys_LG = LogicGatesString2Matrix_Expanded(sys_LG).tolist()
    else:
        sys_LG = [[j for j in range(0, TotalNumberOfGenes)] for i in range(0, TotalNumberOfGenes)]
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

    outfile_check = open(sys_path + '/result/' + sys_output_name + '_CheckLoopCounter.txt'.format(TotalNumberOfGenes), 'a', buffering=1)
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

    BaseMatrixCollector_2 = []
    if type(sys_input_RNAseq_2) != str:
        for keys in sys_WTTP_2:
            for each in sys_WTTP_2[keys][0]:
                if each >= 0:
                    sys_WTTP_2[keys][1][each] = np.nan
                elif each <= -2:
                    sys_WTTP_2[keys][1][-each - 2] = np.nan
                else:
                    pass
            BaseMatrixCollector_2.append(sys_WTTP_2[keys][1])
        BaseMatrixCollector_2 = np.array(BaseMatrixCollector_2)

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
    print('indexes_of_diff_gene: ', indexes_of_diff_gene)
    
    def ParametersInitiations(MutationRate_, which_i, null_net, specific_gene):
        '''Parameters Initiations'''
        # Name
        # mRNA
        ####################################################
        if sys_inital_AM == '':
            RandomStringList = np.random.randint(2, size=(1, TotalNumberOfGenes ** 2))[0]
            RandomString = ''
            if null_net:
                for INT in RandomStringList:
                    RandomString = RandomString + '0'
            else:
                for INT in RandomStringList:
                    RandomString = RandomString + str(INT)

            Configuration = GetRegulatorForGene(RandomString, specific_gene)
        else:
            Configuration = GetRegulatorForGene(sys_inital_AM, specific_gene)
        for each_Ci in range(0, len(Configuration)):
            if each_Ci not in indexes_of_diff_gene:
                Configuration = Configuration[:each_Ci] + '0' + Configuration[each_Ci+1:]
            else:
                pass
        ####################################################

        MutationRate = 1  # must be integer

        TranscriptionRate = sys_promoter_strengths[which_i] * 60 * sys_mRNA_elongation_rate / sys_gene_length.tolist()
        TranscriptionRate = TranscriptionRate[specific_gene]
        LogicGates = sys_LG[specific_gene]

        #######################################################################################
        TranscriptionThreshold = []
        f0_c = []
        t1_t2_t3 = []
        for i in range(0, TotalNumberOfGenes):
            expression_levels = np.array(list(set(BaseMatrixCollector.T[i])))
            expression_levels = expression_levels[~np.isnan(expression_levels)]
            Result_t, Result_f0_c = get_param_by_TPM(np.array(expression_levels))
            TranscriptionThreshold_i = []
            for j in range(0, TotalNumberOfGenes):
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

        if type(sys_input_RNAseq_2) != str:
            expression_meas = np.array(list(BaseMatrixCollector.T[specific_gene]))
            expression_stds = np.array(list(BaseMatrixCollector_2.T[specific_gene]))
            Regulation_HA1 = '(1/(1+{}))*(X**n/(X**n+{}**n)+{}*X**n/(X**n+{}**n))'.format(f0_c[specific_gene][2], t1_t2_t3[specific_gene][0], f0_c[specific_gene][2], t1_t2_t3[specific_gene][1])
            Regulation_HA2 = '(1/(1+{}))*(X**n/(X**n+{}**n)+{}*X**n/(X**n+{}**n))'.format(f0_c[specific_gene][3], t1_t2_t3[specific_gene][0], f0_c[specific_gene][3], t1_t2_t3[specific_gene][1])
            Regulation_HR1 = '(1/(1+{}))*(X**n/(X**n+{}**n)+{}*X**n/(X**n+{}**n))'.format(f0_c[specific_gene][4], t1_t2_t3[specific_gene][0], f0_c[specific_gene][4], t1_t2_t3[specific_gene][2])
            Regulation_HR2 = '(1/(1+{}))*(X**n/(X**n+{}**n)+{}*X**n/(X**n+{}**n))'.format(f0_c[specific_gene][5], t1_t2_t3[specific_gene][0], f0_c[specific_gene][5], t1_t2_t3[specific_gene][2])
            Regulation_F = '{}+(1-{})*({})-({})*({})+({}+{}-1)*({})*({})'.format(f0_c[specific_gene][0], f0_c[specific_gene][0], Regulation_HA1, f0_c[specific_gene][0], Regulation_HR1, f0_c[specific_gene][0], f0_c[specific_gene][1], Regulation_HA2, Regulation_HR2)
            Sigmoid_k_init = 30
            mask = ~ (np.isnan(expression_meas) | np.isnan(expression_stds)) 
            expression_distributions = [list(t) for t in dict.fromkeys(zip(expression_meas[mask].tolist(), expression_stds[mask].tolist()))]
            #print('expression_distributions: ', expression_distributions)
            Sigmoid_k_init = min(30, fit_n_from_distribution_list_return_n(Regulation_F, expression_distributions))
        else:
            Sigmoid_k_init = 30
        print('Sigmoid_k_init: ', Sigmoid_k_init)
        print('Configuration: ', Configuration)
        DegradationRatemRNA = TranscriptionRate/TranscriptionPofileMax[specific_gene]
        Leakage = TranscriptionRate*TranscriptionPofileMin[specific_gene]/TranscriptionPofileMax[specific_gene]

        return Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_k_init, Leakage, f0_c

    GRN_List = []
    Out_List = []
    sys_cache = [{} for i in range(0, len(indexes_of_diff_gene))]
    for i in indexes_of_diff_gene: # the number of GRN and Out means the exact index of the specific gene
        GRN_List.append('GRN_{}'.format(i))
        Out_List.append('Out_{}'.format(i))

    for i in range(0, len(indexes_of_diff_gene)):
        if sys_checkpoint:
            with open(sys_path+'/data/GRN_checkpoint_{}_{}.pkl'.format(sys_output_name, i), 'rb') as f:
                exec('{} = pickle.load(f)'.format(GRN_List[i]))
            with open(sys_path+'/data/GRN_checkpoint_{}_{}.pkl'.format(sys_output_name, i), 'rb') as f:
                exec('{} = pickle.load(f)'.format(Out_List[i]))

            filename = "sys_cache_{}_gene_{}.txt".format(sys_output_name, indexes_of_diff_gene[i])
            filepath = os.path.join(sys_path, "data", filename)

            with open(filepath, "r") as f:
                header_line = f.readline().rstrip("\n")
                keys = header_line.split("|")
                data_dict = {k: [] for k in keys}
                for line in f:
                    row = line.rstrip("\n").split("|")
                    for idx_key, key in enumerate(keys):
                        if idx_key < len(row):
                            data_dict[key].append(row[idx_key])
                        else:
                            data_dict[key].append("")
                sys_cache[i] = copy.deepcopy(data_dict)

        else:
            (Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0_c) = ParametersInitiations(GlobalMutationRate, 0, False, indexes_of_diff_gene[i])
            exec('{} = GRN(\'{}\', {},{},{},{},{},{},{},{},{},{},{})'.format(GRN_List[i], GRN_List[i], '[]', 'Configuration', 'MutationRate', 'TranscriptionRate', 'DegradationRatemRNA', 'TranscriptionThreshold', 'LogicGates', 'Sigmoid_ks', 'Leakages', 'f0_c', 'sys_input_ChIP'))
            exec('{}=copy.deepcopy({})'.format(Out_List[i], GRN_List[i]))

    ###########################
    ### MAIN EVOLUTION LOOP ###
    ###########################
    
    FinishedGenes = []
    Fitness = [1e6 for j in range(0, len(indexes_of_diff_gene))]
    while Loopcounter <= int(sys_iteration_num):
        print('GRN1==>', GRN_1.Configuration)
        # Loop n times
        starting_time = time.time()
        DistanceOnLoop = [[] for m in range(0, len(WTTP))]

        for i in range(0, len(WTTP)):
            '''Re-initiate CurrentDistance'''
            CurrentDistance = copy.deepcopy(Fitness)

            '''Update promoter strength and dependent parameters'''
            TR = sys_promoter_strengths[i] * 60 * sys_mRNA_elongation_rate / sys_gene_length.tolist()
            for j in range(0, len(indexes_of_diff_gene)):
                if j in FinishedGenes:
                    continue
                else:
                    pass
                exec('{}.Update_Transcription_Rate({})'.format(GRN_List[j], TR[indexes_of_diff_gene[j]]))

            '''Setting up mRNA numbers'''
            NewmRNA = np.array([0 for i in range(0, len(WTTP[str(i)][1]))], dtype=float)
            for j in range(0, len(NewmRNA)):
                # The New mRNA list is generated by WTTP*Perturbation
                NewmRNA[j] = WTTP[str(i)][1][j] * (1 + PerturbationPower * (np.random.random_sample([1])[0] * np.random.choice([1, -1])))

            for j in range(0, len(indexes_of_diff_gene)):
                if j in FinishedGenes:
                    continue
                else:
                    pass
                exec('{}.SetmRNA({}, {}, {})'.format(GRN_List[j], 'NewmRNA', WTTP[str(i)][0], 'Overexpression[i]'))
            ########################## Recording GRNs that have been tested to avoid repetitive calculation ##########################
            for j in range(0, len(indexes_of_diff_gene)):
                if j in FinishedGenes:
                    continue
                else:
                    pass
                if eval('{}.f0_c[{}]'.format(GRN_List[j], indexes_of_diff_gene[j]))[:2] in [[0, 0], [1, 1]]:
                    temp_GRN_AM_LG = [eval('{}.Configuration'.format(GRN_List[j])), ''.join(map(str, eval('{}.f0_c[{}]'.format(GRN_List[j], indexes_of_diff_gene[j]))[:2]))+'\t'+",".join(map(str, eval('{}.LogicGates'.format(GRN_List[j]))))]
                else:
                    temp_GRN_AM_LG = [eval('{}.Configuration'.format(GRN_List[j])), ",".join(map(str, eval('{}.LogicGates'.format(GRN_List[j]))))]
                if temp_GRN_AM_LG[1] in sys_cache[j]:
                    if temp_GRN_AM_LG[0] not in sys_cache[j][temp_GRN_AM_LG[1]]:
                        sys_cache[j][temp_GRN_AM_LG[1]].append(temp_GRN_AM_LG[0])
                    else:
                        pass
                else:
                    sys_cache[j][temp_GRN_AM_LG[1]] = [temp_GRN_AM_LG[0]]
            ########################## Recording GRNs that have been tested to avoid repetitive calculation ##########################

            
            '''Run the dynamics model'''
            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            jobs = []

            for each_GRN in range(0, len(indexes_of_diff_gene)):
                if each_GRN in FinishedGenes:
                    continue
                else:
                    pass
                p = multiprocessing.Process(
                    target=Run_Dynamics,
                    args=(eval('{}'.format(GRN_List[each_GRN])), i, each_GRN, TrainingCount, WTTP, Overexpression[i], return_dict, indexes_of_diff_gene))
                jobs.append(p)
                p.start()

            for proc in jobs:
                proc.join()

            for j in range(0, len(indexes_of_diff_gene)):
                if j in FinishedGenes:
                    continue
                else:
                    pass
                IsPointAttractor = False
                npmRNA_continuous = np.array(return_dict[j])
                ########################## Using Jacobian eigenvalue to judge if this is a fixed-point attractor ###################
                y = sp.symbols('y1:%d' % (TotalNumberOfGenes + 1))
                #print('indexes_of_diff_gene: ', indexes_of_diff_gene[j])
                point = {}
                temp_scipystring = eval('{}.Delta_mRNA({}, {}, {}, analytical_expression=1)'.format(GRN_List[j], WTTP[str(i)][0], list(Overexpression[i]), list(WTTP[str(i)][1])))
                for npmRNA_i in range(0, TotalNumberOfGenes):
                    if npmRNA_i == indexes_of_diff_gene[j]:
                        point[y[npmRNA_i]] = npmRNA_continuous[-1]
                    else:
                        point[y[npmRNA_i]] = list(WTTP[str(i)][1])[npmRNA_i]
                #print('list(WTTP[str(i)][1])', list(WTTP[str(i)][1]))
                #print('point: ', point, '\n\n')
                #print('temp_scipystring: ', temp_scipystring)
                
                Derivative_values = []
                for eq_i in range(0, len(temp_scipystring.split('\n'))):
                    eq_string = temp_scipystring.split('\n')[eq_i].replace(" ", "")
                    #print('KO_List: ', WTTP[str(i)][0], 'npmRNA_continuous[-1]: ', list(np.round(npmRNA_continuous[-1])))
                    if eq_string[:10] == 'delta_mRNA':
                        temp_eq_string = eq_string.split('=')[1]
                        #print('temp_eq_string: ', temp_eq_string)
                        exec('eq_of_string = {}'.format(eq_string.split('=')[1]))
                        if npmRNA_continuous[-1]+eq_of_string.subs(point) < 0:
                            Derivative_values.append(-npmRNA_continuous[-1])
                        elif indexes_of_diff_gene[j] in WTTP[str(i)][0]:
                            Derivative_values.append(0)
                        else:
                            Derivative_values.append(eq_of_string.subs(point))
                    else:
                        continue
                if sp.diff(eq_of_string, y[indexes_of_diff_gene[j]]).subs(point) < 0 and (abs(Derivative_values[0]) < max(0.5, 0.01*max(npmRNA_continuous)) or sp.diff(eq_of_string, y[indexes_of_diff_gene[j]]).subs(point)+abs(Derivative_values[0])<0):
                    IsPointAttractor = True
                else:
                    IsPointAttractor = False
                    print('Gene {} Derivative: '.format(indexes_of_diff_gene[j]), abs(Derivative_values[0]), max(0.5, 0.01*max(npmRNA_continuous)), '2nd Derivative: ', sp.diff(eq_of_string, y[indexes_of_diff_gene[j]]).subs(point))
                #print('Derivative_values: ', Derivative_values)
                #print('sp.diff(eq_of_string, y[var_j]).subs(point): ', sp.diff(eq_of_string, y[indexes_of_diff_gene[j]]).subs(point))
                #print('*******************************************************************\n\n\n\n\n')
                #print('max Derivative: ', max(Derivative_values), 'index: ', Derivative_values.index(max(Derivative_values)))
                ########################## Using Jacobian eigenvalue to judge if this is a fixed-point attractor ###################

                if indexes_of_diff_gene[j] in WTTP[str(i)][0]:
                    CurrentDistance[j] = 0
                else:
                    if IsPointAttractor:
                        CurrentDistance[j] = abs(WTTP[str(i)][1][indexes_of_diff_gene[j]]-npmRNA_continuous[-1])/(TranscriptionPofileMax[indexes_of_diff_gene[j]]-TranscriptionPofileMin[indexes_of_diff_gene[j]])
                    else:
                        pass

            DistanceOnLoop[i] = CurrentDistance
            #print(i, '\'s CurrentDistance: ', CurrentDistance)
        print('DistanceOnLoop: ', np.mean(DistanceOnLoop, axis=0))
        #print('sys_cache', sys_cache[3], list(sys_input_ChIP[:, 3]))

        '''Check fitness and generate mutation.'''
        for j in range(0, len(indexes_of_diff_gene)):
            if j in FinishedGenes:
                continue
            else:
                pass
            if np.mean(DistanceOnLoop, axis=0)[j] < float(Fitness[j]):
                exec('{}=copy.deepcopy({})'.format(Out_List[j], GRN_List[j]))
                Fitness[j] = np.mean(DistanceOnLoop, axis=0)[j]
            else:
                exec('{}=copy.deepcopy({})'.format(GRN_List[j], Out_List[j]))
        print('GRN2==>', GRN_1.Configuration)
        manager = multiprocessing.Manager()
        shared_dict = manager.dict()
        jobs = []
        for j in range(0, len(indexes_of_diff_gene)):
            if j in FinishedGenes:
                continue
            else:
                pass
            shared_dict[GRN_List[j]] = eval(GRN_List[j])
            p = multiprocessing.Process(target=Multi_GenerateMutation, args=(shared_dict, GRN_List[j], sys_cache, j, indexes_of_diff_gene, sys_LG, sys_only_TF_DNA, sys_max_attempts))
            jobs.append(p)
            p.start()

        for proc in jobs:
            proc.join()

        for j in range(0, len(indexes_of_diff_gene)):
            if j in FinishedGenes:
                continue
            else:
                pass
            exec('{} = copy.deepcopy(shared_dict[\'{}\'])'.format(GRN_List[j], GRN_List[j]))
            if eval('{}.MutationRate'.format(GRN_List[j])) == 'All tested!' and j not in FinishedGenes:
                FinishedGenes.append(j)
                print('All tested!!!')
            else:
                pass
        print('GRN3==>', GRN_1.Configuration)

        FinishedGenes.extend([index_i for index_i in np.where(np.max(DistanceOnLoop, axis=0) < 1e-10)[0] if index_i not in FinishedGenes])
        print('FinishedGenes: ', FinishedGenes)

        with open(sys_path + '/result/' + sys_output_name + '_GRN_components.txt', 'w') as f:
            f.write('GeneIndex\tConfiguration\tLogicGates\tf0_c\tSigmoid_k\tLeakage\tDegradationRatemRNA\tTranscriptionRate\tTranscriptionThreshold\n')
            for j in range(0, len(indexes_of_diff_gene)):
                f.write(str(indexes_of_diff_gene[j])+'\t')
                f.write(str(eval('{}.Configuration'.format(Out_List[j])))+'\t')
                f.write(",".join(map(str, list(eval('{}.LogicGates'.format(Out_List[j]))) ))+'\t')
                f.write(",".join(map(str, list(eval('{}.f0_c'.format(Out_List[j])))[indexes_of_diff_gene[j]]))+'\t')
                f.write(str(eval('{}.Sigmoid_k'.format(Out_List[j])))+'\t')
                f.write(str(eval('{}.Leakage'.format(Out_List[j])))+'\t')
                f.write(str(eval('{}.DegradationRatemRNA'.format(Out_List[j])))+'\t')
                f.write(str(eval('{}.TranscriptionRate'.format(Out_List[j])))+'\t')
                f.write(str(np.array((eval('{}.TranscriptionThreshold'.format(Out_List[j]))[indexes_of_diff_gene[j]])).tolist())+'\n')

        outfile_check.write('{}\t{}\n'.format(Loopcounter, np.round(time.time() - starting_time, 3)))

        OutputFile3 = open(sys_path + '/result/' + sys_output_name + '_ResultFlow.txt', 'a')
        OutputFile3.write("\t".join(map(str, np.round(Fitness,5))) + "\n")
        OutputFile3.close()

        ### Stop if all genes have finished ###
        Loopcounter = Loopcounter + 1
        if len(FinishedGenes) == len(indexes_of_diff_gene):
            Loopcounter = int(sys_iteration_num) + 1
        else:
            pass
        ### Stop if all genes have finished ###

    outfile_check.close()
    ##################################
    ### END OF MAIN EVOLUTION LOOP ###
    ##################################

    ##################################
    ###  CONSTRUCT THE FINAL GRN   ###
    ##################################

    Final_AM = ''
    for j in range(0, TotalNumberOfGenes):
        for k in range(0, TotalNumberOfGenes):
            if k in indexes_of_diff_gene:
                Final_AM = Final_AM + str(eval('{}.Configuration'.format(Out_List[indexes_of_diff_gene.index(k)]))[j])
            else:
                Final_AM = Final_AM + '0'
    Final_LG = []
    Final_TranscriptionRate = []
    Final_Sigmoid_k = []
    Final_mRNA_DegradationRate = []
    Final_Leakage = []
    Final_f0_c = []
    Final_TranscriptionThreshold = []
    for j in range(0, TotalNumberOfGenes):
        if j in indexes_of_diff_gene:
            Final_LG.append(list(eval('{}.LogicGates'.format(Out_List[indexes_of_diff_gene.index(j)]))))
            Final_TranscriptionRate.append(float(eval('{}.TranscriptionRate'.format(Out_List[indexes_of_diff_gene.index(j)]))))
            Final_Sigmoid_k.append(float(eval('{}.Sigmoid_k'.format(Out_List[indexes_of_diff_gene.index(j)]))))
            Final_mRNA_DegradationRate.append(float(eval('{}.DegradationRatemRNA'.format(Out_List[indexes_of_diff_gene.index(j)]))))
            Final_Leakage.append(float(eval('{}.Leakage'.format(Out_List[indexes_of_diff_gene.index(j)]))))
            Final_f0_c.append(list(eval('{}.f0_c[{}]'.format(Out_List[indexes_of_diff_gene.index(j)], indexes_of_diff_gene.index(j)))))
            Final_TranscriptionThreshold.append(np.array(eval('{}.TranscriptionThreshold[{}]'.format(Out_List[indexes_of_diff_gene.index(j)], indexes_of_diff_gene.index(j)))).tolist())
        else:
            Final_LG.append([unit for unit in range(0, TotalNumberOfGenes)])
            Final_TranscriptionRate.append(0)
            Final_Sigmoid_k.append(0)
            Final_mRNA_DegradationRate.append(0)
            Final_Leakage.append(0)
            Final_f0_c.append([0,0,0,0,0,0])
            Final_TranscriptionThreshold.append([[0,0,0] for _ in range(0, TotalNumberOfGenes)])

    ### Pickle the GRN instances ###
    for j in range(0, len(Out_List)):
        with open(sys_path+'/result/GRN_checkpoint_{}_{}.pkl'.format(sys_output_name, j), 'wb') as f:
            exec('pickle.dump({}, f)'.format(Out_List[j]))
    ### Pickle the GRN instances ###

    print(
        'Adjacency Matrix:\t{}\nLogic Gate:\t{}\nf0_c:\t{}\nHill Coefficient:\t{}\nTranscription Rate:\t{}\nmRNA Degradation Rate:\t{}\nLeakage Rate:\t{}\nTF Effective Threshold:\t{}\n\n'.format(
            Final_AM,
            Matrix2String01_LG_Expanded(Final_LG),
            Final_f0_c,
            Final_Sigmoid_k,
            Final_TranscriptionRate,
            Final_mRNA_DegradationRate,
            Final_Leakage,
            Final_TranscriptionThreshold))

    #######################################
    ### SAVE ADDITIONAL OUTPUT TO FILES ###
    #######################################
    OutputFile = open(sys_path + '/result/' + 'BetaTest_G{}_Result.txt'.format(TotalNumberOfGenes), 'a')
    OutputFile.write('{}\t{}\t{}\n'.format(Final_AM, Matrix2String01_LG_Expanded(Final_LG), np.mean(Fitness)))
    OutputFile.close()

    OutputFile = open(sys_path + '/result/' + 'BetaTest_G{}_Result_f0.txt'.format(TotalNumberOfGenes), 'a')
    OutputFile.write(
        'Adjacency Matrix:\t{}\nLogic Gate:\t{}\nf0:\t{}\nHill Coefficient:\t{}\nTranscription Rate:\t{}\nmRNA Degradation Rate:\t{}\nLeakage Rate:\t{}\nTF Effective Threshold:\t{}\n\n'.format(
            Final_AM,
            Matrix2String01_LG_Expanded(Final_LG),
            Final_f0_c,
            Final_Sigmoid_k,
            Final_TranscriptionRate,
            Final_mRNA_DegradationRate,
            Final_Leakage,
            Final_TranscriptionThreshold))
    OutputFile.close()

    for j in range(0, len(indexes_of_diff_gene)):
        keys = sorted(sys_cache[j].keys())
        with open(sys_path + '/result/' + "sys_cache_{}_gene_{}.txt".format(sys_output_name, indexes_of_diff_gene[j]), "w") as f:
            # Write header row (column names)
            header = "|".join(map(str, keys))
            f.write(header + "\n")
            
            # Use zip_longest to iterate over the lists row-wise.
            for row in itertools.zip_longest(*(sys_cache[j][k] for k in keys), fillvalue=""):
                line = "|".join(map(str, row))
                f.write(line + "\n")
