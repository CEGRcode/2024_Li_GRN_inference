import os, sys, re
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import random
import math
import copy
from scipy.integrate import solve_ivp
from distance_functions import *
from mutation_functions import *
from dynamics import *
from utility_functions import *
from population_functions import *
from GRN_Expanded_Combinatorial_2025 import GRN
from Combine_Redundant_Attractors import *
import sympy as sp
import getopt
from pathlib import Path

def read_refinement_files(dir_path='./result'):
    p = Path(dir_path)
    files = sorted(p.glob('Refinement_for_gene_*.txt'))
    out = {}
    for f in files:
        # try UTF-8, fall back to latin1 if needed
        try:
            text = f.read_text(encoding='utf-8')
        except UnicodeDecodeError:
            text = f.read_text(encoding='latin1')
        # extract gene name from filename (whatever is between prefix and .txt)
        m = re.match(r'Refinement_for_gene_(.+)\.txt$', f.name)
        key = m.group(1) if m else f.stem
        out[key] = text
    return out

if __name__ == '__main__':
    sys_path = os.getcwd()
    sys_input_RNAseq = ""
    sys_promoter_strengths = ""
    sys_gene_length = []
    sys_PerturbationPower = 0
    sys_iteration_num = 800
    sys_output_name = ""
    sys_random_seed = 42
    sys_training_count = 3000
    sys_random_test = 0
    sys_gene_to_focus = 0
    sys_refinement = 0

    #############################
    ### PARSE INPUT ARGUMENTS ###
    #############################
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "hr:i::t:n::p::l:o:k:e:j:m::c::f::",
                                   ["help", "input_RNAseq", "input_LG", "iteration_num", "promoter_strengths", "training_count", "PerturbationPower", "gene_length", "output_name", "randomseed", "test_GRN", "random_mode", "gene_to_focus", "refinement"])
    except getopt.GetoptError:
        print('Usage: GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py [-h] -r <input_RNAseq> [-i <iteration_num>]'
              ' -t <promoter_strengths> [-n <training_count>]'
              ' [-p <PerturbationPower>] [-c gene_to_focus]'
              ' -l <gene_length> -o <output_name> [-e <randomseed>] -j <test_GRN> [-m <random_mode>] [-f <refinement>]')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: GRN_Dynamic_Simulator_Combinatorial_Local_multistate_2025.py [-h] -r <input_RNAseq> [-i <iteration_num>]'
                ' -t <promoter_strengths> [-n <training_count>]'
                ' [-p <PerturbationPower>] [-c gene_to_focus]'
                ' -l <gene_length> -o <output_name> [-e <randomseed>] -j <test_GRN> [-m <random_mode>] [-f <refinement>]')
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
        elif opt in ("-i", "--iteration_num"):
            sys_iteration_num = arg
        elif opt in ("-f", "--refinement"):
            sys_refinement = arg
        elif opt in ("-c", "--gene_to_focus"):
            sys_gene_to_focus = int(arg)
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
        elif opt in ("-n", "--training_count"):
            sys_training_count = arg
        elif opt in ("-p", "--PerturbationPower"):
            sys_PerturbationPower = float(arg)
        elif opt in ("-m", "--random_mode"):
            if arg in ['True', 'TRUE', 1, '1']:
                sys_random_test = 1
            else:
                sys_random_test = 0
        elif opt in ("-o", "--output_name"):
            sys_output_name = arg
        elif opt in ("-j" "--test_GRN"):
            sys_test_GRN = sys_path + '/' + arg
            test_GRN = {}
            with open(sys_test_GRN, 'r') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    if parts:  # make sure the line isn't empty
                        key = parts[0]
                        value = parts[1:]
                        test_GRN[key] = value
        else:
            pass
    StableStatesCollector = {}
    TotalNumberOfGenes = len(sys_gene_length)
    WTTP = sys_WTTP
    TrainingCount = int(sys_training_count)

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
    indexes_of_same_gene = []
    for x in range(0, len(TranscriptionPofileMax)):
        if TranscriptionPofileMax[x] == TranscriptionPofileMin[x]:
            indexes_of_same_gene.append(x)
        else:
            indexes_of_diff_gene.append(x)
    Diff_gene_cannot_be_inferred = []  # Manually exclude differentially expressed genes if needed, the genes will be fixed to their ground truth values
    #indexes_of_diff_gene = [0, 1, 2, 5, 6, 7, 10, 15, 16, 18, 19, 22, 25, 29, 31, 33, 35, 37, 38, 39, 40, 41, 43, 44, 47, 48, 51, 53, 54, 56, 57, 58, 59, 60, 62, 65, 67, 68, 69, 72, 73, 75] # for consistency of transcriptional profile prediction test
    Num_of_genes_diff_TPM = len(indexes_of_diff_gene)
    ### Done Getting the genes having same TPM across all samples ###
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
        if sys_random_test:
            Configuration = String012ToMatrix(RandomString)
            LogicGates = LogicGatesString2Matrix_Expanded(','.join([','.join(map(str, range(TotalNumberOfGenes)))] * TotalNumberOfGenes))
            f0_c = eval(test_GRN['f0:'][0])
            #for f0_c_index in range(0, len(f0_c)):
            #    f0_c[f0_c_index][0] = np.random.rand()
            #    f0_c[f0_c_index][1] = np.random.rand()
        else:
            Configuration = String012ToMatrix(test_GRN['Adjacency Matrix:'][0])
            LogicGates = LogicGatesString2Matrix_Expanded(test_GRN['Logic Gate:'][0])
            f0_c = eval(test_GRN['f0:'][0])
        for each_index_of_same_gene in indexes_of_same_gene:
            Configuration[0][each_index_of_same_gene, :] = 0
            Configuration[1][each_index_of_same_gene, :] = 0
        ####################################################

        MutationRate = 0  # must be integer

        Sigmoid_k_init = eval(test_GRN['Hill Coefficient:'][0])
        
        #LogicGates = LogicGatesString2Matrix_Expanded('')
        TranscriptionThreshold = eval(test_GRN['TF Effective Threshold:'][0])
        
        DegradationRatemRNA = eval(test_GRN['mRNA Degradation Rate:'][0])
        Leakage = eval(test_GRN['Leakage Rate:'][0])
        TranscriptionRate = eval(test_GRN['Transcription Rate:'][0])
        return Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_k_init, Leakage, f0_c

    ########################################################################### Setting systematic Parameters #############################################################################
    (Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0_c) = ParametersInitiations(False)
    this_GRN = GRN('this_GRN', [], Configuration, MutationRate, TranscriptionRate, DegradationRatemRNA, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0_c)
    #######################################################################################################################################################################################
    
    # plug-in the refinement results
    if sys_refinement:
        refinement_data = read_refinement_files('./result')
        for each_refined_gene in refinement_data:
            refinement_AM = refinement_data[each_refined_gene].split(';')[1]
            refinement_LG = refinement_data[each_refined_gene].split(';')[2].split('\t')[0].split(',')
            refinement_LG = [int(x) for x in refinement_LG]
            refinement_f0_c = refinement_data[each_refined_gene].split(';')[3].split('\t')
            refinement_f0_c = [float(x) for x in refinement_f0_c]
            print(refinement_AM, refinement_LG, refinement_f0_c)
            
            this_GRN_configuration = ConfigurationTo012(this_GRN.Configuration)
            for configuration_index in range(0, len(this_GRN_configuration)):
                if configuration_index % TotalNumberOfGenes == int(each_refined_gene):
                    this_GRN_configuration = this_GRN_configuration[:configuration_index] + refinement_AM[configuration_index // TotalNumberOfGenes] + this_GRN_configuration[(configuration_index+1):] 
                else:
                    pass
            this_GRN.SetAM(String012ToMatrix(this_GRN_configuration))
            this_GRN.LogicGates[int(each_refined_gene)] = refinement_LG
            this_GRN.f0_c[int(each_refined_gene)] = refinement_f0_c
    else:
        pass
    
    outfile = open('./GRN_dynamic_local_search_result_{}.txt'.format(sys_output_name), 'a')
    outfile.write('indexes_of_diff_gene\t'+'\t'.join(map(str, indexes_of_diff_gene))+'\n')
    print('indexes_of_diff_gene: ', indexes_of_diff_gene)
    #print('Overexpression: ', Overexpression)
    GRN_Performance = []

    AM_for_the_gene = ''.join(map(str, [ch_ for index_, ch_ in enumerate(ConfigurationTo012(this_GRN.Configuration)) if index_ % TotalNumberOfGenes == sys_gene_to_focus]))
    TF_indices_for_the_gene = [i for i, ch in enumerate(AM_for_the_gene) if ch != '0']
    idx_map = {val: i for i, val in enumerate(indexes_of_diff_gene)}
    TF_indices_for_the_gene_in_diff_gene = [idx_map[x] for x in TF_indices_for_the_gene]
    #print('TF_indices_for_the_gene: ', TF_indices_for_the_gene, '\n')
    #print('TF_indices_for_the_gene_in_diff_gene: ', TF_indices_for_the_gene_in_diff_gene, '\n')
    #print('Configuration: ', AM_for_the_gene, '\n')
    #print('Logic gates: ', this_GRN.LogicGates[sys_gene_to_focus], '\n')
    #print('f0_c: ', this_GRN.f0_c[sys_gene_to_focus], '\n')
    #print('Sigmoid_k: ', this_GRN.Sigmoid_k[sys_gene_to_focus], '\n')
    #print('TR: ', this_GRN.TranscriptionRate[sys_gene_to_focus], '\n')
    #print('DR: ', this_GRN.DegradationRatemRNA[sys_gene_to_focus], '\n')

    for Global_i in range(0, len(WTTP)):

        # Setting up mRNA Protein numbers
        NewmRNA = np.array([0 for i in range(0, len(WTTP[str(Global_i)][1]))], dtype=float)

        for j in range(0, len(NewmRNA)):
            # The New mRNA list is generated by WTTP*Perturbation
            NewmRNA[j] = WTTP[str(Global_i)][1][j] * (1 + random.choice([1, -1]) * sys_PerturbationPower)

        #  Set overexpression or knockout
        this_GRN.SetmRNA(NewmRNA, WTTP[str(Global_i)][0], Overexpression[Global_i])

        # Run the model for a few times
        mRNACheckList = []
        scipystring = this_GRN.Delta_mRNA_Network(WTTP[str(Global_i)][0], Overexpression[Global_i], indexes_of_diff_gene, this_GRN.mRNA, Diff_gene_cannot_be_inferred)  # Calculate delta_mRNA
        #print(scipystring)
        exec(scipystring)
        #print('list(this_GRN.mRNA):', list(this_GRN.mRNA))
        transformed_KO_list = [i_ for i_, val in enumerate(indexes_of_diff_gene) if val in WTTP[str(Global_i)][0]]
        sol = solve_ivp(update_mRNA_protein,
                        [0, TrainingCount],
                        list(this_GRN.mRNA[indexes_of_diff_gene]),
                        args=([WTTP[str(Global_i)][0]]),
                        t_eval=[int(TrainingCount / 250) * tick for tick in range(0, 250 + 1)], method='RK45')
        npmRNA_continuous = sol.y.T

        print(Global_i)
        #print('AM: ', ConfigurationTo012(this_GRN.Configuration))
        #print('initial state:', list(np.round(this_GRN.mRNA[indexes_of_diff_gene], 1)), list(np.round(this_GRN.mRNA[indexes_of_diff_gene], 1)))
        #print('final state:  ', list(np.round(npmRNA_continuous[-1], 1)), list(np.round(npmRNA_continuous[-1], 1)))
        print('initial state:', list(np.round(this_GRN.mRNA, 1))[sys_gene_to_focus], list(np.round(this_GRN.mRNA[TF_indices_for_the_gene], 1)))
        print('final state:  ', list(np.round(npmRNA_continuous[-1], 1))[indexes_of_diff_gene.index(sys_gene_to_focus)], list(np.round(npmRNA_continuous[-1][TF_indices_for_the_gene_in_diff_gene], 1)))
        outfile.write('transcriptional profile #{}'.format(Global_i)+'\n')
        outfile.write('initial state\t'+'\t'.join(map(str, list(np.round(this_GRN.mRNA[indexes_of_diff_gene], 1))))+'\n')
        outfile.write('final state\t'+'\t'.join(map(str, list(np.round(npmRNA_continuous[-1], 1))))+'\n')

        # Use the eigenvalues of Jacobian matrix to determine if this is a fixed-point attractor
        y = sp.symbols('y1:%d' % (Num_of_genes_diff_TPM + 1))
        point = {}
        for npmRNA_i in range(0, Num_of_genes_diff_TPM):
            point[y[npmRNA_i]] = npmRNA_continuous[-1][npmRNA_i]
        Jacobian_matrix = [[] for Jm_i in range(0, Num_of_genes_diff_TPM)]
        Derivative_values = []
        for eq_i in range(0, len(scipystring.split('\n'))):
            eq_string = scipystring.split('\n')[eq_i].replace(" ", "")
            if eq_string[:11] == 'delta_mRNA[':
                exec('eq_of_string = {}'.format(eq_string.split('=')[1]))
                Derivative_values.append(eq_of_string.subs(point))
                for var_j in range(0, len(y)):
                    exec('df{}_dy{}_at_point = {}'.format(eq_i - 3, var_j, (sp.diff(eq_of_string, y[var_j])).subs(point)))
            else:
                continue
        for matrix_i in range(0, Num_of_genes_diff_TPM):
            for matrix_j in range(0, Num_of_genes_diff_TPM):
                Jacobian_matrix[matrix_i].append(eval('df{}_dy{}_at_point'.format(matrix_i, matrix_j)))
        eigenvalues = np.linalg.eigvals(Jacobian_matrix)
        Derivative_values = [abs(x) for x in Derivative_values]
        Derivative_values = [value for idx, value in enumerate(Derivative_values) if idx not in transformed_KO_list] # remove the genes that have been KOed.
        if max(eigenvalues) <= 0 and max(Derivative_values) < max(0.5, 0.01*max(npmRNA_continuous[-1])):
            IsPointAttractor = True
        else:
            IsPointAttractor = False
        Attractor_distance = np.round((1/len(indexes_of_diff_gene))*GetAttractorDistance(this_GRN.mRNA[indexes_of_diff_gene], npmRNA_continuous[-1], np.array(TranscriptionPofileMax)[indexes_of_diff_gene], np.array(TranscriptionPofileMin)[indexes_of_diff_gene]), 5)
        GRN_Performance.append(Attractor_distance)

        #print('Eigenvalues: ', eigenvalues)
        #print('Derivative_values: ', Derivative_values)
        print('IsPointAttractor: ', IsPointAttractor)
        outfile.write('is point attractor\t'+ str(IsPointAttractor) + '\n')
        print('Attractor distance: ', Attractor_distance, '\n')
        outfile.write('attractor distance\t'+ str(Attractor_distance) + '\n')
        #print('Current distance: ', abs(this_GRN.mRNA-npmRNA_continuous[-1])/TranscriptionPofileMax[indexes_of_diff_gene])
        ########################################################################################################################################################################################

        ############################################################################# Collect fixed-point attractor ############################################################################

    print('Overall performance: ', np.mean(GRN_Performance), sum(GRN_Performance))
    outfile.write('overall performance\t'+str(np.mean(GRN_Performance))+'\t'+str(sum(GRN_Performance))+'\n')
    outfile.close()