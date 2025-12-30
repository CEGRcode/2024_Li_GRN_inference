import os
import sys
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

sys_path = os.getcwd()
sys_promoter_strengths = ""
sys_gene_length = []
sys_PerturbationPower = 0
sys_iteration_num = 800
sys_output_name = ""
sys_random_seed = 42
sys_training_count = 3000
sys_random_test = 0

#############################
### PARSE INPUT ARGUMENTS ###
#############################
try:
    opts, args = getopt.getopt(sys.argv[1:],
                                "hi::t:n::p::l:o:k:e:j:m::",
                                ["help", "input_LG", "iteration_num", "promoter_strengths", "training_count", "PerturbationPower", "gene_length", "output_name", "randomseed", "test_GRN", "random_mode"])
except getopt.GetoptError:
    print('Usage: GRN_Dynamic_Simulator_Combinatorial_Global_multistate_2025.py [-h] [-i <iteration_num>]'
            ' -t <promoter_strengths> [-n <training_count>]'
            ' [-p <PerturbationPower>]'
            ' -l <gene_length> -o <output_name> [-e <randomseed>] -j <test_GRN> [-m <random_mode>]')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print('Usage: GRN_Dynamic_Simulator_Combinatorial_Global_multistate_2025.py [-h] [-i <iteration_num>]'
            ' -t <promoter_strengths> [-n <training_count>]'
            ' [-p <PerturbationPower>]'
            ' -l <gene_length> -o <output_name> [-e <randomseed>] -j <test_GRN> [-m <random_mode>]')
        sys.exit()
    elif opt in ("-i", "--iteration_num"):
        sys_iteration_num = arg
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
    
########################################################################### Setting systematic Parameters #############################################################################
TrainingCount = sys_training_count
sys_PerturbationPower = sys_PerturbationPower
sys_random_seed = 42
TotalNumberOfGenes = len(sys_gene_length)
#######################################################################################################################################################################################
StableStatesCollector = {}

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
        for f0_c_index in range(0, len(f0_c)):
            f0_c[f0_c_index][0] = np.random.rand()
            f0_c[f0_c_index][1] = np.random.rand()
    else:
        Configuration = String012ToMatrix(test_GRN['Adjacency Matrix:'][0])
        LogicGates = LogicGatesString2Matrix_Expanded(test_GRN['Logic Gate:'][0])
        f0_c = eval(test_GRN['f0:'][0])

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
indexes_of_diff_gene = [i for i, x in enumerate(DegradationRatemRNA) if x != 0]
Num_of_genes_diff_TPM = len(indexes_of_diff_gene)
#######################################################################################################################################################################################
print('indexes_of_diff_gene: ', indexes_of_diff_gene)

max_values = np.array([TranscriptionRate[i]/(DegradationRatemRNA[i]+1e-99) for i in range(0, TotalNumberOfGenes)])
min_values = np.array([Leakages[i]/(DegradationRatemRNA[i]+1e-99) for i in range(0, TotalNumberOfGenes)])
unique_counts = np.array([2 for i in range(0, TotalNumberOfGenes)])
#unique_counts = df.nunique()
# Example usage:
total_iterations = 8000
range_limits = [(min_values.tolist()[i], max_values.tolist()[i]) for i in range(0, TotalNumberOfGenes)]
proportions = list(((unique_counts-1)/sum(unique_counts-1)))  # Example: more resolution on the first dimension
#print('range_limits: ', range_limits)
#print('proportions: ', proportions)
combinations = generate_limited_combinations(total_iterations, range_limits, proportions)
#print('combinations: ', combinations)
# Convert combinations to WTTPs:
WTTP = {}
KO_list = []
for combo_i in range(0, len(combinations)):
    WTTP[str(combo_i)] = [KO_list, np.array(combinations[combo_i])]
Overexpression = np.zeros((len(WTTP), TotalNumberOfGenes))

for Global_i in range(0, len(WTTP)):

    '''Setting up mRNA Protein numbers'''
    NewmRNA = np.array([0 for i in range(0, len(WTTP[str(Global_i)][1]))], dtype=float)

    for j in range(0, len(NewmRNA)):
        # The New mRNA list is generated by WTTP*Perturbation
        NewmRNA[j] = WTTP[str(Global_i)][1][j] * (1 + random.choice([1, -1]) * sys_PerturbationPower)

    #  Set overexpression or knockout
    this_GRN.SetmRNA(NewmRNA, WTTP[str(Global_i)][0], Overexpression[Global_i])

    # Run the model for a few times
    mRNACheckList = []
    scipystring = this_GRN.Delta_mRNA_Network(WTTP[str(Global_i)][0], Overexpression[Global_i], indexes_of_diff_gene, this_GRN.mRNA)  # Calculate delta_mRNA
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
    #print('initial state:', list(np.round(this_GRN.mRNA[indexes_of_diff_gene],1)))
    #print('final state:  ', list(np.round(npmRNA_continuous[-1],1)[indexes_of_diff_gene]))

    '''Use the eigenvalues of Jacobian matrix to determine if this is a fixed-point attractor'''
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
    #print('Eigenvalues: ', eigenvalues)
    #print('Derivative_values: ', Derivative_values)
    #print('IsPointAttractor: ', IsPointAttractor)
    #print('\n')
    ########################################################################################################################################################################################

    ############################################################################# Collect fixed-point attractor ############################################################################
    mRNA_list_for_compare = []
    for each in npmRNA_continuous[-1]:
        mRNA_list_for_compare.append(int(each))
        #print('mRNA_list_for_compare: ', mRNA_list_for_compare)

    if IsPointAttractor == True:
        #print(str(mRNA_list_for_compare))
        if str(mRNA_list_for_compare) in StableStatesCollector:
            StableStatesCollector[str(mRNA_list_for_compare)] = StableStatesCollector[str(mRNA_list_for_compare)] + 1
        else:
            StableStatesCollector[str(mRNA_list_for_compare)] = 1
    elif IsPointAttractor == False:
        pass
    else:
        raise Exception('Error!')

#print(StableStatesCollector)

StableStates = []
Frequency_Threshold = 0
for keys in StableStatesCollector:
    print(keys)
    if StableStatesCollector[keys] <= Frequency_Threshold:
        pass
    else:
        StableStates.append(eval(keys))

print('Total Length:', len(StableStates))

Combined = Combine_Redundant_Attractors(arrays = np.array(StableStates), cutoff=0.05)

Attractor_Distance_matrix = np.zeros((len(Combined), len(Combined)))
for novel_i in range(0, len(Combined)):
    for novel_j in range(0, len(Combined)):
        Attractor_Distance_matrix[novel_i][novel_j] = GetAttractorDistance(Combined[novel_i],Combined[novel_j],np.max(Combined, axis=0),np.min(Combined, axis=0))
Unique_attractor_N = (1 + np.count_nonzero(np.mean(Attractor_Distance_matrix, axis=0) > 1))
Combined = np.round(Combined,1).tolist()

# Insert 0s for non-diff genes.
for i in range(0, len(Combined)):
    temp_list = [0]*TotalNumberOfGenes
    for j in range(0, TotalNumberOfGenes):
        if j in indexes_of_diff_gene:
            temp_list[j] = Combined[i][indexes_of_diff_gene.index(j)]
        else:
            pass
    Combined[i] = temp_list

print('Unique: ', Unique_attractor_N)
print(Combined)

'''
with open('/Users/rl884/Downloads/Test_attractors_99.txt', 'w') as file:
    for row in Combined:
        file.write("\t".join(map(str, row)) + "\n")'''
