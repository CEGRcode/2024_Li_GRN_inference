import sys
import getopt
import os
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

sys_path = os.getcwd()
sys_input_RNAseq = ""
sys_input_ChIP = ""
sys_input_Path = ""
sys_training_count = 2500
sys_PerturbationPower = 0
sys_output_name = ""
Order_of_genes = ''
num_chunks = 1
sys_start_AM_i = 0
sys_start_AM = ''

try:
    opts, args = getopt.getopt(sys.argv[1:],
                                "hr:c::i:n::p::o:g:k:",
                                ["help", "input_RNAseq", "input_ChIP", "input_path",
                                "training_count", "PerturbationPower",
                                "output_name", "gene_names", "start_point"])
except getopt.GetoptError:
    print('Usage: remove_dispensable_edges.py [-h] -r <input_RNAseq> [-c <input_ChIP>] -i <input_path> '
            ' [-n <training_count>]'
            ' [-p <PerturbationPower>]'
            ' -o <output_name> -g <gene_names> -s <start_point>')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print('Usage: remove_dispensable_edges.py [-h] -r <input_RNAseq> [-c <input_ChIP>] -i <input_path> '
            ' [-n <training_count>]'
            ' [-p <PerturbationPower>]'
            ' -o <output_name> -g <gene_names> -s <start_point>')
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
    elif opt in ("-s", "--start_point"):
        sys_start_input = open(sys_path + '/' + arg, 'r')
        last_line = None
        for last_line in sys_start_input:
            pass
        sys_start_AM_i = int(last_line.split()[0])
        sys_start_AM = last_line.split()[1]
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
print('Number of profiles: ', len(WTTP))
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
#print('sys_input_ChIP: ', sys_input_ChIP)

Fitness_Record = []
if sys_start_AM == '':
    pass
else:
    Final_AM = sys_start_AM

for AM_i in range(sys_start_AM_i, len(Final_AM)):
    print('AM_i: ', AM_i, flush=True)
    Attractor_Distance_for_One_Loop = 0

    if AM_i == sys_start_AM_i:
        this_GRN.SetAM(String012ToMatrix(Final_AM))
    elif Final_AM[AM_i-1] == '0' or [item for sublist in sys_input_ChIP for item in sublist][AM_i-1] == 1:
        Fitness_Record.append(1e99)
        continue
    else:
        this_GRN.SetAM(String012ToMatrix(Final_AM[:AM_i-1]+'0'+Final_AM[AM_i:]))

    #print('AM: ', ConfigurationTo012(this_GRN.Configuration), flush=True)
    for Global_i in range(0, len(WTTP)):
        print('AM_i: ', AM_i, '\t', 'Global_i: ', Global_i, flush=True)
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
        transformed_KO_list = [i_ for i_, val in enumerate(indexes_of_diff_gene) if val in WTTP[str(Global_i)][0]]
        #print('list(this_GRN.mRNA):', list(this_GRN.mRNA))
        sol = solve_ivp(update_mRNA_protein,
                        [0, TrainingCount],
                        list(this_GRN.mRNA[indexes_of_diff_gene]),
                        args=([WTTP[str(Global_i)][0]]),
                        t_eval=[int(TrainingCount / 250) * tick for tick in range(0, 250 + 1)], method='RK45')

        npmRNA_continuous = sol.y.T

        #print(Global_i, flush=True)
        #print('initial state:', list(np.round(this_GRN.mRNA, 1)[indexes_of_diff_gene]), flush=True)
        #print('final state:  ', list(np.round(npmRNA_continuous[-1], 1)), flush=True)

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
        #print('Attractor distance: ', np.round((1/len(indexes_of_diff_gene))*GetAttractorDistance(this_GRN.mRNA[indexes_of_diff_gene], npmRNA_continuous[-1], np.array(TranscriptionPofileMax)[indexes_of_diff_gene]), 5), '\n')
        ########################################################################################################################################################################################

        ############################################################################# Collect fixed-point attractor ############################################################################
        if IsPointAttractor == True:
            Attractor_Distance_for_One_Loop = Attractor_Distance_for_One_Loop + np.round((1/len(indexes_of_diff_gene))*GetAttractorDistance(this_GRN.mRNA[indexes_of_diff_gene], npmRNA_continuous[-1], np.array(TranscriptionPofileMax)[indexes_of_diff_gene], np.array(TranscriptionPofileMin)[indexes_of_diff_gene]), 5)
        elif IsPointAttractor == False:
            Attractor_Distance_for_One_Loop = Attractor_Distance_for_One_Loop + TotalNumberOfGenes * 5
        else:
            raise Exception('Error!')

    Fitness_Record.append(Attractor_Distance_for_One_Loop)
    #print('Fitness_Record: ', Fitness_Record, flush=True)
    #print('\n\n')
    if Attractor_Distance_for_One_Loop <= Fitness_Record[0] and AM_i != 0:
        Final_AM = Final_AM[:AM_i-1] + '0' + Final_AM[AM_i:]
    else:
        pass

    f = open('./result/temp_remove_dispensible_edge.txt', 'w')
    print(str(AM_i)+'\t'+Final_AM, file=f, flush=True)
    f.close()