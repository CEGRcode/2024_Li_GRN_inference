import os
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
from GRN import GRN

sys_path = os.getcwd()

########################################################################### Setting systematic Parameters #############################################################################
sys_input_RNAseq = ""
sys_input_ChIP = ""
sys_promoter_strengths = ""
sys_protein_degradation_rate = ""
sys_mRNA_elongation_rate = 4.8
sys_aa_elongation_rate = 8
sys_gene_length = []
sys_samplesize = 1
sys_training_count = 250
sys_event = {100:{'mRNA': {3: 3000}}, 200:{'DegradationRatemRNA': {6: 50}}, 300:{}} # variables are mRNA, Protein, TranscriptionRate, TranslationRate, DegradationRatemRNA, DegradationRateProtein
sys_PerturbationPower = 0.1
sys_iteration_num = 800
sys_output_name = ""
#######################################################################################################################################################################################

############################################################################ Data from an existing model ##############################################################################
sys_input_RNAseq = np.array(pd.read_csv(
    '/Users/ruihaoli/Downloads/in_silico_GRN7/Alpha_Attractors_7.txt',
    header=None, delimiter='\t').iloc[1:, 1:], dtype=float)
#print(sys_input_RNAseq)
sys_WTTP = {}
for i in range(0, sys_input_RNAseq.shape[0]):
    if sys_input_RNAseq[i][0] == -1:
        sys_WTTP[str(i)] = [[], sys_input_RNAseq[i][1:]]
    else:
        sys_WTTP[str(i)] = [[int(sys_input_RNAseq[i][0])], sys_input_RNAseq[i][1:]]

sys_gene_length = np.array(pd.read_csv(
    '/Users/ruihaoli/Downloads/in_silico_GRN7/Alpha_GeneLength_7.txt',
    header=None, delimiter='\t')[1:], dtype=int)[0]

sys_promoter_strengths = np.array(pd.read_csv(
    '/Users/ruihaoli/Downloads/in_silico_GRN7/Alpha_PromoterStrength_7.txt',
    header=None, delimiter='\t')[1:], dtype=float)

sys_protein_degradation_rate = [(0.00796) * 60 for i in range(0, len(sys_gene_length))]
#######################################################################################################################################################################################

########################################################################### Setting training Parameters ###############################################################################
Loopcounter = 0
SampleSize = sys_samplesize
TotalNumberOfGenes = len(sys_gene_length)
minDistance = []
meanDistance = []
Total_Time_Span = max(sys_event.keys())
minDistance_Single = 1.01e6
Target_Distance = 0
# Recombination_Frequency = 1
SelectionPower = 0.25
Natural_Selection_Memory = 0
TimeToMakeALongJump = 0
InitialGlobalMutationRate = 2
GlobalMutationRate = InitialGlobalMutationRate
PerturbationPower = sys_PerturbationPower
GlobalMutationDirection_LG = InitiateGlobalMutationDirection_LG(TotalNumberOfGenes)
#######################################################################################################################################################################################

##################################################### Setting up population proportions for each individual GRN instance ##############################################################
################################################## Calculating the minimal, average, and maximal expression levels of a gene ##########################################################
InitialProportions = 1 / SampleSize
InitialTranscriptionProfile = []
WTTP = copy.deepcopy(sys_WTTP)
BaseMatrixCollector = []
TranscriptionPofileMax = []
TranscriptionPofileMin = []
TranscriptionPofileAve = []
Overexpression = np.zeros((len(sys_WTTP), TotalNumberOfGenes))
Overexpression[:] = np.nan
'''if overexpressed, the maximum mRNA/protein level shouldn't be used in parameter estimation. So, we store the overexpression mRNA/protein levels
in Overexpression and mask those values in the mRNA/protein vectors. The indicator will be -2 if first gene is overexpressed and -3 if second gene.'''
Knockout = np.zeros((len(sys_WTTP), TotalNumberOfGenes))
Knockout[:] = np.nan
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

for x in range(0, SampleSize):
    InitialTranscriptionProfile.append([float(random.randrange(1, 100)) for i in range(0, TotalNumberOfGenes)])
np.array(InitialTranscriptionProfile)
#######################################################################################################################################################################################

############################################################################ Scaling promoter Strengths ###############################################################################
dt_scaler = (max(TranscriptionPofileMax / (np.min(sys_promoter_strengths, axis=0) * 60 * sys_mRNA_elongation_rate / sys_gene_length)) / 8)
print('dt_scaler: ', dt_scaler)
sys_promoter_strengths = sys_promoter_strengths * dt_scaler
#######################################################################################################################################################################################

############################################################################ Setting up model parameters ##############################################################################
def ParametersInitiations():
    '''Parameters initiations'''
    # Name
    # mRNA
    Protein = [0.0 for i in range(0, TotalNumberOfGenes)]

    # Configuration = String012ToMatrix('0112002121110112020202200')
    Configuration = String012ToMatrix('0000000000000001102002001021000010010000100001001')
    ########################################################################################################
    '''RandomStringList = np.random.randint(3,size = (1, TotalNumberOfGenes**2))[0]
    RandomString = ''
    for INT in RandomStringList:
        RandomString = RandomString + str(INT)
    Configuration = String012ToMatrix(RandomString)'''
    ########################################################################################################

    if type(sys_input_ChIP) == str:
        pass
    else:
        for ChIP_i in range(0, sys_input_ChIP.shape[0]):
            for ChIP_j in range(0, sys_input_ChIP.shape[1]):
                if sys_input_ChIP[ChIP_i][ChIP_j] == 0:
                    Configuration[0][ChIP_i][ChIP_j] = 0
                    Configuration[1][ChIP_i][ChIP_j] = 0
                else:
                    pass

    # Need to make sure the activators/repressors proportion is ~10%, and the activators and repressors do not present simutaneously.

    MutationRate = 2  # must be integer

    #TranscriptionRate = sys_promoter_strengths[0]*60*sys_mRNA_elongation_rate/sys_gene_length.tolist()
    TranscriptionRate = ((np.max(sys_promoter_strengths, axis=0) * 60 * sys_mRNA_elongation_rate) / sys_gene_length).tolist()

    TranslationRate = (sys_aa_elongation_rate * 60 / (sys_gene_length / 3)).tolist()

    DegradationRatemRNA = [0 for i in range(0, TotalNumberOfGenes)]
    for i in range(0, len(DegradationRatemRNA)):
        DegradationRatemRNA[i] = TranscriptionRate[i] / TranscriptionPofileMax[i]

    DegradationRateProtein = sys_protein_degradation_rate

    DilutionRate = [0.0 for i in range(0, TotalNumberOfGenes)]

    # TranscriptionThreshold = np.random.rand(TotalNumberOfGenes, TotalNumberOfGenes)
    '''TranscriptionThreshold = np.random.randint(6,size=(1, TotalNumberOfGenes, TotalNumberOfGenes))[0]'''
    TranscriptionThreshold = np.array(
        [[5.0 for i in range(0, TotalNumberOfGenes)] for i in range(0, TotalNumberOfGenes)])
    for i in range(0, TranscriptionThreshold.shape[0]):
        for j in range(0, TranscriptionThreshold.shape[1]):
            TranscriptionThreshold[i][j] = (TranscriptionPofileAve[i] * TranslationRate[i]) / DegradationRateProtein[i]

    Sigmoid_k_init = []
    for i in range(0, TotalNumberOfGenes):
        Sigmoid_k_init.append(np.log(10 ** -3 / (1 - 10 ** -3)) / np.log(TranscriptionThreshold[i][0] / (
                    TranscriptionPofileMax[i] * TranslationRate[i] / DegradationRateProtein[i])))
        # Sigmoid_k_init.append(np.log((1-10**-3)/10**-3)/np.log(TranscriptionThreshold[i][0]/(TranscriptionPofileMin[i]*TranslationRate[i]/DegradationRateProtein[i])))

    LogicGates = LogicGatesString2Matrix('00000010000010')
    # LogicGates = [np.random.randint(2, size=2).tolist() for i in range(0, TotalNumberOfGenes)]

    Leakage = []
    # PromoterMinList = ((np.min(sys_promoter_strengths, axis=0)*60*sys_mRNA_elongation_rate)/sys_gene_length).tolist()
    for i in range(0, TotalNumberOfGenes):
        # Leakage.append(min(TranscriptionPofileMin[i]*DegradationRatemRNA[i], PromoterMinList[i]))
        Leakage.append(TranscriptionPofileMin[i] * DegradationRatemRNA[i])
    # f0 = np.random.uniform(0,1,TotalNumberOfGenes)
    # f0 = [0.2260398,  0.72477163 ,0.18961756, 0.53594393, 0.24441876]
    # f0 = []
    # for i in range(0, TotalNumberOfGenes):
    #    f0.append((DegradationRatemRNA[i]*TranscriptionPofileAve[i]-Leakage[i])/(TranscriptionRate[i]-Leakage[i]))
    f0 = [0.1 for f0i in range(0, TotalNumberOfGenes)]
    # print(f0)
    return (Protein, Configuration, MutationRate, TranscriptionRate, TranslationRate, DegradationRatemRNA,
            DegradationRateProtein, DilutionRate, TranscriptionThreshold, LogicGates, Sigmoid_k_init, Leakage, f0)


GRN_List = []
Out_List = []
for i in range(0, SampleSize):
    GRN_List.append('GRN_{}'.format(i + 1))
    Out_List.append('Out_{}'.format(i + 1))

for i in range(0, SampleSize):
    Protein, Configuration, MutationRate, TranscriptionRate, TranslationRate, DegradationRatemRNA, DegradationRateProtein, DilutionRate, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0 = ParametersInitiations()
    exec('{} = GRN(\'{}\', {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})'.format(GRN_List[i], GRN_List[i],
                                                                                 'InitialTranscriptionProfile[i]',
                                                                                 'Protein', 'Configuration',
                                                                                 'InitialProportions', 'MutationRate',
                                                                                 'TranscriptionRate', 'TranslationRate',
                                                                                 'DegradationRatemRNA',
                                                                                 'DegradationRateProtein',
                                                                                 'DilutionRate',
                                                                                 'TranscriptionThreshold', 'LogicGates',
                                                                                 'Sigmoid_ks', 'Leakages', 'f0', 'sys_input_ChIP'))
#######################################################################################################################################################################################

############################################################################ Setting up utility parameters ############################################################################
Proportions = []
for i in range(0, SampleSize):
    Proportions.append(eval('{}.Proportion'.format(GRN_List[i])))

DistanceToShow = []
ProportionToShow = []
mRNAMonitor = [[] for i in range(0, TotalNumberOfGenes)]
ProteinMonitor = [[] for i in range(0, TotalNumberOfGenes)]
# DistanceOnLoop = [[] for m in range(0, len(WTTP))]
Loopcounter = 0
minDistance = []
meanDistance = []
minDistance_Single = 1.01e6
StableStatesCollector = {}
RunTime_ = [1000, 1000, 1600, 1800, 600, 1050]
Matched_Number = 0
Calculate_AttractorDistance = []

Scaler = []
for i in range(0, TotalNumberOfGenes):
    Scaler.append([4000, 0])
Scaler = np.array(Scaler)
#######################################################################################################################################################################################

############################################################################ Test example #######################################################################################################
# WTTP = {'0':[[4], np.array([603.8380791,202.7548367,20.45564933,46.73061065,0,0.382960406,4.433457027])], '1':[[2], [25.22811239,195.8378067,0,24.29728989,195.094524,4.856189477,18.59902051
# ]], '2':[[5], np.array([589.9956061,384.4656287,776.7735358,36.09745911,34.73911524,0,123.0572403])]}
#################################################################################################################################################################################################


for Global_i in range(0, len(WTTP)):
    fig = plt.figure(figsize=(10, 10))
    for Number_of_initial_positions in range(0, 1):
        ################################## Updates the rate of translation and promoter strengths and make sure they satisfy the steady state assumption ######################################
        '''for j in range(0, SampleSize):
            # exec('{}.SetmRNA({},{})'.format(GRN_List[j], 'NewmRNA', knocklist[Global_i]))
            exec('{}.SetTranslationRate({})'.format(GRN_List[j], 'TranslationRate'))
            exec('{}.Update_Promoter_Strength({},{},{},{})'.format(GRN_List[j], Global_i, 'sys_promoter_strengths', 'sys_mRNA_elongation_rate', 'sys_gene_length'))'''
        #######################################################################################################################################################################################

        ########################################################################## Setting up initial mRNA and proteins #######################################################################
        NewmRNA = np.array([0 for i in range(0, len(WTTP[str(Global_i)][1]))], dtype=float)
        NewProtein = np.array([0 for i in range(0, len(WTTP[str(Global_i)][1]))], dtype=float)
        for j in range(0, len(NewmRNA)):
            # The New mRNA list is generated by WTTP*Perturbation
            NewmRNA[j] = WTTP[str(Global_i)][1][j] * (
                        1 + PerturbationPower * (np.random.random_sample([1])[0] * np.random.choice([1, -1])))
        # print(NewmRNA)
        for j in range(0, len(NewmRNA)):
            NewProtein[j] = eval(
                '({}.TranslationRate[j]*NewmRNA[j])/({}.DegradationRateProtein[j]+{}.DilutionRate[j])'.format(
                    GRN_List[0], GRN_List[0], GRN_List[0]))
        #######################################################################################################################################################################################

        ############################################################################## Setting up overexpression ##############################################################################
        for j in range(0, SampleSize):
            # exec('{}.SetmRNA({},{})'.format(GRN_List[j], 'NewmRNA', knocklist[Global_i]))
            exec('{}.SetmRNA({},{},Overexpression[Global_i])'.format(GRN_List[0], 'NewmRNA', WTTP[str(Global_i)][0]))
            # exec('{}.SetProtein({},{})'.format(GRN_List[j], 'NewProtein', knocklist[Global_i]))
            exec('{}.SetProtein({},{},Overexpression[Global_i])'.format(GRN_List[0], 'NewProtein', WTTP[str(Global_i)][0]))
        ########################################################################################################################################################################################

        ########################################################################## Solving the differential equations ##########################################################################
        npmRNA_continuous = []
        npmRNA = []
        temp_previous_timepoint = 0
        previous_state = list(NewmRNA) + list(NewProtein)

        ########################################################################## Parameter re-initiation for each GRN #####################################################################
        for j in range(0, SampleSize):
            Protein, Configuration, MutationRate, TranscriptionRate, TranslationRate, DegradationRatemRNA, DegradationRateProtein, DilutionRate, TranscriptionThreshold, LogicGates, Sigmoid_ks, Leakages, f0 = ParametersInitiations()
            exec('{} = GRN(\'{}\', {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})'.format(GRN_List[j], GRN_List[j],'InitialTranscriptionProfile[j]','Protein', 'Configuration','InitialProportions','MutationRate','TranscriptionRate','TranslationRate','DegradationRatemRNA','DegradationRateProtein','DilutionRate','TranscriptionThreshold','LogicGates','Sigmoid_ks', 'Leakages', 'f0','sys_input_ChIP'))
        #####################################################################################################################################################################################

        for events_timepoints in sys_event:
            time_interval = events_timepoints - temp_previous_timepoint
            # solve the DEs separately in the time intervals of each event
            scipystring = eval('{}.Delta_mRNA(WTTP[str(Global_i)][0], Overexpression[Global_i])'.format(GRN_List[0]))  # Construct the differential equations
            #print(events_timepoints, '/n', scipystring)
            exec(scipystring)
            sol = solve_ivp(update_mRNA_protein, [0, time_interval], previous_state,
                            args=([WTTP[str(Global_i)][0]]), t_eval=[tick for tick in range(0, time_interval)])
            if len(npmRNA_continuous) == len(npmRNA) == 0:
                npmRNA_continuous = sol.y[:TotalNumberOfGenes].T
                npmRNA = np.round(sol.y[:TotalNumberOfGenes]).T
            else:
                npmRNA_continuous = np.concatenate((npmRNA_continuous, sol.y[:TotalNumberOfGenes].T), axis=0)
                npmRNA = np.concatenate((npmRNA, np.round(sol.y[:TotalNumberOfGenes]).T), axis=0)
            exec('{}.SetmRNA({},{},Overexpression[Global_i])'.format(GRN_List[0], 'sol.y.T[-1][:TotalNumberOfGenes]', WTTP[str(Global_i)][0]))
            exec('{}.SetProtein({},{},Overexpression[Global_i])'.format(GRN_List[0], 'sol.y.T[-1][TotalNumberOfGenes:]', WTTP[str(Global_i)][0]))
            ###################################################################### Update the GRN variables according to the event ##############################################################
            for keys_ in sys_event[events_timepoints]:
                # keys_ are the variables that need to be altered, such as TranscriptionRate, mRNA, and Protein
                for indices in sys_event[events_timepoints][keys_]:
                    # indices are the genes whose corresponding variables need to be altered
                    exec('{}.{}[{}] = max(0, {}.{}[{}] + {})'.format(GRN_List[0], keys_, indices, GRN_List[0], keys_, indices, sys_event[events_timepoints][keys_][indices]))
            #####################################################################################################################################################################################
            temp_previous_timepoint = events_timepoints
            previous_state = list(eval('{}.mRNA'.format(GRN_List[0]))) + list(eval('{}.Protein'.format(GRN_List[0])))
            #print('previous_state: {}'.format(previous_state), '\n', 'npmRNA_continuous: {}'.format(npmRNA_continuous[-1]))
        #npmRNA_continuous = np.array(npmRNA_continuous)
        #npmRNA = np.array(npmRNA)
        ########################################################################################################################################################################################

        ############################################################################### Determine attractor type ###############################################################################
        IsPointAttractor = True
        # print(np.std(npmRNA[150:251,4]))
        for o in range(0, TotalNumberOfGenes):
            if np.std(npmRNA[npmRNA.shape[0] - 100:npmRNA.shape[0] + 1, o]) < 1:
                continue
            else:
                arithmetic_sequence = np.argwhere(npmRNA[Total_Time_Span - 100:Total_Time_Span + 1, o] == np.amax(
                    npmRNA[Total_Time_Span - 100:Total_Time_Span + 1, o]))
                diff_seq = []
                for i in range(0, len(arithmetic_sequence) - 1):
                    diff_seq.append(arithmetic_sequence[i + 1][0] - arithmetic_sequence[i][0])
                if (np.std(diff_seq) < 1) and (len(diff_seq) >= 2):
                    continue
                else:
                    IsPointAttractor = False
        #print('IsPointAttractor: ', IsPointAttractor)
        ########################################################################################################################################################################################

        ############################################################################# Collect fixed-point attractor ############################################################################
        if IsPointAttractor == True:
            if str(eval('{}.mRNA'.format(GRN_List[0]))) in StableStatesCollector:
                StableStatesCollector[str(eval('{}.mRNA'.format(GRN_List[0])))] = StableStatesCollector[str(eval(
                    '{}.mRNA'.format(GRN_List[0])))] + 1
            else:
                StableStatesCollector[str(eval('{}.mRNA'.format(GRN_List[0])))] = 1
        elif IsPointAttractor == False:
            pass
        else:
            raise Exception('Error 007!')
        ########################################################################################################################################################################################

        ####################################################################################### Make plots #####################################################################################
        Plot_Curve = [[] for i in range(0, TotalNumberOfGenes)]
        Average_level = []

        for i in range(0, TotalNumberOfGenes):
            for j in range(0, Total_Time_Span):
                Plot_Curve[i].append(npmRNA_continuous[j][i])

        for i in range(0, TotalNumberOfGenes):
            Average_level.append(np.mean(Plot_Curve[i][-50:]))

        x = range(0, Total_Time_Span)
        colormap = ['red', 'blue', 'green', 'orange', 'black', 'yellow', 'purple', 'gold', 'silver']
        genelabels = ['Gene{}'.format(Genei) for Genei in range(0, TotalNumberOfGenes)]
        '''
        for i in range(0, TotalNumberOfGenes):
            plt.plot([50 * repeat for repeat in range(0, 6)], [WTTP[str(Global_i)][1][i] for repeat in range(0, 6)],
                     color=colormap[i], linestyle=':', marker='o')
        '''
        for i in range(0, TotalNumberOfGenes):
            if Number_of_initial_positions == 0:
                #pass
                plt.plot(x, Plot_Curve[i], c=colormap[i], label=genelabels[i], linewidth = 0.5)
            else:
                #pass
                plt.plot(x, Plot_Curve[i], c=colormap[i], linewidth = 0.5)
        ########################################################################################################################################################################################

        Loopcounter = Loopcounter + 1
            # plt.axhline(y=sys_WTTP[str(Global_i)][1][i], color=colormap[i], linestyle=':', marker = 'o')

            # plt.plot([ xasdf for xasdf in range(0, len(sol.y[5+i]))],sol.y[5+i] )
            # plt.plot(x, Plot_Curve[i])
        # plt.ylim((0,850))
    # plt.ylim((0,1.1*max(sys_WTTP[str(Global_i)][1])))
    plt.legend()
    plt.xlabel('Iteration', fontsize=25)
    plt.ylabel('Expression level', fontsize=25)
    plt.title('Gene expression simulation, N={}'.format(Number_of_initial_positions+1), fontsize=25)
    plt.show()
    #print('Experimental Attractor:', list(np.round(WTTP[str(Global_i)][1], 2)))
    # print('Last state:', mRNACheckList[-1])
    #print('Model Attractor:       ', list(np.round(np.array((Average_level)), 2)))

    # print('MSE\t', mean_squared_error(WTTP[str(Global_i)][1], np.array((Average_level))), '\t', mean_squared_error(WTTP[str(Global_i)][1], sys_WTTP_C[str(Global_i)][1]))
    if GetAttractorDistance(WTTP[str(Global_i)][1], np.array((Average_level)),
                            TranscriptionPofileMax) / TotalNumberOfGenes < 0.2:
        Matched_Number = Matched_Number + 1
    else:
        pass
    Calculate_AttractorDistance.append(np.mean(npmRNA_continuous[-50:], axis=0))

CurrentDistance = []
for G1 in range(0, len(WTTP)):
    CurrentDistance.append(
        GetAttractorDistance(Calculate_AttractorDistance[G1], WTTP[str(G1)][1], TranscriptionPofileMax))
# print(np.sum(CurrentDistance))

Attractor_Distance_matrix = np.zeros((len(WTTP), len(WTTP)))
for novel_i in range(0, len(Calculate_AttractorDistance)):
    for novel_j in range(0, len(Calculate_AttractorDistance)):
        Attractor_Distance_matrix[novel_i][novel_j] = GetAttractorDistance(Calculate_AttractorDistance[novel_i],
                                                                           Calculate_AttractorDistance[novel_j],
                                                                           TranscriptionPofileMax)
Unique_attractor_N = (1 + np.count_nonzero(np.mean(Attractor_Distance_matrix, axis=0) > 1))

# print(np.sum(Calculate_AttractorDistance, axis=0))
#print(np.sum(CurrentDistance))
#print('# Unique Attractors: ', Unique_attractor_N)

CurrentDistance = np.array(CurrentDistance)
#print('Mean Attractor Distance: ', np.mean(CurrentDistance) / Unique_attractor_N)
# print(len(WTTP), '\t', Matched_Number)
# print('\n', sorted(StableStatesCollector.items(), key=lambda StableStatesCollector:StableStatesCollector[1]))
#print(GRN_1.TranscriptionRate)