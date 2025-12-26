import copy
import numpy as np
import math
import random
import time
from distance_functions import *
from mutation_functions import *
from dynamics import *
import sympy as sp
# non-overlapping generations in a haploid infinite population, the change in relative frequency is proportional to the product of fitness and frequency of a given type divided by the mean fitness.


class GRN:

    def __init__(self, Name, mRNA, Configuration,
                MutationRate, TranscriptionRate,
                  DegradationRatemRNA,
                 TranscriptionThreshold, LogicGates, Sigmoid_k,
                 Leakage, f0_c, sys_input_ChIP='', mutation_phase=1):
        '''Parameters Initiation'''
        self.Name = Name
        self.mRNA = mRNA
        self.mutation_phase = mutation_phase
        # Configuration.shape = (2,m,n); (0,m,n) for Activation and (1,m,n) for Repression
        self.Configuration = Configuration
        self.MutationRate = MutationRate
        self.TranscriptionRate = TranscriptionRate
        self.DegradationRatemRNA = DegradationRatemRNA
        # For each GRN, its shape should be (7000,7000), Genei is an activator of Genej.
        self.TranscriptionThreshold = TranscriptionThreshold
        self.LogicGates = LogicGates  # For each gene, it has 2 parameters
        self.Sigmoid_k = Sigmoid_k
        self.f0_c = f0_c
        self.Leakage = Leakage
        self.Memo_mRNA = []
        self.sol = []
        self.sys_input_ChIP = sys_input_ChIP

    def ResetMemo_mRNA(self):
        self.Memo_mRNA = []
        self.sol = []
        return

    def UpdateProportion(self, NewProportion):
        '''Replace the old proportion by the new one'''
        self.Proportion = NewProportion
        if self.Proportion < 0 or self.Proportion > 1:
            raise Exception('Proportion Error!')
        return

    def GetDistance(self, WT_mRNA, transcriptionalprofilemax):
        '''Reads are equal'''
        # return(np.linalg.norm(self.mRNA - WT_mRNA))
        '''Genes are equal'''
        outdistance = 0.0
        for y in range(0, len(WT_mRNA)):
            outdistance = outdistance + \
                abs((self.mRNA[y]-WT_mRNA[y])/transcriptionalprofilemax[y])
        return outdistance

    def GetActivatorsRepressors(self):
        '''Return the indices of activators[0]/Repressors[1] for all genes'''
        arlist = [[i for i, c in enumerate(self.Configuration) if c == '1'], [i for i, c in enumerate(self.Configuration) if c == '2']]
        return arlist

    def Delta_mRNA(self, knockoutlist_in, TRM, initial_values, analytical_expression=0):
        '''Calculating the Delta mRNA'''
        specific_gene = int(self.Name.split('_')[1])
        ARList = self.GetActivatorsRepressors()
        scipystring = 'def update_mRNA_protein(t, y, knockoutlist):'

        # Current Activators are ARList[0]; Current Repressors are ARList[1]

        '''Form the differential equations'''
        TempString_CA1_sci = ''  # Scipy
        TempString_CA2_sci = ''  # Scipy
        TempString_CR1_sci = ''  # Scipy
        TempString_CR2_sci = ''  # Scipy

        '''Determine CA_1'''
        CA_complexes = {}
        for j in range(0, len(ARList[0])):
            if self.LogicGates[ARList[0][j]] in CA_complexes:
                CA_complexes[self.LogicGates[ARList[0][j]]].append(ARList[0][j])
            else:
                CA_complexes[self.LogicGates[ARList[0][j]]] = [ARList[0][j]]

        TempString_CSA_List = []
        for each_complex in CA_complexes:
            TempString_CA_complexes = ''
            for each_subunit in CA_complexes[each_complex]:
                if analytical_expression:
                    TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y[{}]'.format(each_subunit),
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][1],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][2]) + TempString_CA_complexes  # Scipy
                elif each_subunit == specific_gene:
                    TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y',
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][1],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][2]) + TempString_CA_complexes  # Scipy
                else:
                    TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        initial_values[each_subunit],
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][1],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][2]) + TempString_CA_complexes  # Scipy
            TempString_CA_complexes = TempString_CA_complexes[:-1]
            TempString_CSA_List.append(TempString_CA_complexes)

        for each in TempString_CSA_List:
            TempString_CA1_sci = TempString_CA1_sci + '(1-{})*'.format(each)

        TempString_CA1_sci = TempString_CA1_sci[:-1]

        if TempString_CA1_sci == '':
            TempString_CA1_sci = '0'
        else:
            TempString_CA1_sci = '(1-' + TempString_CA1_sci + ')'

        '''Determine CA_2'''
        CA_complexes = {}
        for j in range(0, len(ARList[0])):
            if self.LogicGates[ARList[0][j]] in CA_complexes:
                CA_complexes[self.LogicGates[ARList[0][j]]].append(ARList[0][j])
            else:
                CA_complexes[self.LogicGates[ARList[0][j]]] = [ARList[0][j]]

        TempString_CSA_List = []
        for each_complex in CA_complexes:
            TempString_CA_complexes = ''
            for each_subunit in CA_complexes[each_complex]:
                if analytical_expression:
                    TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y[{}]'.format(each_subunit),
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][1],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][3]) + TempString_CA_complexes  # Scipy
                elif each_subunit == specific_gene:
                    TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y',
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][1],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][3]) + TempString_CA_complexes  # Scipy
                else:
                    TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        initial_values[each_subunit],
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][1],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][3]) + TempString_CA_complexes  # Scipy
            TempString_CA_complexes = TempString_CA_complexes[:-1]
            TempString_CSA_List.append(TempString_CA_complexes)

        for each in TempString_CSA_List:
            TempString_CA2_sci = TempString_CA2_sci + '(1-{})*'.format(each)

        TempString_CA2_sci = TempString_CA2_sci[:-1]

        if TempString_CA2_sci == '':
            TempString_CA2_sci = '0'
        else:
            TempString_CA2_sci = '(1-' + TempString_CA2_sci + ')'

        '''Determine CR1'''
        CR_complexes = {}
        for j in range(0, len(ARList[1])):
            if self.LogicGates[ARList[1][j]] in CR_complexes:
                CR_complexes[self.LogicGates[ARList[1][j]]].append(ARList[1][j])
            else:
                CR_complexes[self.LogicGates[ARList[1][j]]] = [ARList[1][j]]

        TempString_CSR_List = []
        for each_complex in CR_complexes:
            TempString_CR_complexes = ''
            for each_subunit in CR_complexes[each_complex]:
                if analytical_expression:
                    TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y[{}]'.format(each_subunit),
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][2],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][4]) + TempString_CR_complexes  # Scipy
                elif each_subunit == specific_gene:
                    TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y',
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][2],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][4]) + TempString_CR_complexes  # Scipy
                else:
                    TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        initial_values[each_subunit],
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][2],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][4]) + TempString_CR_complexes  # Scipy
            TempString_CR_complexes = TempString_CR_complexes[:-1]
            TempString_CSR_List.append(TempString_CR_complexes)

        for each in TempString_CSR_List:
            TempString_CR1_sci = TempString_CR1_sci + '(1-{})*'.format(each)
        TempString_CR1_sci = TempString_CR1_sci[:-1]
        if TempString_CR1_sci == '':
            TempString_CR1_sci = '0'
        else:
            TempString_CR1_sci = '(1-' + TempString_CR1_sci + ')'

        '''Determine CR2'''
        CR_complexes = {}
        for j in range(0, len(ARList[1])):
            if self.LogicGates[ARList[1][j]] in CR_complexes:
                CR_complexes[self.LogicGates[ARList[1][j]]].append(ARList[1][j])
            else:
                CR_complexes[self.LogicGates[ARList[1][j]]] = [ARList[1][j]]

        TempString_CSR_List = []
        for each_complex in CR_complexes:
            TempString_CR_complexes = ''
            for each_subunit in CR_complexes[each_complex]:
                if analytical_expression:
                    TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y[{}]'.format(each_subunit),
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][2],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][5]) + TempString_CR_complexes  # Scipy
                elif each_subunit == specific_gene:
                    TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        'y',
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][2],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][5]) + TempString_CR_complexes  # Scipy
                else:
                    TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                        initial_values[each_subunit],
                        self.TranscriptionThreshold[specific_gene][each_subunit][0],
                        self.TranscriptionThreshold[specific_gene][each_subunit][2],
                        self.Sigmoid_k,
                        self.f0_c[specific_gene][5]) + TempString_CR_complexes  # Scipy
            TempString_CR_complexes = TempString_CR_complexes[:-1]
            TempString_CSR_List.append(TempString_CR_complexes)

        for each in TempString_CSR_List:
            TempString_CR2_sci = TempString_CR2_sci + '(1-{})*'.format(each)
        TempString_CR2_sci = TempString_CR2_sci[:-1]
        if TempString_CR2_sci == '':
            TempString_CR2_sci = '0'
        else:
            TempString_CR2_sci = '(1-' + TempString_CR2_sci + ')'

        if not analytical_expression:
            # f0 + (1-f0)*CA1 - CR1* f0 +(f0p-1+f0)*CA2 *CR2
            TempString_sci = ('{}+({}-{})*('.format(self.Leakage,
                                                    self.TranscriptionRate,
                                                    self.Leakage)
                            + '{}+(1-{})*{}-{}*{}+({}-1+{})*{}*{}'.format(
                self.f0_c[specific_gene][0],
                self.f0_c[specific_gene][0],
                TempString_CA1_sci,
                TempString_CR1_sci,
                self.f0_c[specific_gene][0],
                self.f0_c[specific_gene][1],
                self.f0_c[specific_gene][0],
                TempString_CA2_sci,
                TempString_CR2_sci)
                + ')-{}*y'.format(self.DegradationRatemRNA))

            scipystring = (scipystring
                        + '\n    delta_mRNA='
                        + TempString_sci)

            scipystring = scipystring + ('\n    if (y<=0 and delta_mRNA<=0) or ({} in knockoutlist):'
                                            '\n        outlist = -y'
                                            '\n    elif {} in knockoutlist:'
                                            '\n        outlist = {}-y'
                                            '\n    else:'
                                            '\n        outlist = delta_mRNA').format(specific_gene, -2-specific_gene, TRM[specific_gene])

            scipystring = scipystring + '\n    return outlist'
        else:
            # f0 + (1-f0)*CA1 - CR1* f0 +(f0p-1+f0)*CA2 *CR2
            TempString_sci = ('{}+({}-{})*('.format(self.Leakage,
                                                    self.TranscriptionRate,
                                                    self.Leakage)
                            + '{}+(1-{})*{}-{}*{}+({}-1+{})*{}*{}'.format(
                self.f0_c[specific_gene][0],
                self.f0_c[specific_gene][0],
                TempString_CA1_sci,
                TempString_CR1_sci,
                self.f0_c[specific_gene][0],
                self.f0_c[specific_gene][1],
                self.f0_c[specific_gene][0],
                TempString_CA2_sci,
                TempString_CR2_sci)
                + ')-{}*y[{}]'.format(self.DegradationRatemRNA, specific_gene))

            scipystring = (scipystring
                        + '\n    delta_mRNA='
                        + TempString_sci)

            scipystring = scipystring + ('\n    if (y<=0 and delta_mRNA<=0) or ({} in knockoutlist):'
                                            '\n        outlist = -y'
                                            '\n    elif {} in knockoutlist:'
                                            '\n        outlist = {}-y'
                                            '\n    else:'
                                            '\n        outlist = delta_mRNA').format(specific_gene, -2-specific_gene, TRM[specific_gene])

            scipystring = scipystring + '\n    return outlist'

        return(scipystring)

    def GenerateMutation(self, tested_AMLG, index_of_diff_gene, sys_LG_, only_TF_DNA=False, max_iter=1000):
        '''Mutations occur upon 1.Configuration; 2.LogicGates'''
        specific_gene = int(self.Name.split('_')[1])
        cache_index = index_of_diff_gene.index(specific_gene)
        self_LG_string = ",".join(map(str, self.LogicGates))
        tested_LG_strings = list(tested_AMLG[cache_index].keys())
        if self_LG_string in tested_LG_strings:
            tesetd_AM_strings = list(tested_AMLG[cache_index][self_LG_string])
        else:
            tesetd_AM_strings = []

        if type(self.sys_input_ChIP) == str:
            mutated_Configuration = HammingMutation_3(len(self.Configuration), tesetd_AM_strings,  list(self.sys_input_ChIP[:, specific_gene]), max_iter)
            if mutated_Configuration == False:
                for _ in range(0, max_iter):
                    mutated_LogicGates = Mutation_LG_Simple(self.MutationRate, sys_LG_[specific_gene], index_of_diff_gene)
                    if ",".join(map(str, mutated_LogicGates)) not in tested_LG_strings:
                        self.LogicGates = mutated_LogicGates
                        self.Configuration = HammingMutation_3(len(self.Configuration), [],  list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                        return
                    else:
                        mutated_Configuration = HammingMutation_3(len(self.Configuration), tested_AMLG[cache_index][",".join(map(str, mutated_LogicGates))],  list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                        if mutated_Configuration != False:
                            self.LogicGates = mutated_LogicGates
                            self.Configuration = mutated_Configuration
                            return
                        else:
                            continue
                self.MutationRate = 'All tested!'
                return
            else:
                self.Configuration = mutated_Configuration
                return
        else:
            mutated_Configuration = HammingMutation_1(len(self.Configuration), tesetd_AM_strings, list(self.sys_input_ChIP[:, specific_gene]), max_iter)
            if mutated_Configuration == False and (not only_TF_DNA):
                mutated_Configuration = HammingMutation_2(len(self.Configuration), tesetd_AM_strings, list(self.sys_input_ChIP[:, specific_gene]),max_iter)
                if mutated_Configuration == False:
                    mutated_Configuration = HammingMutation_3(len(self.Configuration), tesetd_AM_strings, list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                    if mutated_Configuration == False:
                        for _ in range(0, max_iter):
                            mutated_LogicGates = Mutation_LG_Simple(self.MutationRate, sys_LG_[specific_gene], index_of_diff_gene)
                            if ",".join(map(str, mutated_LogicGates)) not in tested_LG_strings:
                                self.LogicGates = mutated_LogicGates
                                self.Configuration = HammingMutation_1(len(self.Configuration), [], list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                                return
                            else:
                                mutated_Configuration = HammingMutation_1(len(self.Configuration), tested_AMLG[cache_index][",".join(map(str, mutated_LogicGates))], list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                                if mutated_Configuration != False:
                                    self.LogicGates = mutated_LogicGates
                                    self.Configuration = mutated_Configuration
                                    return
                                else:
                                    mutated_Configuration = HammingMutation_2(len(self.Configuration), tested_AMLG[cache_index][",".join(map(str, mutated_LogicGates))], list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                                    if mutated_Configuration != False:
                                        self.LogicGates = mutated_LogicGates
                                        self.Configuration = mutated_Configuration
                                        return
                                    else:
                                        mutated_Configuration = HammingMutation_3(len(self.Configuration), tested_AMLG[cache_index][",".join(map(str, mutated_LogicGates))], list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                                        if mutated_Configuration != False:
                                            self.LogicGates = mutated_LogicGates
                                            self.Configuration = mutated_Configuration
                                            return
                                        else:
                                            continue
                        self.MutationRate = 'All tested!'
                        return
                    else:
                        self.Configuration = mutated_Configuration
                        return
                else:
                    self.Configuration = mutated_Configuration
                    return
            else:
                if mutated_Configuration == False:
                    for _ in range(0, max_iter):
                        mutated_LogicGates = Mutation_LG_Simple(self.MutationRate, sys_LG_[specific_gene], index_of_diff_gene)
                        if ",".join(map(str, mutated_LogicGates)) not in tested_LG_strings:
                            self.LogicGates = mutated_LogicGates
                            self.Configuration = HammingMutation_1(len(self.Configuration), [], list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                            return
                        else:
                            mutated_Configuration = HammingMutation_1(len(self.Configuration), tested_AMLG[cache_index][",".join(map(str, mutated_LogicGates))], list(self.sys_input_ChIP[:, specific_gene]), max_iter)
                            if mutated_Configuration != False:
                                self.LogicGates = mutated_LogicGates
                                self.Configuration = mutated_Configuration
                                return
                            else:
                                continue
                    self.MutationRate = 'All tested!'
                    return
                else:
                    self.Configuration = mutated_Configuration
                    return

    def SetmRNA(self, NewmRNAList, knockout_list, Overexpression):
        #NewmRNAList = NewmRNAList.tolist()
        self.mRNA = copy.deepcopy(NewmRNAList)
        for i in range(0, len(knockout_list)):
            if knockout_list[i] >= 0:
                self.mRNA[knockout_list[i]] = 0
            else:
                self.mRNA[-knockout_list[i]-2] = Overexpression[i]
        return

    def SetMemo(self, NewMemo):
        self.Memo_mRNA = NewMemo
        return

    def SetMutationRate(self, NewMutationRate):
        self.MutationRate = copy.deepcopy(NewMutationRate)
        return

    def Setf0(self, Newf0):
        my_gene = int(self.Name.split('_')[1])
        self.f0_c[my_gene] = copy.deepcopy(Newf0)
        return

    def SetAM(self, NewAM):
        self.Configuration = copy.deepcopy(NewAM)
        return

    def Update_Transcription_Rate(self, TranscriptionRate):
        self.TranscriptionRate = TranscriptionRate
        return

    def Update_Promoter_Strength(self, which_i, sys_promoter_strengths, sys_mRNA_elongation_rate, sys_gene_length):
        TranscriptionRate = sys_promoter_strengths[which_i]*60*sys_mRNA_elongation_rate/sys_gene_length.tolist()
        self.TranscriptionRate = TranscriptionRate
        return

    def GetActivatorsRepressors_Network(self):
        '''Return the indices of activators[0]/Repressors[1] for all genes'''
        arlist = [[] for i in range(0, self.Configuration.shape[0])]
        for i in range(0, self.Configuration.shape[0]):
            temp_list = []
            for j in range(0, len(self.f0_c)):
                temp_list.append(
                    np.where(self.Configuration[i][:, j] == 1)[0].tolist())
            arlist[i] = temp_list
        return arlist

    def Delta_mRNA_Network(self, knockoutlist_in, TRM, indexes_of_diff_gene, initial_values):
        '''Calculating the Delta mRNA'''
        ARList = self.GetActivatorsRepressors_Network()
        scipystring = 'def update_mRNA_protein(t, y, knockoutlist):' + \
            '\n    delta_mRNA = [0]*len(y)' + '\n    outlist = [0]*len(y)'  # Scipy
        for i_of_diff in range(0, len(indexes_of_diff_gene)):
            i = indexes_of_diff_gene[i_of_diff]
            '''loop the mRNA list (For each gene)'''
            # Current Activators are ARList[0][i]; Current Repressors are ARList[1][i]

            '''Form the differential equations'''
            TempString_CA1_sci = ''  # Scipy
            TempString_CA2_sci = ''  # Scipy
            TempString_CR1_sci = ''  # Scipy
            TempString_CR2_sci = ''  # Scipy

            '''Determine CA_1'''
            CA_complexes = {}
            for j in range(0, len(ARList[0][i])):
                if self.LogicGates[i][ARList[0][i][j]] in CA_complexes:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]].append(ARList[0][i][j])
                else:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]] = [ARList[0][i][j]]

            TempString_CSA_List = []
            for each_complex in CA_complexes:
                TempString_CA_complexes = ''
                for each_subunit in CA_complexes[each_complex]:
                    if each_subunit in indexes_of_diff_gene:
                        TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            'y[{}]'.format(indexes_of_diff_gene.index(each_subunit)),
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][1],
                            self.Sigmoid_k[i],
                            self.f0_c[i][2]) + TempString_CA_complexes  # Scipy
                    else:
                        TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            initial_values[each_subunit],
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][1],
                            self.Sigmoid_k[i],
                            self.f0_c[i][2]) + TempString_CA_complexes  # Scipy
                TempString_CA_complexes = TempString_CA_complexes[:-1]
                TempString_CSA_List.append(TempString_CA_complexes)

            for each in TempString_CSA_List:
                TempString_CA1_sci = TempString_CA1_sci + '(1-{})*'.format(each)

            TempString_CA1_sci = TempString_CA1_sci[:-1]

            if TempString_CA1_sci == '':
                TempString_CA1_sci = '0'
            else:
                TempString_CA1_sci = '(1-' + TempString_CA1_sci + ')'

            '''Determine CA_2'''
            CA_complexes = {}
            for j in range(0, len(ARList[0][i])):
                if self.LogicGates[i][ARList[0][i][j]] in CA_complexes:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]].append(ARList[0][i][j])
                else:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]] = [ARList[0][i][j]]

            TempString_CSA_List = []
            for each_complex in CA_complexes:
                TempString_CA_complexes = ''
                for each_subunit in CA_complexes[each_complex]:
                    if each_subunit in indexes_of_diff_gene:
                        TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            'y[{}]'.format(indexes_of_diff_gene.index(each_subunit)),
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][1],
                            self.Sigmoid_k[i],
                            self.f0_c[i][3]) + TempString_CA_complexes  # Scipy
                    else:
                        TempString_CA_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            initial_values[each_subunit],
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][1],
                            self.Sigmoid_k[i],
                            self.f0_c[i][3]) + TempString_CA_complexes  # Scipy
                TempString_CA_complexes = TempString_CA_complexes[:-1]
                TempString_CSA_List.append(TempString_CA_complexes)

            for each in TempString_CSA_List:
                TempString_CA2_sci = TempString_CA2_sci + '(1-{})*'.format(each)

            TempString_CA2_sci = TempString_CA2_sci[:-1]

            if TempString_CA2_sci == '':
                TempString_CA2_sci = '0'
            else:
                TempString_CA2_sci = '(1-' + TempString_CA2_sci + ')'

            '''Determine CR1'''
            CR_complexes = {}
            for j in range(0, len(ARList[1][i])):
                if self.LogicGates[i][ARList[1][i][j]] in CR_complexes:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]].append(ARList[1][i][j])
                else:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]] = [ARList[1][i][j]]

            TempString_CSR_List = []
            for each_complex in CR_complexes:
                TempString_CR_complexes = ''
                for each_subunit in CR_complexes[each_complex]:
                    if each_subunit in indexes_of_diff_gene:
                        TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            'y[{}]'.format(indexes_of_diff_gene.index(each_subunit)),
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][2],
                            self.Sigmoid_k[i],
                            self.f0_c[i][4]) + TempString_CR_complexes  # Scipy
                    else:
                        TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            initial_values[each_subunit],
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][2],
                            self.Sigmoid_k[i],
                            self.f0_c[i][4]) + TempString_CR_complexes  # Scipy
                TempString_CR_complexes = TempString_CR_complexes[:-1]
                TempString_CSR_List.append(TempString_CR_complexes)

            for each in TempString_CSR_List:
                TempString_CR1_sci = TempString_CR1_sci + '(1-{})*'.format(each)
            TempString_CR1_sci = TempString_CR1_sci[:-1]
            if TempString_CR1_sci == '':
                TempString_CR1_sci = '0'
            else:
                TempString_CR1_sci = '(1-' + TempString_CR1_sci + ')'

            '''Determine CR2'''
            CR_complexes = {}
            for j in range(0, len(ARList[1][i])):
                if self.LogicGates[i][ARList[1][i][j]] in CR_complexes:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]].append(ARList[1][i][j])
                else:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]] = [ARList[1][i][j]]

            TempString_CSR_List = []
            for each_complex in CR_complexes:
                TempString_CR_complexes = ''
                for each_subunit in CR_complexes[each_complex]:
                    if each_subunit in indexes_of_diff_gene:
                        TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            'y[{}]'.format(indexes_of_diff_gene.index(each_subunit)),
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][2],
                            self.Sigmoid_k[i],
                            self.f0_c[i][5]) + TempString_CR_complexes  # Scipy
                    else:
                        TempString_CR_complexes = 'modified_sigmoid({}, {}, {}, {}, {})*'.format(
                            initial_values[each_subunit],
                            self.TranscriptionThreshold[i][each_subunit][0],
                            self.TranscriptionThreshold[i][each_subunit][2],
                            self.Sigmoid_k[i],
                            self.f0_c[i][5]) + TempString_CR_complexes  # Scipy
                TempString_CR_complexes = TempString_CR_complexes[:-1]
                TempString_CSR_List.append(TempString_CR_complexes)

            for each in TempString_CSR_List:
                TempString_CR2_sci = TempString_CR2_sci + '(1-{})*'.format(each)
            TempString_CR2_sci = TempString_CR2_sci[:-1]
            if TempString_CR2_sci == '':
                TempString_CR2_sci = '0'
            else:
                TempString_CR2_sci = '(1-' + TempString_CR2_sci + ')'
            # f0 + (1-f0)*CA1 - CR1* f0 +(f0p-1+f0)*CA2 *CR2
            TempString_sci = ('{}+({}-{})*('.format(self.Leakage[i],
                                                    self.TranscriptionRate[i],
                                                    self.Leakage[i])
                            + '{}+(1-{})*{}-{}*{}+({}-1+{})*{}*{}'.format(
                self.f0_c[i][0],
                self.f0_c[i][0],
                TempString_CA1_sci,
                TempString_CR1_sci,
                self.f0_c[i][0],
                self.f0_c[i][1],
                self.f0_c[i][0],
                TempString_CA2_sci,
                TempString_CR2_sci)
                + ')-{}*y[{}]'.format(self.DegradationRatemRNA[i], i_of_diff))

            scipystring = (scipystring
                        + '\n    delta_mRNA[{}]='.format(i_of_diff)
                        + TempString_sci)

        for i in range(0, len(indexes_of_diff_gene)):
            scipystring = scipystring + ('\n    if (y[{}]<=0 and delta_mRNA[{}]<=0) or ({} in knockoutlist):'
                                         '\n        outlist[{}] = -y[{}]'
                                         '\n    elif {} in knockoutlist:'
                                         '\n        outlist[{}] = {}-y[{}]'
                                         '\n    else:'
                                         '\n        outlist[{}] = delta_mRNA[{}]').format(
                i, i, i, i, i,
                -2-i, i, TRM[i], i, i, i)

        scipystring = scipystring + '\n    return outlist'

        return(scipystring)