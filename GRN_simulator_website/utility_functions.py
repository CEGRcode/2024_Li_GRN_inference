import numpy as np
import copy
import json
import os
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
from sklearn import preprocessing
import colorsys
import math
import h5py
import itertools
import random

def n_smallest_indices(lst, n):
    if n > len(lst):
        raise ValueError("n cannot be greater than the number of elements in the list.")
    indexed_list = list(enumerate(lst))
    sorted_list = sorted(indexed_list, key=lambda x: x[1])
    indices = [sorted_list[i][0] for i in range(n)]
    return indices

def Sort_LG(_Matrix_LG):
    Out_Matrix_LG = []
    for j in range(0, len(_Matrix_LG)):
        _Matrix = _Matrix_LG[j]
        Complexes = {}
        for i in range(0, len(_Matrix)):
            if _Matrix[i] not in Complexes:
                Complexes[_Matrix[i]] = [i]
            else:
                Complexes[_Matrix[i]].append(i)
        Temp_LG_List = ['a' for i in range(0, len(_Matrix))]
        To_Count = 0
        To_fill = 0
        Done_complexes = {}
        while len(Done_complexes) < len(Complexes):
            for keys in Complexes:
                if keys not in Done_complexes:
                    if To_Count in Complexes[keys]:
                        for each_index in Complexes[keys]:
                            Temp_LG_List[each_index] = To_fill
                        To_fill = To_fill + 1
                        Done_complexes[keys] = ''
                    else:
                        pass
                else:
                    pass
            To_Count = To_Count + 1
        Out_Matrix_LG.append(Temp_LG_List)
    return Out_Matrix_LG

def Matrix2String01_LG(_Matrix):
    '''Convert numpy array into 1-D string representation of Logic Gate Matrix'''
    _Matrix = np.array(_Matrix)
    OutString01 = ''
    for i in range(0, _Matrix.shape[0]):
        for j in range(0, _Matrix.shape[1]):
            OutString01 = OutString01 + str(_Matrix[i][j])
    return OutString01

def Matrix2String01_LG_Expanded(_Matrix):
    '''Convert numpy array into 1-D string representation of Logic Gate Matrix'''
    _Matrix = np.array(_Matrix)
    OutString01 = ''
    for i in range(0, _Matrix.shape[0]):
        for j in range(0, _Matrix.shape[1]):
            OutString01 = OutString01 + str(_Matrix[i][j]) + ','
    OutString01 = OutString01[:-1]
    return OutString01

def LogicGatesString2Matrix(string):
    '''Convert 1-D string representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    outmatrix = np.random.randint(3, 4, (int(len(string)/2), 2))
    stringindex = 0
    for i in range(0, outmatrix.shape[0]):
        for j in range(0, outmatrix.shape[1]):
            outmatrix[i][j] = int(string[stringindex])
            stringindex = stringindex + 1
    return outmatrix

def LogicGatesString2Matrix_Expanded(string):
    '''Convert 1-D string representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    string_matrix = string.split(',')
    SQRTLEN = int(np.sqrt(len(string_matrix)))
    outmatrix = np.random.randint(8, 9, (SQRTLEN, SQRTLEN))
    stringindex = 0
    for i in range(0, outmatrix.shape[0]):
        for j in range(0, outmatrix.shape[1]):
            #print(stringindex)
            outmatrix[i][j] = int(string_matrix[stringindex])
            stringindex = stringindex + 1
    return outmatrix

def GetCorespondingMatrix(List012):
    '''Convert 1-D array representation of Weighted Adjacency Matrix into numpy array for fast matrix manipulations'''
    SQRTLEN = int(np.sqrt(len(List012)))
    OutMatrix = np.random.randint(8, 9, (2, SQRTLEN, SQRTLEN))
    for i in range(0, len(List012)):
        if List012[i] == 0:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == 1:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 1
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == 2:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 1
        else:
            raise Exception('List012 must just have 012!\n{}'.format(List012))
    return OutMatrix


def String012ToMatrix(List012):
    '''Convert 1-D string representation of Weighted Adjacency Matrix into numpy array for fast matrix manipulations'''
    # String '012' to Matrix
    SQRTLEN = int(np.sqrt(len(List012)))
    OutMatrix = np.random.randint(8, 9, (2, SQRTLEN, SQRTLEN))
    for i in range(0, len(List012)):
        if List012[i] == '0':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == '1':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 1
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == '2':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 1
        else:
            raise Exception('List012 must just have 012!\n{}'.format(List012))
    return OutMatrix


def GetCorespondingMatrix_LG(List01):
    '''Convert 1-D array representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    SQRTLEN = int(len(List01)/2)
    OutMatrix = np.random.randint(8, 9, (SQRTLEN, 2))
    for i in range(0, len(List01)):
        OutMatrix[i//2][i % 2] = int(List01[i])
    return OutMatrix


def ConfigurationTo012(ConfigurationMatrix):
    '''Convert configuration matrix to string'''
    OutString = ''
    for j in range(0, ConfigurationMatrix.shape[1]):
        for z in range(0, ConfigurationMatrix.shape[2]):
            if ConfigurationMatrix[0][j][z] == 0 and ConfigurationMatrix[1][j][z] == 0:
                OutString = OutString + '0'
            elif ConfigurationMatrix[0][j][z] == 1 and ConfigurationMatrix[1][j][z] == 0:
                OutString = OutString + '1'
            elif ConfigurationMatrix[0][j][z] == 0 and ConfigurationMatrix[1][j][z] == 1:
                OutString = OutString + '2'
            else:
                raise Exception('Configuration (1,1)!')
    return OutString


def maskmean(array_):
    '''mean of array_ with entries called 'mask' removed'''
    out_mean = []
    for each_array in array_.T:
        temp_mean = []
        for each in each_array[0]:
            if each == 'mask':
                pass
            else:
                temp_mean.append(float(each))
        out_mean.append(np.nanmean(temp_mean))
    return out_mean

def GetmRNASearchSpace(TotalGeneNum, RunNum, BoundList):
    '''obtain the coordinates of evenly-distributed vertices in a n-dimensional space'''
    BoundDeltaList = []
    for i in range(0, len(BoundList)):
        BoundDeltaList.append(abs(BoundList[i][0] - BoundList[i][1]))
    if len(BoundDeltaList) == TotalGeneNum:
        pass
    else:
        raise Exception('len(BoundDeltaList) != TotalGeneNum')

    NumLevels = 0
    for i in range(2, 999):
        if math.log(RunNum, i) >= TotalGeneNum:
            NumLevels = i
        else:
            pass

    if NumLevels > 1:
        if np.floor(math.log(RunNum / ((NumLevels - 1) ** TotalGeneNum), 10) / math.log(NumLevels / (NumLevels - 1), 10)) > TotalGeneNum:
            NumLevels = NumLevels + 1
            NumHigh = int(np.floor(
                math.log(RunNum / ((NumLevels - 1) ** TotalGeneNum), 10) / math.log(NumLevels / (NumLevels - 1), 10)))
            NumLow = int(TotalGeneNum - NumHigh)
        else:
            NumLow = 0
    else:
        NumLow = 0

    Delta = []
    CounterListEndState = []
    OutList = []
    CounterList = []
    if NumLevels >= 2:
        IndexOfLowLevelGenes = n_smallest_indices(BoundDeltaList, NumLow)
        for i in range(0, len(BoundDeltaList)):
            if i in IndexOfLowLevelGenes:
                Delta.append(BoundDeltaList[i] / (NumLevels - 1))
                CounterListEndState.append(BoundDeltaList[i] * ((NumLevels - 2) / (NumLevels - 1)))
                CounterList.append(0)
            else:
                Delta.append(BoundDeltaList[i] / (NumLevels))
                CounterList.append(0)
                CounterListEndState.append(BoundDeltaList[i] * ((NumLevels - 1) / NumLevels))
    else:
        NumOf1LevelGenes = int(np.round(TotalGeneNum - math.log(RunNum, 2)))
        IndexOf1LevelGenes = n_smallest_indices(BoundDeltaList, NumOf1LevelGenes)
        for i in range(0, len(BoundDeltaList)):
            if i in IndexOf1LevelGenes:
                Delta.append(0)
                CounterListEndState.append(BoundDeltaList[i] * 0.5)
                CounterList.append(BoundDeltaList[i] * 0.5)
            else:
                Delta.append(BoundDeltaList[i])
                CounterListEndState.append(BoundDeltaList[i])
                CounterList.append(0)
    CounterListInitialState = copy.deepcopy(CounterList)

    while True:
        OutList.append(copy.deepcopy(CounterList))
        if CounterList == CounterListEndState:
            break
        else:
            pass
        for i in range(0, len(CounterList)):
            if CounterList[i] >= CounterListEndState[i]:
                CounterList[i] = CounterListInitialState[i]
            else:
                CounterList[i] = CounterList[i] + Delta[i]
                break

    for i in range(0, len(OutList)):
        for j in range(0, len(OutList[i])):
            OutList[i][j] = np.round(OutList[i][j], 2)

    return np.array(OutList)

def json2ea(json_data):
    '''Convert .json to input format for EA'''
    # Parse JSON data
    data = json.loads(json_data)

    # Extract nodes and edges from JSON
    nodes = {node["id"]: {k: v for k, v in node.items() if k != "id"} for node in data["nodes"]}
    edges = [(edge["source"], edge["target"], edge["label"], edge["style"]) for edge in data["edges"]]

    # Identify unique node IDs
    Order_of_genes = list(nodes.keys())
    num_nodes = len(Order_of_genes)

    # Initialize adjacency matrix
    adjacency_matrix = [[0] * num_nodes for _ in range(num_nodes)]

    # Get the adjacency matrix
    for edge in edges:
        source_index = Order_of_genes.index(edge[0])
        target_index = Order_of_genes.index(edge[1])
        if edge[3][0] == 'triangle':
            adjacency_matrix[source_index][target_index] = 1
        elif edge[3][0] == 'tee':
            adjacency_matrix[source_index][target_index] = 2
        else:
            adjacency_matrix[source_index][target_index] = 0
    AM = ''
    for each_x in adjacency_matrix:
        for each_y in each_x:
            AM = AM + str(each_y)

    # Get the f0
    f0 = []
    f0p = []
    k = []
    t1 = []
    t2 = []
    t3 = []
    c1 = []
    c2 = []
    c3 = []
    c4 = []
    TR = []
    Deg = []
    Lk = []
    for node in Order_of_genes:
        f0.append(nodes[node]['f0'])
        f0p.append(nodes[node]['f0p'])
        k.append(nodes[node]['k'])
        t1.append(nodes[node]['t1'])
        t2.append(nodes[node]['t2'])
        t3.append(nodes[node]['t3'])
        c1.append(nodes[node]['c1'])
        c2.append(nodes[node]['c2'])
        c3.append(nodes[node]['c3'])
        c4.append(nodes[node]['c4'])
        TR.append(nodes[node]['TR'])
        Deg.append(nodes[node]['Deg'])
        Lk.append(nodes[node]['Lk'])

    f0 = str(f0)
    f0p = str(f0p)
    k = str(k)
    t1 = str(t1)
    t2 = str(t2)
    t3 = str(t3)
    c1 = str(c1)
    c2 = str(c2)
    c3 = str(c3)
    c4 = str(c4)
    TR = str(TR)
    Deg = str(Deg)
    Lk = str(Lk)

    temp_LG = {}

    for edge in edges:
        for each in Order_of_genes:
            if edge[1] == each and edge[3][1] == 'dashed':
                if edge[1] not in temp_LG:
                    temp_LG[edge[1]] = {edge[0]: edge[2]}
                else:
                    temp_LG[edge[1]].update({edge[0]: edge[2]})

    LG = []
    for each in Order_of_genes:
        if each not in temp_LG:
            LG.append([i for i in range(0, num_nodes)])
        else:
            LG.append(['a' for i in range(0, num_nodes)])
            counter = 0
            copy_temp_LG = copy.deepcopy(temp_LG[each])
            for gene in Order_of_genes:
                if gene in copy_temp_LG:
                    keys_to_remove = []
                    for keys in temp_LG[each]:
                        if (gene in temp_LG[each]) and (temp_LG[each][keys] == temp_LG[each][gene]):
                            LG[-1][Order_of_genes.index(keys)] = counter
                            keys_to_remove.append(keys)
                        else:
                            pass
                    for key_to_remove in keys_to_remove:
                        del temp_LG[each][key_to_remove]
                    if keys_to_remove == []:
                        pass
                    else:
                        counter = counter + 1
                else:
                    LG[-1][Order_of_genes.index(gene)] = counter
                    counter = counter + 1

    LG_ = ''
    for each_i in LG:
        for each_j in each_i:
            LG_ = LG_ + str(each_j) + ','
    LG_ = LG_[:-1]

    return AM, f0, f0p, k, t1, t2, t3, c1, c2, c3, c4, TR, Deg, Lk, LG_, '\t'.join(Order_of_genes)

def ea2json(AM, LG, f0, Order_of_genes, Vmax=''):
    '''Convert EA output to .json'''
    if Vmax == '':
        Vmax = [1 for i in range(0, len(Order_of_genes))]
    else:
        Vmax = eval(Vmax)

    json = ''
    # Add the nodes
    f0 = eval(f0)

    LG = LG.split(',')

    json = json + '{\n  "nodes": [\n'
    for i in range(0, len(Order_of_genes)):
        json = json + '\t{\n\t  ' + '"id": "{}",\n'.format(Order_of_genes[i]) + '\t  "label": "{}",\n'.format(Order_of_genes[i]) + '\t  "sua7Occupancy": {},\n'.format(Vmax[i]) + '\t  "f0": {}'.format(f0[i])+ '\n  \t},\n'
    json = json[:-2]

    # Add the edges
    json = json + '\n  ],\n  "edges": [\n'
    num_nodes = len(f0)
    for i in range(0, len(AM)):
        if AM[i] == '0':
            continue
        elif AM[i] == '1':
            LG_slide = LG[num_nodes*(i%num_nodes):num_nodes*(1+i%num_nodes)]
            #print(i, LG_slide, i//num_nodes)
            if LG_slide.count(LG_slide[i//num_nodes]) > 1:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle"\n\t  ]'+ '\n  \t},\n'
            else:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"triangle"\n\t  ]'+ '\n  \t},\n'
        else:
            LG_slide = LG[num_nodes*(i%num_nodes):num_nodes*(1+i%num_nodes)]
            #print(i, LG_slide, i//num_nodes)
            if LG_slide.count(LG_slide[i//num_nodes]) > 1:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"tee"\n\t  ]'+ '\n  \t},\n'
            else:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"tee"\n\t  ]'+ '\n  \t},\n'
    json = json[:-2]
    json = json + '\n  ]\n}'
    return json

# Generate colors for the composite plot
def generate_distinct_colors(n):
    colors = []
    for i in range(n):
        hue = i / n
        lightness = (50 + 20 * (i % 2)) / 100
        saturation = 90 / 100
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
        colors.append(rgb)
    return colors

# Read in data for the composite plot
def read_cdt(forward_file, reverse_file):
    WindowSize = 500 # compensate the removed double tails in return. Adjust the x-axis accordingly.
    length = 500
    forward_CDT = [0 for i in range(0, length)]
    reverse_CDT = [0 for i in range(0, length)]

    Skip_first_row = True
    forward_input = open(forward_file, 'r')
    for line in forward_input:
        if Skip_first_row:
            Skip_first_row = False
            continue
        else:
            pass
        for i in range(0, len(line.split()[2:])):
            forward_CDT[i] = forward_CDT[i] + float(line.split()[i+2])
    forward_input.close()

    Skip_first_row = True
    reverse_input = open(reverse_file, 'r')
    for line in reverse_input:
        if Skip_first_row:
            Skip_first_row = False
            continue
        else:
            pass
        for i in range(0, len(line.split()[2:])):
            reverse_CDT[i] = reverse_CDT[i] + float(line.split()[i+2])
    reverse_input.close()
    
    forward_CDT = preprocessing.normalize([np.array(forward_CDT[1:-1])])[0]
    reverse_CDT = preprocessing.normalize([np.array(reverse_CDT[1:-1])])[0]
    
    forward_CDT = forward_CDT[(int(length/2)-int(WindowSize/2)):(int(length/2)+int(WindowSize/2))]
    reverse_CDT = reverse_CDT[(int(length/2)-int(WindowSize/2)):(int(length/2)+int(WindowSize/2))]
    
    return forward_CDT, reverse_CDT

# Read in data in a h5 file for the composite plot
def read_cdt_h5(BED_ID, sample_ID):
    WindowSize = 500 # compensate the removed double tails in return. Adjust the x-axis accordingly.
    length = 500
    forward_CDT = [0 for i in range(0, length)]
    reverse_CDT = [0 for i in range(0, length)]

    if sample_ID == 'IgG':
        with h5py.File('./static/data/IgG_BED_combined.h5', 'r') as h5file:
            for each in h5file['{}_BED_{}_BAM_sense'.format(BED_ID, sample_ID)][:]:
                if binary_to_string(each).startswith('YORF'):
                    continue
                else:
                    for i in range(0, len(each.split()[2:])):
                        forward_CDT[i] = forward_CDT[i] + float(each.split()[i+2])
            for each in h5file['{}_BED_{}_BAM_anti'.format(BED_ID, sample_ID)][:]:
                if binary_to_string(each).startswith('YORF'):
                    continue
                else:
                    for i in range(0, len(each.split()[2:])):
                        reverse_CDT[i] = reverse_CDT[i] + float(each.split()[i+2])  
    else:
        with h5py.File('./static/data/{}_BED_combined.h5'.format(BED_ID), 'r') as h5file:
            for each in h5file['{}_BED_{}_BAM_sense_SCALE'.format(BED_ID, sample_ID)][:]:
                if binary_to_string(each).startswith('YORF'):
                    continue
                else:
                    for i in range(0, len(each.split()[2:])):
                        forward_CDT[i] = forward_CDT[i] + float(each.split()[i+2])
            for each in h5file['{}_BED_{}_BAM_anti_SCALE'.format(BED_ID, sample_ID)][:]:
                if binary_to_string(each).startswith('YORF'):
                    continue
                else:
                    for i in range(0, len(each.split()[2:])):
                        reverse_CDT[i] = reverse_CDT[i] + float(each.split()[i+2])
    
    return forward_CDT, reverse_CDT


def binary_to_string(binary_data):
    # Check if binary_data is a numpy.bytes_ object
    if isinstance(binary_data, np.bytes_):
        # Convert numpy.bytes_ to bytes
        binary_data = bytes(binary_data)

    # Convert binary data to a normal string
    # Convert bytes to a binary string
    binary_string = ''.join(format(byte, '08b') for byte in binary_data)

    # Split binary string into chunks of 8 bits (1 byte)
    byte_size = 8
    binary_values = [binary_string[i:i + byte_size] for i in range(0, len(binary_string), byte_size)]

    # Convert each byte to an integer and then to a character
    characters = [chr(int(bv, 2)) for bv in binary_values]

    # Join all characters to form the final string
    return ''.join(characters)

def GetRegulatorForGene(List012, i):
    if np.sqrt(len(List012)) < i:
        raise Exception('Gene index out of range.')
    else:
        outstring = ''
        for j in range(0, len(List012)):
            if j%np.sqrt(len(List012)) == i:
                outstring = outstring + List012[j]
            else:
                pass
    return outstring

def n_smallest_indices(lst, n):
    if n > len(lst):
        raise ValueError("n cannot be greater than the number of elements in the list.")
    indexed_list = list(enumerate(lst))
    sorted_list = sorted(indexed_list, key=lambda x: x[1])
    indices = [sorted_list[i][0] for i in range(n)]
    return indices

def split_and_average(lst):
    # Split list into 3 groups as evenly as possible
    if len(lst) == 1:
        return [lst[0], lst[0]]
    elif len(lst) == 2:
        return [lst[0], lst[1]]
    else:
        pass

    split_lists = np.array_split(lst, 3)
    group1, group2, group3 = [list(group) for group in split_lists]

    # Compute required averages
    avg1 = (max(group1) + min(group2)) / 2
    avg2 = (max(group2) + min(group3)) / 2

    return [avg1, avg2]

def allocate_samples(total_iterations, proportions):
    """
    Determines the number of samples per range given total_iterations and proportions.
    Uses binary search to maximize the counts while keeping product <= total_iterations.
    """
    # Normalize proportions in case they don't sum to 1.
    proportions = np.array(proportions) / np.sum(proportions)
    print(proportions)
    n = len(proportions)
    # Binary search for scaling factor k.
    lo, hi = 1, total_iterations  # k must be at least 1.
    best_counts = [1] * n  # fallback: at least one sample per range.
    
    while lo <= hi:
        mid = (lo + hi) // 2
        # Calculate counts: ensure at least 1 sample per range.
        counts = [max(round(p * mid), 1) for p in proportions]
        prod = math.prod(counts)
        if prod <= total_iterations:
            best_counts = counts  # valid, try for a larger k.
            lo = mid + 1
        else:
            hi = mid - 1
            
    index_non_zero_proportions = [i for i, val in enumerate(proportions) if val != 0]
    last_index_to_add = 0
    while math.prod(best_counts) < 0.9*total_iterations:
        min_value = min(best_counts)
        min_index = [i for i, val in enumerate(best_counts) if val == min_value]
        random_min_index_non_zero_in_proportions = random.choice(list(set(index_non_zero_proportions) & set(min_index)))
        best_counts[random_min_index_non_zero_in_proportions] = best_counts[random_min_index_non_zero_in_proportions] + 1
        last_index_to_add = random_min_index_non_zero_in_proportions

    if math.prod(best_counts) > total_iterations:
        best_counts[last_index_to_add] = best_counts[last_index_to_add] - 1
    return best_counts
    
def generate_limited_combinations(total_iterations, range_limits, proportions):
    """
    Generates a limited number of combinations by sampling each range non-uniformly.
    """
    # Normalize proportions in case they don't sum to 1.
    proportions = np.array(proportions) / np.sum(proportions)
    print(proportions)
    # Determine number of samples for each range.
    sample_counts = allocate_samples(total_iterations, proportions)
    
    # Generate the sampled values for each range.
    sampled_ranges = [
        np.linspace(start, end, num, dtype=int) if num > 1 else np.array([start])
        for (start, end), num in zip(range_limits, sample_counts)
    ]
    
    # Generate Cartesian product of the sampled ranges.
    combinations = list(itertools.product(*sampled_ranges))
    return combinations