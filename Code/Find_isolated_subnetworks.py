import numpy as np
import pandas as pd

def Is_intersection(list_of_sets):
    '''Return true/false if some sets have intersection.'''
    is_inter = False
    for i in range(0, len(list_of_sets)):
        for j in range(0, len(list_of_sets)):
            if i == j:
                continue
            else:
                if len(list_of_sets[i].intersection(list_of_sets[j])) != 0:
                    is_inter = True
                else:
                    pass
    return is_inter

def Find_Isolated_Network(AM):
    '''Take the adjacency matrix and output the indice of the isolated networks.'''
    AM_subnetworks = []
    for i in range(0, AM.shape[0]):
        temp_subnetwork = []
        if i in temp_subnetwork or sum(abs(AM[i])) == 0:
            continue
        else:
            gene_to_search = i
            temp_subnetwork.append(gene_to_search)
        while not all([j[0] in temp_subnetwork for j in np.argwhere((AM[gene_to_search] != 0))]):
            '''Run untill all the non-zero entries are included in temp_subnetwork.'''
            gene_to_append = []
            #print('to_search', gene_to_search)
            #print('where !=0: ', list(np.argwhere(AM[gene_to_search] != 0)))
            for each in np.argwhere(AM[gene_to_search] != 0):
                if each[0] in temp_subnetwork:
                    continue
                else:
                    gene_to_append.append(each[0])
            #print('to_append: ', gene_to_append)
            temp_subnetwork.append(gene_to_append[0])
            gene_to_search = gene_to_append[0]
            
        temp_subnetwork = set(temp_subnetwork)
        AM_subnetworks.append(temp_subnetwork)
    #print('AM_subnetworks (',len(AM_subnetworks),') :',AM_subnetworks,'\n')
    '''Combining intersected sets.'''
    while Is_intersection(AM_subnetworks):
        output_subnetworks = []
        while len(AM_subnetworks) > 0:
            i = np.random.randint(len(AM_subnetworks))
            Is_isolated = True
            for j in range(0, len(AM_subnetworks)):
                if j == i:
                    continue
                else:
                    pass
                if len(AM_subnetworks[i].intersection(AM_subnetworks[j])) != 0:
                    output_subnetworks.append(AM_subnetworks[i].union(AM_subnetworks[j]))
                    Is_isolated = False
                    if i > j:
                        del AM_subnetworks[i]
                        del AM_subnetworks[j]
                    else:
                        del AM_subnetworks[j]
                        del AM_subnetworks[i]
                    break
                else:
                    pass
            if Is_isolated:
                output_subnetworks.append(AM_subnetworks[i])
                del AM_subnetworks[i]
            else:
                pass
        AM_subnetworks = output_subnetworks
    return output_subnetworks
