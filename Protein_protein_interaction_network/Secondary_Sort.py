import numpy as np
import pandas as pd
import math
import os
import sys 
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
from sklearn import preprocessing
from scipy import linalg
from scipy.spatial import distance
from scipy.stats import norm, skew, ks_2samp, zscore
from scipy.integrate import simpson
import zipfile
import shutil
    
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

def calculate_auc(data):
    data = np.array(data)
    significant = data > 0
    regions = np.diff(np.concatenate(([0], significant.astype(int), [0])))
    start_indices = np.where(regions == 1)[0]
    end_indices = np.where(regions == -1)[0]
    total_area = 0
    for start, end in zip(start_indices, end_indices):
        region_area = simpson(data[start:end], dx=1)
        total_area += region_area
    return total_area

def get_largest_peak_width(data):
    peaks_width = []
    peak_width = 0
    for each in data:
        if each == 0:
            peaks_width.append(peak_width)
            peak_width = 0
            continue
        else:
            peak_width = peak_width + 1
    return max(peaks_width)

def get_num_of_intersec(folder_path, ID_, BED_ID_):
    forward_bg, reverse_bg = read_cdt(folder_path+'/{}_BED_{}_BAM_sense_SCALE.cdt'.format(BED_ID_, ID_), folder_path+'/{}_BED_{}_BAM_anti_SCALE.cdt'.format(BED_ID_, ID_))
    x = np.linspace(0, len(forward_bg), len(forward_bg))
    new_x = np.linspace(0, len(forward_bg), len(forward_bg)*100)
    new_forward_bg = np.interp(new_x, x, forward_bg)
    num_of_intersec_forward = []
    locs_of_intersec_forward = []
    for i in np.arange(0, max(new_forward_bg), max(new_forward_bg)/1000):
        temp_locs = np.argwhere((new_forward_bg >= i) & (new_forward_bg <= i+max(new_forward_bg)/1000)).tolist()
        num_of_intersec_forward.append(len(temp_locs))
        locs_of_intersec_forward.append(temp_locs)
    x = np.linspace(0, len(reverse_bg), len(reverse_bg))
    new_x = np.linspace(0, len(reverse_bg), len(reverse_bg)*100)
    new_reverse_bg = np.interp(new_x, x, reverse_bg)
    num_of_intersec_reverse = []
    locs_of_intersec_reverse = []
    for i in np.arange(0, max(new_reverse_bg), max(new_reverse_bg)/1000):
        temp_locs = np.argwhere((new_reverse_bg >= i) & (new_reverse_bg <= i+max(new_reverse_bg)/1000)).tolist()
        num_of_intersec_reverse.append(len(temp_locs))
        locs_of_intersec_reverse.append(temp_locs)
    for i in range(0, len(locs_of_intersec_forward)):
        locs_of_intersec_forward[i] = sum(locs_of_intersec_forward[i], [])
    for i in range(0, len(locs_of_intersec_reverse)):
        locs_of_intersec_reverse[i] = sum(locs_of_intersec_reverse[i], [])
    num_of_intersec_forward = (num_of_intersec_forward - np.min(num_of_intersec_forward)) / (np.max(num_of_intersec_forward) - np.min(num_of_intersec_forward))
    num_of_intersec_reverse = (num_of_intersec_reverse - np.min(num_of_intersec_reverse)) / (np.max(num_of_intersec_reverse) - np.min(num_of_intersec_reverse))
    forward_bg = (forward_bg - np.min(forward_bg)) / (np.max(forward_bg) - np.min(forward_bg))
    reverse_bg = (reverse_bg - np.min(reverse_bg)) / (np.max(reverse_bg) - np.min(reverse_bg))
    return forward_bg, reverse_bg, num_of_intersec_forward, num_of_intersec_reverse, locs_of_intersec_forward, locs_of_intersec_reverse

def main_loop(BED_ID):
    folder_path = './SCALE_cdt/Results/{}_BED'.format(BED_ID)

    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        files = os.listdir(folder_path)
    else:
        return
    matching_files = [file for file in files if file.endswith('.txt')]
    
    infile = open(folder_path+'/'+matching_files[0], 'r')
    Sample_IDs = []
    for line in infile:
        if line.split()[0] == 'YEP_ID':
            continue
        else:
            Sample_IDs.append(line.split()[0])
    infile.close()
    
    folder_path = './Temp_cdt/SCALE_cdt'
    
    forward_bg, reverse_bg, intersec_forward_bg, intersec_reverse_bg, locs_intersec_forward_bg, locs_of_intersec_reverse_bg = get_num_of_intersec(folder_path, BED_ID, BED_ID)
    Overall_Scores = []
    
    for each_ID in Sample_IDs:
        forward, reverse, intersec_forward, intersec_reverse, locs_intersec_forward, locs_intersec_reverse = get_num_of_intersec(folder_path, each_ID, BED_ID)
        peaks, _ = find_peaks(intersec_forward, prominence=0.05)

        if len(peaks) == 0:
            Peak_width_forward = 0
            Complexity_forward = 0
            Peak_height_forward = 0
            Num_Of_Peaks_forward = len(peaks)
        else:
            forward_filtered = [0 for i in range(0, len(forward))]
            for i in range(0, len(forward)):
                if forward[i] >= (peaks[-1]/1000)*max(forward):
                    forward_filtered[i] = forward[i]
                else:
                    continue
            peaks_expand = []
            for each in peaks:
                peaks_expand.append(each)
                for i in range(1, 1):
                    peaks_expand.append(each+i)
                    peaks_expand.append(each-i)
            
            intersec_forward_trim = []
            for i in range(0, len(intersec_forward)):
                if i in peaks_expand:
                    continue
                else:
                    intersec_forward_trim.append(intersec_forward[i])
            Peak_width_forward = get_largest_peak_width(forward_filtered)
            Complexity_forward = calculate_auc(intersec_forward_trim)
            Num_Of_Peaks_forward = len(peaks)

            if np.argwhere(intersec_forward == max(intersec_forward))[0][0]:
                Peak_height_forward = 1/(np.argwhere(intersec_forward == max(intersec_forward))[0][0])
            else:
                outfile = open('./Second_sort/{}_Second_Sort.txt'.format(BED_ID), 'a')
                outfile.write(str(each_ID)+'\t'+'Bad Quality'+'\n')
                outfile.close()
                continue
    
        peaks, _ = find_peaks(intersec_reverse, prominence=0.05)
        if len(peaks) == 0:
            Peak_width_reverse = 0
            Complexity_reverse = 0
            Num_Of_Peaks_reverse = len(peaks)
            Peak_height_reverse = 0
        else:
            reverse_filtered = [0 for i in range(0, len(reverse))]
            for i in range(0, len(reverse)):
                if reverse[i] >= (peaks[-1]/1000)*max(reverse):
                    reverse_filtered[i] = reverse[i]
                else:
                    continue
            peaks_expand = []
            for each in peaks:
                peaks_expand.append(each)
                for i in range(1, 1):
                    peaks_expand.append(each+i)
                    peaks_expand.append(each-i)
            
            intersec_reverse_trim = []
            for i in range(0, len(intersec_reverse)):
                if i in peaks_expand:
                    continue
                else:
                    intersec_reverse_trim.append(intersec_reverse[i])
            Peak_width_reverse = get_largest_peak_width(reverse_filtered)
            Complexity_reverse = calculate_auc(intersec_reverse_trim)
            Num_Of_Peaks_reverse = len(peaks)
    
            if np.argwhere(intersec_reverse == max(intersec_reverse))[0][0]:
                Peak_height_reverse = 1/(np.argwhere(intersec_reverse == max(intersec_reverse))[0][0])
            else:
                outfile = open('./Second_sort/{}_Second_Sort.txt'.format(BED_ID), 'a')
                outfile.write(str(each_ID)+'\t'+'Bad Quality'+'\n')
                outfile.close()
                continue

        Over_all_Score_forward = Complexity_forward*Peak_width_forward*Peak_height_forward
        Over_all_Score_reverse = Complexity_reverse*Peak_width_reverse*Peak_height_reverse
        Overall_Score = np.max([Over_all_Score_forward, Over_all_Score_reverse])
        Overall_Scores.append(Overall_Score)

        outfile = open('./Second_sort/{}_Second_Sort.txt'.format(BED_ID), 'a')
        outfile.write(str(each_ID)+'\t'+str(Num_Of_Peaks_forward)+'\t'+str(Num_Of_Peaks_reverse)+'\t'+str(np.round(Complexity_forward, 1))+'\t'+str(np.round(Complexity_reverse, 1))+'\t'+str(np.round(Peak_width_forward, 1))+'\t'+str(np.round(Peak_width_reverse, 1))+'\t'+str(np.round(Peak_height_forward, 6))+'\t'+str(np.round(Peak_height_reverse, 6))+'\t'+str(np.round(Overall_Score, 1))+'\n')
        outfile.close()
    return

directory = './SCALE_cdt/'

for filename in os.listdir(directory):
    # Check if the file is a zip file
    if filename.endswith('zipped_BED_BAM_sense.zip'):
        zip_file_path = os.path.join(directory, filename)
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            # Extract all contents to the specified directory
            zip_ref.extractall('./Temp_cdt/')
    else:
        continue
    main_loop(filename.split('_')[0])
    if os.path.exists('./Temp_cdt/') and os.path.isdir('./Temp_cdt/'):
        shutil.rmtree('./Temp_cdt/')
    else:
        pass
    os.makedirs('./Temp_cdt/')

