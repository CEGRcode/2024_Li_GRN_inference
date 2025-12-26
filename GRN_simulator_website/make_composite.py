import os
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
import numpy as np
from sklearn import preprocessing
import colorsys

def generate_distinct_colors(n):
    colors = []
    for i in range(n):
        hue = i / n
        lightness = (50 + 20 * (i % 2)) / 100
        saturation = 90 / 100
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
        colors.append(rgb)
    return colors

# Generate 10 distinct colors
colors = generate_distinct_colors(10)


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

YEP_GO = {}
YEP_GO_file = open('./static/data/YEP_GO.txt', 'r')
for line in YEP_GO_file:
    YEP_GO[line.split('\t')[0].upper()] = line.split('\t')[1][:-1]
YEP_GO_file.close()

SampleKeyFile = open('./static/data/SupplementaryData-Table4_Sample-Key_tabular.tab', 'r')
Sample_Keys = {}
for line in SampleKeyFile:
    Sample_Keys[int(line.split()[4])] = line.split()[5]
SampleKeyFile.close()

Selected_Sample_Codes = ['8599', '18459']
BED_ID = '14619'
colors = generate_distinct_colors(len(Selected_Sample_Codes))
forward_bg, reverse_bg = read_cdt('./static/data/unzipped/{}_BED_{}_BAM_sense_SCALE.cdt'.format(BED_ID, BED_ID), './static/data/unzipped/{}_BED_{}_BAM_anti_SCALE.cdt'.format(BED_ID, BED_ID))

plt.plot(range(-250, len(forward_bg)-250), gaussian_filter(forward_bg, sigma=0), c='grey', alpha=0.5, linestyle='--', label=BED_ID+'_{}'.format(Sample_Keys[int(BED_ID)])+'_sense')
plt.plot(range(-250, len(reverse_bg)-250), -gaussian_filter(reverse_bg, sigma=0), c='grey', alpha=0.5, linestyle='--', label=BED_ID+'_{}'.format(Sample_Keys[int(BED_ID)])+'_anti')
i = 0
for each in Selected_Sample_Codes:
    file_path = os.path.join('./static/data/unzipped', '{}_BED_{}_BAM_sense_SCALE.cdt'.format(BED_ID, each))
    forward_CDT, reverse_CDT = read_cdt(file_path, file_path.replace("sense", "anti"))
    plt.plot(range(-250, len(forward_CDT)-250), gaussian_filter(forward_CDT, sigma=0), c=colors[i], label=each+'_{}'.format(Sample_Keys[int(each)])+'_sense')
    plt.plot(range(-250, len(reverse_CDT)-250), -gaussian_filter(reverse_CDT, sigma=0), c=colors[i], label=each+'_{}'.format(Sample_Keys[int(each)])+'_anti')
    i=i+1
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=len(Selected_Sample_Codes)+1, fontsize=7.5)
plt.savefig('./static/data/BED_{}_BAM_{}.png'.format(str(BED_ID)+'_BED', each))
plt.close()