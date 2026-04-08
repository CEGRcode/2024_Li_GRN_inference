import os
import sys
import h5py
import numpy as np
from sklearn import preprocessing

input_dir = sys.argv[1] if len(sys.argv) > 1 else "./"

def binary_to_string(binary_data):
    if isinstance(binary_data, np.bytes_):
        binary_data = bytes(binary_data)
    binary_string = ''.join(format(byte, '08b') for byte in binary_data)
    binary_values = [binary_string[i:i + 8] for i in range(0, len(binary_string), 8)]
    characters = [chr(int(bv, 2)) for bv in binary_values]
    return ''.join(characters)

for filename in os.listdir(input_dir):
    if not filename.endswith('_BED.h5'):
        continue

    h5_file_path = os.path.join(input_dir, filename)

    with h5py.File(h5_file_path, 'r') as h5file:
        for keys in h5file:
            print(keys)
            WindowSize = 500
            length = 500
            any_CDT = [0 for _ in range(length)]

            for each in h5file[keys][:]:
                line = binary_to_string(each)
                if line.startswith('YORF'):
                    continue
                parts = line.split()
                for i in range(len(parts[2:])):
                    any_CDT[i] += float(parts[i + 2])

            any_CDT = preprocessing.normalize([np.array(any_CDT[1:-1])])[0]
            any_CDT = any_CDT[(length // 2 - WindowSize // 2):(length // 2 + WindowSize // 2)]

            output_h5 = os.path.join(input_dir, keys.split('_')[0] + '_BED_combined.h5')
            with h5py.File(output_h5, 'a') as h5file_write:
                data_as_lines = ['\t'.join(map(str, any_CDT))]
                data_array = np.array(data_as_lines, dtype='S')
                if keys in h5file_write:
                    del h5file_write[keys]
                h5file_write.create_dataset(keys, data=data_array, compression="gzip")
