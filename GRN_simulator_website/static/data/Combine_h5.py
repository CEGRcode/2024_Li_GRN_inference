import os
import h5py
import numpy as np
import zipfile
import shutil
from sklearn import preprocessing

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

for filename in os.listdir('./'):
    if filename.endswith('_BED.h5'):
        # Define the full path to the h5 file
        h5_file_path = os.path.join('./', filename)
    else:
        continue

    with h5py.File(h5_file_path, 'r') as h5file:

        for keys in h5file:
            print(keys)
            WindowSize = 500 # compensate the removed double tails in return. Adjust the x-axis accordingly.
            length = 500
            any_CDT = [0 for i in range(0, length)]
            for each in h5file[keys][:]:
                if binary_to_string(each).startswith('YORF'):
                    pass
                else:
                    for i in range(0, len(each.split()[2:])):
                        any_CDT[i] = any_CDT[i] + float(each.split()[i+2])
            any_CDT = preprocessing.normalize([np.array(any_CDT[1:-1])])[0]
            any_CDT = any_CDT[(int(length/2)-int(WindowSize/2)):(int(length/2)+int(WindowSize/2))]

            output_h5 = keys.split('_')[0] + '_BED_combined.h5'
            with h5py.File(output_h5, 'a') as h5file_write:
                data_as_lines = ['\t'.join(map(str, any_CDT))]  # Creates a list with one line where elements are tab-separated

                # Convert to a NumPy array, similar to how you handled 'content'
                data_array = np.array(data_as_lines, dtype='S')  # Store it as byte strings

                # Create the dataset in the HDF5 file
                h5file_write.create_dataset(keys, data=data_array, compression="gzip")
