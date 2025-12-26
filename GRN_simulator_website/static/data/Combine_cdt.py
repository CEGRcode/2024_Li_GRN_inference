import os
import h5py
import numpy as np
import zipfile
import shutil

for filename in os.listdir('./'):
    if filename.endswith('.zip'):
        # Define the full path to the zip file
        zip_file_path = os.path.join('./', filename)
        
        # Define the extraction directory
        extract_dir = './SCALE_cdt/'
        
        # Create the extraction directory if it does not exist
        if not os.path.exists(extract_dir):
            os.makedirs(extract_dir)
        
        # Unzip the file
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

        # Create a new HDF5 file
        output_h5 = filename.split('_')[0] + '_BED.h5'
        with h5py.File(output_h5, 'w') as h5file:
            # Loop through each .cdt file in the directory
            for cdt_file in os.listdir('./SCALE_cdt/SCALE_cdt/'):
                if cdt_file.endswith('.cdt'):
                    file_path = os.path.join('./SCALE_cdt/SCALE_cdt/', cdt_file)
                    # Read the content of the .cdt file
                    with open(file_path, 'r') as file:
                        content = file.readlines()
                        
                        # Create a dataset in the HDF5 file for this .cdt file
                        dataset_name = os.path.splitext(cdt_file)[0]
                        
                        # Store the content as a dataset in the HDF5 file
                        h5file.create_dataset(dataset_name, data=np.array(content, dtype='S'), compression="gzip")

                        # Add metadata to identify the source file
                        h5file[dataset_name].attrs['source_file'] = cdt_file

        print(f"All .cdt files have been concatenated into {output_h5}.")
    else:
        pass
    if os.path.isdir('./SCALE_cdt/'):
        shutil.rmtree('./SCALE_cdt/')
    else:
        pass
