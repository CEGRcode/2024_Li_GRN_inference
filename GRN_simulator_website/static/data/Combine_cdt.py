import os
import sys
import h5py
import numpy as np
import zipfile
import shutil

input_dir = sys.argv[1] if len(sys.argv) > 1 else "./"
output_dir = sys.argv[2] if len(sys.argv) > 2 else "./"

for filename in os.listdir(input_dir):
    if filename.endswith('.zip'):
        zip_file_path = os.path.join(input_dir, filename)
        extract_dir = os.path.join(output_dir, 'SCALE_cdt')

        if os.path.exists(extract_dir):
            shutil.rmtree(extract_dir)
        os.makedirs(extract_dir)

        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

        output_h5 = os.path.join(output_dir, filename.split('_')[0] + '_BED.h5')
        with h5py.File(output_h5, 'w') as h5file:
            for cdt_file in os.listdir(os.path.join(extract_dir, 'SCALE_cdt')):
                if cdt_file.endswith('.cdt'):
                    file_path = os.path.join(extract_dir, 'SCALE_cdt', cdt_file)
                    with open(file_path, 'r') as file:
                        content = file.readlines()

                    dataset_name = os.path.splitext(cdt_file)[0]
                    h5file.create_dataset(
                        dataset_name,
                        data=np.array(content, dtype='S'),
                        compression="gzip"
                    )
                    h5file[dataset_name].attrs['source_file'] = cdt_file

        print(f"All .cdt files have been concatenated into {output_h5}.")

        if os.path.isdir(extract_dir):
            shutil.rmtree(extract_dir)
