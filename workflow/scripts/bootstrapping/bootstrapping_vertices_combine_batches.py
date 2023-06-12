import os, sys, time, h5py, utilities
import numpy as np
import pandas as pd
import nibabel as nib
from itertools import product
from joblib import Parallel, delayed, dump, load
from brainspace.utils.parcellation import reduce_by_labels

def main(batch_size, batches):
    wd       = '/project/6050199/rhaast/03_Ongoing/hippocampal_perfusion'
    out_file = wd + f'/results/bootstrapping/bootstrap_matrix_full-fit_Nboot-1000_batch-avg_new_vertices.hdf5'
    
    f = h5py.File(out_file, "w")
    
    for i, ii in enumerate(batches):
        data = np.load(
            wd + f'/results/bootstrapping/bootstrap_matrix_full-fit_Nboot-1000_batch-{ii}_new_vertices.npz'
        )

        batch_data = data['bootstrap_matrix_full']
        
        if i == 0:       
            bootstrap_matrix = f.create_dataset(
                "bootstrap_matrix",
                (
                    batch_data.shape[0],
                    batch_data.shape[1],
                    batch_data.shape[2],
                    batch_data.shape[3],
                    batch_size*len(batches)), dtype='f'
            )

        print(f'Loading bootstrap_matrix_full-fit_Nboot-1000_batch-{ii}_new_vertices.npz ...')            
            
        start = i*batch_size
        end   = start+batch_size
        bootstrap_matrix[:,:,:,:,start:end] = batch_data

if __name__ == "__main__":
    import sys
    
    batch_size = int(sys.argv[1])
    batches    = [int(x) for x in sys.argv[2].split(',')]
    
    print(f'Combining {batch_size*len(batches)} samples from {len(batches)} batches')
    
    main(batch_size, batches)