import os, sys, time, utilities
import numpy as np
import pandas as pd
import nibabel as nib
from itertools import product
from joblib import Parallel, delayed, dump, load
from brainspace.utils.parcellation import reduce_by_labels

def main(nboot, batch_size, batch_nr):
    # Load data
    wd = '/project/6050199/rhaast/03_Ongoing/hippocampal_perfusion'

    # Surface
    unfolded   = wd + '/resources/midthickness.L.unfolded.surf.gii'
    coords     = nib.load(unfolded).get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
    nvertices  = len(coords)

    # Atlas
    atlas      = wd + '/resources/BigBrain_ManualSubfieldsUnfolded_254x126.shape.gii'
    subfields  = nib.load(atlas).darrays[0].data
    labels     = ['Sub', 'CA1', 'CA2', 'CA3', 'CA4/DG']
    nsubfields = len(np.unique(subfields))

    # Subjects
    df = pd.read_csv(wd + '/config/participants.txt', dtype=str)
    subjects = df.participant_id.to_list() 
    subjects = [ s for s in subjects if s not in ['10'] ]

    # Hemispheres
    hemis = ['Lflip','R']

    # Runs
    runs = [str(r) for r in range(1,9)]

    # Load input, per subject
    fname = wd + '/results/surface_maps/sub-{0}/full_fit/sub-{0}_DIFF_{1}.native.shape.gii'
    nvolumes = len(nib.load(fname.format(subjects[0], hemis[0])).darrays)

    data_in = np.zeros((
        nvertices,
        len(subjects),
        len(hemis),
        nvolumes
    ))

    for s, subject in enumerate(subjects):
        for h, hemi in enumerate(hemis):
            data = nib.load(fname.format(subject, hemi))
            data = np.array([ data.darrays[x].data for x in range(0,len(data.darrays)) ]).T            
            data_in[:,s,h,:] = data

    nvolumes_per_run = nvolumes/len(runs)
    run_index = np.zeros((nvolumes))
    for i in range(0,len(runs)):
        start = int(i*nvolumes_per_run)
        end   = int(start+nvolumes_per_run)
        run_index[start:end] = i

    # Bootstrapping
    # Set seed to make sure we can reproduce
    np.random.seed(2)
    nboot = int(nboot)

    idx, idy = [], []
    b = 0
    while b < nboot:
        x = list(np.random.choice(range(len(runs)), size=len(runs), replace=False))

        if x not in idx:
            idx.append(x)     
            b += 1

    b = 0
    while b < nboot:
        y = list(np.random.choice(range(len(subjects)), size=len(subjects), replace=False))

        if y not in idy:
            idy.append(y)       
            b += 1

    # Temporary location for dumping data
    tmp_dir = os.environ['SLURM_TMPDIR']
    
    # Dump input to file for memory-efficiency
    input_filename_memmap = os.path.join(tmp_dir, f'data_memmap_{batch_nr}')
    dump(data_in, input_filename_memmap)
    input_memmap = load(input_filename_memmap, mmap_mode='r')

    # Dump output to file for memory-efficiency
    output_filename_memmap = os.path.join(tmp_dir, f'output_memmap_{batch_nr}')
    output_memmap = np.memmap(
        output_filename_memmap, dtype=float, mode='w+',
        shape=(len(subjects), len(runs), len(subfields), 4, int(batch_size))
    )

    # Start parallelization
    start = int(batch_size)*(int(batch_nr)-1)
    end   = start+int(batch_size)
    print(f'Processing bootstrap sample {start} to {end}')
    
    Parallel(n_jobs=-1)(
        delayed(utilities.bootstrap_iteration_full)(
            input_memmap, output_memmap, i, idx[j], idy[j], subfields, run_index
        ) for i, j in enumerate(np.arange(start,end))
    )   

    # Save to Numpy file
    bootstrap_matrix_full = np.array(output_memmap)

    np.savez(
        wd + f'/results/bootstrapping/bootstrap_matrix_full-fit_Nboot-{nboot}_batch-{batch_nr}_new_vertices.npz',
        bootstrap_matrix_full=bootstrap_matrix_full, 
        idx=idx, idy=idy
    )  


if __name__ == "__main__":
    import sys
    
    nboot      = sys.argv[1]
    batch_size = sys.argv[2]
    batch_nr   = sys.argv[3]
    
    print(f'Running bootstrap analyses using {nboot} samples, in batches of {batch_size}. Current batch: {batch_nr}')
    
    main(nboot, batch_size, batch_nr)