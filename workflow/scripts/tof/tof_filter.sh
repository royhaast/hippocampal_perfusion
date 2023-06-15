#!/bin/bash
in_tof=$1
in_file="$(basename "${in_tof}" .nii.gz)"  
out_dir="$(dirname "${in_tof}")" 

source /home/rhaast/venv/segmentator/bin/activate

pushd $out_dir
segmentator_filters --smoothing STEDI --nr_iterations 10 ${in_file}.nii.gz 
mv ${in_file}_STEDI_n10_s0pt5_r0pt5_g1.nii.gz ${in_file}_stedi.nii.gz
popd