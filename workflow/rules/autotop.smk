# Prepare manual segmentations to use as input to the autotop 
# command. This includes merging L and R hemispheres labels,
# relabeling to 'expected' labels, aligning to first T2w volume
# and minor label smoothing to remove lonely voxels.
rule combine_lr_manual_seg:
    input:
        lh = join(config['data']['manual_segs'],'sub-{subject}_hemi-L_hipp.nii.gz'),
        rh = join(config['data']['manual_segs'],'sub-{subject}_hemi-R_hipp.nii.gz'),
    output: join(config['data']['manual_segs'],'sub-{subject}_hemi-LR_hipp.nii.gz'),
    group: 'manual_seg'
    run: 
        import numpy as np
        import nibabel as nib

        L  = nib.load(input.lh)
        R  = nib.load(input.rh)
        LR = L.get_fdata()+R.get_fdata()

        img = nib.Nifti1Image(LR.astype(np.int8), header=L.header, affine=L.affine)
        nib.save(img, output[0])

rule relabel_lr_manual_seg:
    input: rules.combine_lr_manual_seg.output
    output: join(config['data']['manual_segs'],'sub-{subject}_hemi-LR_hipp_relabeled.nii.gz')
    group: 'manual_seg'    
    run:
        import numpy as np
        import nibabel as nib

        LR = nib.load(input[0])
        LR_labels = LR.get_fdata()
        LR_relabeled = np.zeros(LR_labels.shape)

        relabeling = [[1, 1],[2, 2],[20, 3],
                      [21, 4],[22, 5],[23, 6]]

        for labels in relabeling:
            LR_relabeled[LR_labels==labels[0]] = labels[1]

        img = nib.Nifti1Image(LR_relabeled.astype(np.int8), header=LR.header, affine=LR.affine)
        nib.save(img, output[0])

rule get_xfm_manual_seg_to_t2w:
    input:
        flo = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0.nii.gz'),
        ref = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-TSE_run-01_T2w.nii.gz')
    output:
        xfm_ras = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-ras_xfm.txt'),
        xfm_itk = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt'),
        warped = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0_space-refT2w.nii.gz')
    group: 'manual_seg'        
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000     
    shell:
        'reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras} -rigOnly -nac && '
        'c3d_affine_tool {output.xfm_ras} -oitk {output.xfm_itk}'

rule apply_xfm_manual_seg_to_t2w:
    input:
        nii = rules.relabel_lr_manual_seg.output,
        ref = rules.get_xfm_manual_seg_to_t2w.input.flo,
        xfm = rules.get_xfm_manual_seg_to_t2w.output.xfm_itk
    output: join(config['data']['manual_segs'],'sub-{subject}_desc-manualseg.nii.gz')
    group: 'manual_seg'    
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output} -r {input.ref} -t {input.xfm}'     

rule smooth_labels:
    input: rules.apply_xfm_manual_seg_to_t2w.output
    output: join(config['data']['manual_segs'],'sub-{subject}_desc-manualseg_dseg.nii.gz')
    group: 'manual_seg'    
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "c3d {input} -split -foreach -smooth 0.3mm -endfor -merge -o {output}"

# # Run autotop using preprocessed manual segmentations and T2w
# # images. The pipeline will average across different runs, 
# # crop to input space (i.e., coronal oblique 0.3 mm iso) and 
# # then unfold the hippocampus.
# rule autotop:
#     input:
#         bids = 'data/bids', #'/project/6050199/rhaast/99_tmp_scratch_backup/MSTRCHT_BIDSify/gradcorrect',
#         dseg = expand(join(config['data']['manual_segs'],'sub-{subject}_desc-manualseg_dseg.nii.gz'), subject=subjects), #rules.smooth_labels.output,
#         t1w = expand(join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-MP2RAGEdepadded_run-01_T1w.nii.gz'), subject=subjects)
#     output:
#         dseg = expand('results/autotop-dev/results/sub-{subject}/anat/sub-{subject}_space-T1w_desc-subfields_modality-segT2w_template-CITI168_dseg.nii.gz', subject=subjects),
#         #Lflip = directory(expand('results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-Lflip_space-CITI168corobl_desc-cropped_modality-segT2w_autotop', subject=subjects)),
#         #R = directory(expand('results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-R_space-CITI168corobl_desc-cropped_modality-segT2w_autotop', subject=subjects))
#     params:
#         out_dir = 'results/autotop-dev',
#         autotop = config['autotop'],
#         autotop_opts = '--derivatives /scratch/rhaast/HPC_perfusion/data/manual_segs --modality segT2w --filter_T1w acquisition=MP2RAGEdepadded',
#         snakemake_opts = '--cores 8 --use-singularity -p' #&> logs/autotop/sub-{subject}.log
#     group: 'autotop'        
#     threads: 8
#     resources:
#         time = 600,
#         mem_mb = 64000
#     shell:
#         "{params.autotop} {input.bids} {params.out_dir} participant {params.autotop_opts} {params.snakemake_opts}" #--participant_label {wildcards.subject} 