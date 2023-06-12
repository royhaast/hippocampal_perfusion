# Convert FreeSurfer T1 to Nifti
rule mri_convert_freesurfer:
    input: 
        unpack(collect_input),
        fs = 'results/freesurfer/sub-{subject}/mri/T1.mgz'
    output: 'results/freesurfer/sub-{subject}/mri/T1.nii.gz'
    group: 'map_b1map'
    singularity: config['singularity_freesurfer']   
    shell:
        "mri_convert {input.t1w} {output} -nc -rl {input.fs}" 

# Coregister lores MP2RAGE to hires MP2RAGE
rule coregister_lores_mp2rage:
    input:
        src = 'results/freesurfer/sub-{subject}/mri/T1.nii.gz',
        ref = lambda wildcards: 'data/preprocessed/sub-{subject}/MP2RAGE/sub-{subject}_acq-hiresMP2RAGE_{file}'.format(
            subject=wildcards.subject, file=config['mp2rage_data']['T1w']),
    output:
        xfm = 'results/maps/sub-{subject}/sub-{subject}_from-loresMP2RAGE-to-hiresMP2RAGE_type-ras_xfm.txt',
        # nii = 'results/maps/sub-{subject}/sub-{subject}_acq-loresMP2RAGE_space-hiresMP2RAGE.nii.gz'
    group: 'map_b1map'
    singularity: config['singularity_prepdwi'] 
    shell:   
        "reg_aladin -flo {input.src} -ref {input.ref} -aff {output.xfm} -rigOnly -nac"

# Convert transforms to correct format
rule convert_ras_xfm_sa2rage:
    input:
        xfm = 'results/maps/sub-{subject}/sub-{subject}_from-loresMP2RAGE-to-hiresMP2RAGE_type-ras_xfm.txt',    
        src = 'results/freesurfer/sub-{subject}/mri/T1.nii.gz',
        ref = lambda wildcards: 'data/preprocessed/sub-{subject}/MP2RAGE/sub-{subject}_acq-hiresMP2RAGE_{file}'.format(
            subject=wildcards.subject, file=config['mp2rage_data']['T1w']),
    output: 'results/maps/sub-{subject}/sub-{subject}_from-loresMP2RAGE-to-hiresMP2RAGE_type-itk_xfm.txt',
    group: 'map_b1map'    
    singularity: config['singularity_prepdwi'] 
    shell:
        "c3d_affine_tool {input.xfm} -oitk {output}"

# Single transform B1+ map to coronal oblique space
rule warp_b1map_to_corobl_crop:
    input:
        nii = lambda wildcards: 'data/preprocessed/sub-{subject}/SA2RAGE/sub-{subject}_acq-SA2RAGE_{file}'.format(
            subject=wildcards.subject, file=config['sa2rage_data']['B1map']),
        xfm2hires = 'results/maps/sub-{subject}/sub-{subject}_from-loresMP2RAGE-to-hiresMP2RAGE_type-itk_xfm.txt',
        xfm2tse = 'results/maps/sub-{subject}/sub-{subject}_from-hiresMP2RAGE-to-TSE_type-itk_xfm.txt',
        xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = config['template']['hires']
    output: 'results/maps/sub-{subject}/sub-{subject}_B1map_{hemi}.nii.gz'
    group: 'map_b1map' 
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref} -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse} -t {input.xfm2hires}" 

# Flip left hemisphere
rule lr_flip_b1map:
    input: 'results/maps/sub-{subject}/sub-{subject}_B1map_L.nii.gz',
    output: 'results/maps/sub-{subject}/sub-{subject}_B1map_Lflip.nii.gz',
    group: 'map_b1map'     
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

# Sample onto unfolded hippocampal surface
rule sample_b1map_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_B1map_{H}.nii.gz',
        ribbon = rules.extract_gm_ribbon.output,
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_B1map_{H}.native.shape.gii'
    group: 'map_b1map'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"
