# Prepare T2w data for mapping on unfolded hippocampus
rule warp_t2w_to_corobl_crop:
    input:
        nii = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-TSEavg_T2w.nii.gz'),
        init = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt'),
        xfm = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = config['template']['hires']
    output: 'results/maps/sub-{subject}/sub-{subject}_T2w_{hemi}.nii.gz'
    group: 'map_t2w'
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm} -t {input.init}"

# Flip left hemisphere
rule lr_flip_t2w:
    input: 'results/maps/sub-{subject}/sub-{subject}_T2w_L.nii.gz'
    output: 'results/maps/sub-{subject}/sub-{subject}_T2w_Lflip.nii.gz'
    group: 'map_t2w'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

# Sample T2w volume onto unfolded hippocampal surface
rule sample_t2w_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_T2w_{H}.nii.gz',
        ribbon = rules.extract_gm_ribbon.output,
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_T2w_{H}.native.shape.gii'
    group: 'map_t2w'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"