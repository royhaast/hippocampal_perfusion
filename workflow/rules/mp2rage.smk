# Correct MP2RAGE data for B1+
rule mp2rage_correction:
    input: 
        gradcorrect = config['data']['gradcorrect']
    output: 'results/mp2rage_correction/sub-{subject}/anat/sub-{subject}_acq-MP2RAGE_run-01_corrT1_clean.nii.gz'
    params:
        mp2rage_correction = config['data']['mp2rage_correction'],
        script = 'workflow/scripts/mp2rage_correction/wrapper.sh'
    group: 'preprocessing'    
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/mp2rage_correction/sub-{subject}.log'
    shell:
        "bash ./{params.script} {wildcards.subject} {input.gradcorrect} {params.mp2rage_correction}"     

# Remove salt-n-peper background noise from MP2RAGE uni image
rule mprageise:
    input: unpack(collect_input)
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_run-01_corrUNI_clean_unbiased_clean.nii.gz'
    params:
        script = 'workflow/scripts/skullstripping/mprageise.sh'    
    group: 'preprocessing'
    singularity: config['singularity_fmriprep']    
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/mprageise/sub-{subject}.log'
    shell:
        "bash ./{params.script} -i {input.inv2} -u {input.t1w} -o `dirname {output}`" 

# Remove pads from images as these affect skullstripping efficicacy
rule get_pads_mask:
    input: rules.mprageise.output
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_run-01_corrUNI_clean_unbiased_clean_dePadded.nii.gz'
    params:
        script = 'scripts/skullstripping/PadsOff'
    group: 'preprocessing'
    singularity: config['singularity_fmriprep']    
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/padsoff/sub-{subject}.log'    
    shell:
        "bash {params.script} -i {input}"

# Copy to gradient corrected data folder to be used by BIDS apps
rule copy_depadded:
    input:
        nii = rules.get_pads_mask.output,
        json = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-MP2RAGE_run-01_UNI.json')
    output:
        nii = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-MP2RAGEdepadded_run-01_T1w.nii.gz'),
        json = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-MP2RAGEdepadded_run-01_T1w.json')
    group: 'preprocessing'        
    shell:
        "cp {input.nii} {output.nii} && cp {input.json} {output.json}"

# Remove non-brain tissue
rule skull_stripping:
    input: rules.get_pads_mask.output
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_mask.nii.gz'
    params:
        script = 'scripts/skullstripping/run_cat12.m'
    group: 'preprocessing'
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/skullstripping/sub-{subject}.log'
    shell:
        "bash scripts/skullstripping/skullstrip.sh {params.script} {input} `realpath {output}` " #&> {log} {params.out_dir}  

# Apply the brain mask
rule apply_brain_mask:
    input: 
        t1w = rules.mprageise.output,
        brain_mask = rules.skull_stripping.output
    output: 'results/freesurfer/sub-{subject}/mri/orig/001.mgz'
    group: 'preprocessing'
    singularity: config['singularity_freesurfer']
    shell:
        "mri_mask {input.t1w} {input.brain_mask} {output}"    

# Run FreeSurfer recon-all pipeline
rule freesurfer:
    input: rules.apply_brain_mask.output
    output: 
        reconall = 'results/freesurfer/sub-{subject}/scripts/recon-all.done',
        T1 = 'results/freesurfer/sub-{subject}/mri/T1.mgz',
        sphere = expand('results/freesurfer/sub-{{subject}}/surf/{hemi}.sphere.reg', hemi=['lh','rh']),
        pial = expand('results/freesurfer/sub-{{subject}}/surf/{hemi}.pial', hemi=['lh','rh']),
        white = expand('results/freesurfer/sub-{{subject}}/surf/{hemi}.white', hemi=['lh','rh']),
    params:
        sd = 'results/freesurfer'
    group: 'freesurfer'
    singularity: config['singularity_freesurfer']
    threads: 8
    resources:
        time = 300,
        mem_mb = 32000
    shell:
        "export SUBJECTS_DIR={params.sd} && "
        "recon-all -all -s sub-{wildcards.subject} -no-wsgcaatlas -notal-check -threads 8" 

# Prepare MP2RAGE data for mapping on unfolded hippocampus
rule convert_fsl_xfm_mp2rage:
    input:
        xfm = lambda wildcards: 'data/preprocessed/sub-{subject}/MP2RAGE/sub-{subject}_acq-hiresMP2RAGE_{file}'.format(
            subject=wildcards.subject, file=config['xfms']['hiresMP2RAGE']),
        ref = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0.nii.gz'),
        src = lambda wildcards: 'data/preprocessed/sub-{subject}/MP2RAGE/sub-{subject}_acq-hiresMP2RAGE_{file}'.format(
            subject=wildcards.subject, file=config['mp2rage_data']['T1'])
    output: 'results/maps/sub-{subject}/sub-{subject}_from-hiresMP2RAGE-to-TSE_type-itk_xfm.txt'
    group: 'map_t1'    
    singularity: config['singularity_prepdwi'] 
    shell:
        "c3d_affine_tool -ref {input.ref} -src {input.src} {input.xfm} -fsl2ras -oitk {output}"

# Resample MP2RAGE data to coronal oblique space
rule warp_t1_to_corobl_crop:
    input:
        nii = lambda wildcards: 'data/preprocessed/sub-{subject}/MP2RAGE/sub-{subject}_acq-hiresMP2RAGE_{file}'.format(
            subject=wildcards.subject, file=config['mp2rage_data'][wildcards.mp2rage_parameter]),
        xfm2tse = 'results/maps/sub-{subject}/sub-{subject}_from-hiresMP2RAGE-to-TSE_type-itk_xfm.txt',
        xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = config['template']['hires']
    output: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{hemi}.nii.gz'
    group: 'map_t1' 
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref} -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse}" 

# Flip left hemisphere data
rule lr_flip_t1:
    input: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_L.nii.gz',
    output: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_Lflip.nii.gz',
    group: 'map_t1'     
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

# Map to unfolded surface
rule sample_t1_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{H}.nii.gz',
        ribbon = rules.extract_gm_ribbon.output,
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{H}.native.shape.gii'
    group: 'map_t1'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"

# Calculate myelin map based on T1w (UNI) and T2w surface maps
rule calculate_myelin_map:
    input:
        t1w = 'results/surface_maps/sub-{subject}/sub-{subject}_T1w_{H}.native.shape.gii',
        t2w = 'results/surface_maps/sub-{subject}/sub-{subject}_T2w_{H}.native.shape.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_myelin_{H}.native.shape.gii'    
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -metric-math '(t1w/t2w)' {output} -var t1w {input.t1w} -var t2w {input.t2w} -fixnan 0"