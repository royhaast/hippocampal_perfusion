# # Correct MP2RAGE data for B1+ using Sriranga's set up
# rule mp2rage_correction:
#     input: 
#         gradcorrect = config['data']['gradcorrect'],
#         mp2rage_correction = config['data']['mp2rage_correction'],
#     output: 'results/mp2rage_correction/sub-{subject}/anat/sub-{subject}_acq-MP2RAGE_run-01_corrT1_clean.nii.gz'
#     params:
#         script = '/project/6050199/rhaast/03_Ongoing/hippocampal_perfusion/workflow/scripts/mp2rage_correction/wrapper.sh'
#     group: 'preprocessing'    
#     threads: 8
#     resources:
#         mem_mb = 32000
#     log: 'logs/mp2rage_correction/sub-{subject}.log'
#     shell:
#         "bash {params.script} {wildcards.subject} {input.gradcorrect} {input.mp2rage_correction}"     

# rule tar_mp2rage_correction:
#     input: 'results/mp2rage_correction/sub-{subject}'
#     output: 'results/mp2rage_correction/sub-{subject}_mp2rage_corrected.tar.gz'
#     shell:
#         "cd results/mp2rage_correction && tar -czf `basename {output}` sub-{wildcards.subject}"

# Remove salt-n-peper background noise from MP2RAGE uni
# image
rule mprageise:
    input: unpack(collect_input)
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_run-01_corrUNI_clean_unbiased_clean.nii.gz'
    params:
        script = '/project/6050199/rhaast/03_Ongoing/hippocampal_perfusion/workflow/scripts/skullstripping/mprageise.sh'    
    group: 'preprocessing'
    singularity: config['singularity_fmriprep']    
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/mprageise/sub-{subject}.log'
    shell:
        "bash {params.script} -i {input.inv2} -u {input.t1w} -o `dirname {output}`" 

# Remove pads from images as these affect skullstripping
# efficicacy
rule get_pads_mask:
    input: rules.mprageise.output #unpack(collect_input)
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_run-01_corrUNI_clean_unbiased_clean_dePadded.nii.gz'
    params:
        script = '/project/6050199/rhaast/03_Ongoing/hippocampal_perfusion/workflow/scripts/skullstripping/PadsOff'
    group: 'preprocessing'
    singularity: config['singularity_fmriprep']    
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/padsoff/sub-{subject}.log'    
    shell:
        "bash {params.script} -i {input}"

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

# Remove non-brain tissue from MP2RAGE images using Sriranga's
# set up (SPM12 + CAT)
rule skull_stripping:
    input: rules.get_pads_mask.output #unpack(collect_input)
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
    output: 'results/freesurfer/sub-{subject}/scripts/recon-all.done'
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

# Run FreeSurfer's hippocampal segmentation pipeline
rule extract_freesurfer_hippocampi:
    input:
        freesurfer = rules.freesurfer.output,
        t2w = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-TSE_rec-avg.nii.gz')
    output:
        lh = 'results/freesurfer/sub-{subject}/lh.hippoAmygLabels-T1-T2w.v21.mgz',
        rh = 'results/freesurfer/sub-{subject}/rh.hippoAmygLabels-T1-T2w.v21.mgz'
    params:
        sd = '/scratch/rhaast/HPC_perfusion/results/freesurfer'    
    group: 'freesurfer'
    singularity: config['singularity_freesurfer']
    threads: 8
    resources:
        time = 300,
        mem_mb = 32000
    shell:
        "segmentHA_T2.sh sub-{wildcards.subject} {input.t2w} T2w 1 {params.sd}"   

# Prepare MP2RAGE data for mapping on unfolded hippocampus
# rule warp_t1_to_corobl_crop:
#     input:
#         nii = lambda wildcards: 'results/t1/sub-{s}/sub-{s}_{f}'.format(s=wildcards.subject, f=config['mp2rage_data'][wildcards.mp2rage_parameter]),
#         init = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
#         xfm = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
#         ref = join(config['hippunfold_dir'],config['template'])
#     output: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{hemi}.nii.gz'
#     group: 'map_t1' 
#     singularity: config['singularity_prepdwi']
#     threads: 8
#     resources:
#         mem_mb = 32000    
#     shell:
#         "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
#         "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm} -t {input.init}"  

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

rule lr_flip_t1:
    input: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_L.nii.gz',
    output: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_Lflip.nii.gz',
    group: 'map_t1'     
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

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

# Calculate myelin map based on T1w (UNI) and T2w surface
# maps
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

# Get retinotopic maps
rule retinotopy:
    input: 'results/freesurfer/sub-{subject}/scripts/recon-all.done'
    output: 
        lh = 'results/freesurfer/sub-{subject}/surf/lh.wang15_mplbl.mgz',
        rh = 'results/freesurfer/sub-{subject}/surf/rh.wang15_mplbl.mgz'
    shell: 
        """
        python -m neuropythy atlas --verbose sub-{wildcards.subject}
        """   

rule retinotopy_to_labels:
    input:
        lh = 'results/freesurfer/sub-{subject}/surf/lh.wang15_mplbl.mgz',
        rh = 'results/freesurfer/sub-{subject}/surf/rh.wang15_mplbl.mgz'    
    output:
    params:
        sd = 'results/freesurfer'    
    shell:
        """
        export roiname_array=(1 "V1v" "V1d" "V2v" "V2d" "V3v" "V3d" "hV4" "VO1" "VO2" "PHC1" "PHC2" \
            "TO2" "TO1" "LO2" "LO1" "V3B" "V3A" "IPS0" "IPS1" "IPS2" "IPS3" "IPS4" \
            "IPS5" "SPL1" "FEF")
        export subjid={wildcards.subject}

        for i in {1..25} ; do
            mri_cor2label --i {params.sd}/{wildcards.subject}/surf/rh.wang15_mplbl.mgz --id ${i} --l rh.wang15atlas.${roiname_array[${i}]}.label --surf ${subjid} rh inflated
        done
        """