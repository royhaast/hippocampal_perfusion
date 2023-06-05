# Prepare ASL data for mapping on unfolded hippocampus
ants_dimensions = { 'MEAN': 3, 'PWI': 3, 'CBF': 3, 'tSNR': 3, 'DIFF': 4 }

rule convert_fsl_xfm_tse:
    input:
        xfm = lambda wildcards: 'data/tar/2022_03/ASL/FLIRT_MATS/sub-{subject}_{file}'.format(
            subject=wildcards.subject, file=config['xfms']['TSE']),            
        ref = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0.nii.gz'),
        src = lambda wildcards: config['perfusion_data']['per_run']['PWI'].format(subject=wildcards.subject, run=1)
    output: 'results/maps/sub-{subject}/sub-{subject}_from-ASL-to-TSE_type-itk_xfm.txt'
    group: 'map_perfusion_hpc'    
    singularity: config['singularity_prepdwi'] 
    shell:
        "c3d_affine_tool -ref {input.ref} -src {input.src} {input.xfm} -fsl2ras -oitk {output}"

# Proces individual runs first
rule warp_perf_to_corobl_crop_per_run:
    input:
        nii = lambda wildcards: config['perfusion_data']['per_run'][wildcards.asl_parameter].format(subject=wildcards.subject, run=wildcards.run),
        xfm2tse = 'results/maps/sub-{subject}/sub-{subject}_from-ASL-to-TSE_type-itk_xfm.txt',
        xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = lambda wildcards: config['template'][wildcards.tpl_res]
    output: 'results/maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_tpl-{tpl_res}_{asl_parameter}_{hemi}.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']   
    shell:
        """
        ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}
        antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref} \
            -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse}
        """

rule lr_flip_perf_per_run:
    input: 'results/maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_tpl-{tpl_res}_{asl_parameter}_L.nii.gz'
    output: 'results/maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_tpl-{tpl_res}_{asl_parameter}_Lflip.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o {output}"

rule sample_perf_hippocampus_per_run:
    input:
        nii = 'results/maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_tpl-{tpl_res}_{asl_parameter}_{H}.nii.gz',
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
        ribbon = rules.extract_gm_ribbon.output,
    output: 'results/surface_maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_tpl-{tpl_res}_{asl_parameter}_{H}.native.shape.gii'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_connectomewb']   
    shell:
        """
        if [ {wildcards.tpl_res} = 'hires' ] ; then 
            wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} \
                -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}
        else
            wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} \
                -trilinear
        fi
        """
        
# rule reduce_nii_across_runs:
#     input: expand('results/maps/sub-{{subject}}/run-{run}/sub-{{subject}}_run-{run}_tpl-{tpl_res}_{{asl_parameter}}_{{H}}.nii.gz', run=range(1,9)),
#     output: 'results/maps/sub-{subject}/sub-{subject}_tpl-{tpl_res}_{asl_parameter}_{H}.nii.gz'
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_prepdwi']   
#     shell:
#         "AverageImages 3 {output} 0 {input}"  

# Then process full fit first 
rule warp_perf_to_corobl_crop_full_fit:
    input:
        nii = lambda wildcards: config['perfusion_data']['full_fit'][wildcards.asl_parameter].format(subject=wildcards.subject),
        xfm2tse = 'results/maps/sub-{subject}/sub-{subject}_from-ASL-to-TSE_type-itk_xfm.txt',
        xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = lambda wildcards: config['template'][wildcards.tpl_res]
    output: 'results/maps/sub-{subject}/full_fit/sub-{subject}_tpl-{tpl_res}_{asl_parameter}_{hemi}.nii.gz'  
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']   
    shell:
        """
        ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}
        antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref} \
            -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse}
        """

rule lr_flip_perf_full_fit:
    input: 'results/maps/sub-{subject}/full_fit/sub-{subject}_tpl-{tpl_res}_{asl_parameter}_L.nii.gz'
    output: 'results/maps/sub-{subject}/full_fit/sub-{subject}_tpl-{tpl_res}_{asl_parameter}_Lflip.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o {output}"

rule sample_perf_hippocampus_full_fit:
    input:
        nii = 'results/maps/sub-{subject}/full_fit/sub-{subject}_tpl-{tpl_res}_{asl_parameter}_{H}.nii.gz',
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
        ribbon = rules.extract_gm_ribbon.output,
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_tpl-{tpl_res}_{asl_parameter}_{H}.native.shape.gii'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_connectomewb']   
    shell:
        """
        if [ {wildcards.tpl_res} = 'hires' ] ; then 
            wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} \
                -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}
        else
            wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} \
                -trilinear
        fi
        """





# # Map onto surface per run, then average across runs in surface space
# rule sample_perf_hippocampus_per_run:
#     input:
#         nii = 'results/maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_{asl_parameter}_{H}.nii.gz',
#         inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
#         midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
#         outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
#         ribbon = rules.extract_gm_ribbon.output,
#     output: 'results/surface_maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_{asl_parameter}_{H}.native.shape.gii'
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_connectomewb']   
#     shell:
#         "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"

# rule merge_gii_across_runs:
#     input: expand('results/surface_maps/sub-{{subject}}/run-{run}/sub-{{subject}}_run-{run}_{{asl_parameter}}_{{H}}.native.shape.gii', run=['01','02','03','04','05','06','07','08']),
#     output: 'results/surface_maps/sub-{subject}/sub-{subject}_run-concat_{asl_parameter}_{H}.native.shape.gii'
#     params:
#         merge_cmd = construct_merge_cmd    
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_connectomewb']   
#     shell:
#         "{params.merge_cmd}"

# rule reduce_gii_across_runs:
#     input: 'results/surface_maps/sub-{subject}/sub-{subject}_run-concat_{asl_parameter}_{H}.native.shape.gii'
#     output: 'results/surface_maps/sub-{subject}/sub-{subject}_run-avg_{asl_parameter}_{H}.native.shape.gii'
#     singularity: config['singularity_connectomewb']
#     group: 'map_perfusion_hpc'      
#     shell:
#         "wb_command -metric-reduce {input} MEAN {output} -only-numeric"

# Now for ASL difference timeseries, per run and concatenate
rule split_timeseries:
    input: lambda wildcards: config['perfusion_data']['timeseries']['DIFF'].format(subject=wildcards.subject),
    output: 'results/maps/sub-{subject}/full_fit/native/sub-{subject}_vol0000.nii.gz' 
    params:
        prefix = 'results/maps/sub-{subject}/full_fit/native/sub-{subject}_vol'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']         
    shell:
        "fslsplit {input} {params.prefix} -t"

rule warp_perf_to_corobl_crop_timeseries:
    input:
        nii = 'results/maps/sub-{subject}/full_fit/native/sub-{subject}_vol0000.nii.gz',
        xfm2tse = 'results/maps/sub-{subject}/sub-{subject}_from-ASL-to-TSE_type-itk_xfm.txt',
        xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = config['template']['lores']
    output: 'results/maps/sub-{subject}/full_fit/tse/{hemi}/sub-{subject}_vol0000.nii.gz' 
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']   
    shell:
        """
        hemi={wildcards.hemi}
        in_dir=`dirname {input.nii}`
        out_dir=`dirname {output}`
        for file in $in_dir/*.nii.gz ; do
            out_file=`basename $file`
            antsApplyTransforms -d 3 --interpolation Linear -i $file -o $out_dir/$out_file -r {input.ref}  -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse}
        done
        """

rule lr_flip_perf_timeseries:
    input: 'results/maps/sub-{subject}/full_fit/tse/L/sub-{subject}_vol0000.nii.gz' 
    output: 'results/maps/sub-{subject}/full_fit/tse/Lflip/sub-{subject}_vol0000.nii.gz' 
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']   
    shell:
        """
        in_dir=`dirname {input}`
        out_dir=`dirname {output}`
        for file in $in_dir/*.nii.gz ; do
            out_file=`basename $file`
            c4d $file -flip x -o $out_dir/$out_file
        done
        """

checkpoint sample_perf_hippocampus_timeseries:
    input:
        nii = 'results/maps/sub-{subject}/full_fit/tse/{H}/sub-{subject}_vol0000.nii.gz',
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
        ribbon = rules.extract_gm_ribbon.output,
    output: 'results/surface_maps/sub-{subject}/full_fit/{H}/sub-{subject}_vol0000.shape.gii'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_connectomewb']     
    shell:
        """
        in_dir=`dirname {input.nii}`
        out_dir=`dirname {output}`
        exclude="0098 0197 0296 0395 0494 0593 0692 0791"
        for file in $in_dir/*.nii.gz ; do
            fname=`basename $file .nii.gz`
            vol=${{fname:10:4}}   

            if [[ $exclude != *$vol* ]] ; then
                out_file=${{fname}}.shape.gii
                wb_command -volume-to-surface-mapping $file {input.midthickness} $out_dir/$out_file -trilinear
            fi
        done
        """    

rule concatenate_perf_hippocampus_timeseries_gii:
    input:
        'results/surface_maps/sub-{subject}/full_fit/{H}/sub-{subject}_vol0000.shape.gii',
        lambda wildcards: sorted(glob('results/surface_maps/sub-{subject}/full_fit/{H}/sub-{subject}_vol*.shape.gii'.format(subject=wildcards.subject, H=wildcards.H)))
    output: 'results/surface_maps/sub-{subject}/full_fit/sub-{subject}_DIFF_{H}.native.shape.gii'
    params:
        merge_cmd = construct_merge_cmd    
    group: 'map_perfusion_hpc'
    singularity: config['singularity_connectomewb']   
    shell:
        "{params.merge_cmd}"    


# rule merge_timeseries:
#     input: 'results/maps/sub-{subject}/full_fit/tse/{hemi}/sub-{subject}_vol0000.nii.gz' 
#     output: 'results/maps/sub-{subject}/full_fit/tse/{hemi}/sub-{subject}_DIFF_{hemi}.nii.gz' 
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_prepdwi']
#     threads: 8
#     resources:
#         mem_mb = 32000              
#     shell:
#         """
#         out_dir=`dirname {output}`
#         fslmerge -t {output} `ls $out_dir/*.nii.gz`
#         """  

# rule warp_perf_to_corobl_crop_timeseries:
#     input:
#         nii = 'results/maps/sub-{subject}/full_fit/native/sub-{subject}_run-{run}_vol0000.nii.gz', 
#         xfm2tse = 'results/maps/sub-{subject}/run-{run}/sub-{subject}_run-{run}_from-ASL-to-TSE_type-itk_xfm.txt',
#         xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
#         xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
#         ref = join(config['hippunfold_dir'],config['template'])
#     output: 'results/maps/sub-{subject}/run-{run}/tse/{hemi}/sub-{subject}_run-{run}_vol0000.nii.gz'
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_prepdwi']       
#     shell:
#         """
#         hemi={wildcards.hemi}
#         in_dir=`dirname {input.nii}`
#         out_dir=`dirname {output}`
#         for file in $in_dir/*.nii.gz ; do
#             out_file=`basename $file`
#             antsApplyTransforms -d 3 --interpolation Linear -i $file -o $out_dir/$out_file -r {input.ref}  -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse} ; 
#         done
#         """

 

# # Prepare individual ASL runs for mapping on unfolded hippocampus
# rule warp_run_mean_to_corobl_crop:
#     input:
#         nii = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_acq-ASL_run-{run}_moCorr_sDC_PWI_Tmean.nii.gz',
#         xfm2anat = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_M0-to-MP2RAGE_BBR_fsl2greedy2ants_Apply1st.txt',
#         xfm2tse = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_MP2RAGE-to-TSE_final_greedy2ants_Apply2nd.txt',
#         xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
#         xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
#         ref = join(config['hippunfold_dir'],config['template'])
#     output: 'results/perfusion_preprocessing/sub-{subject}/sub-{subject}_run-{run}_mean_{hemi}.nii.gz'
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_prepdwi']   
#     shell:
#         "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
#         "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse} -t {input.xfm2anat}"

# rule lr_flip_run_mean:
#     input: 'results/perfusion_preprocessing/sub-{subject}/sub-{subject}_run-{run}_mean_L.nii.gz'
#     output: 'results/perfusion_preprocessing/sub-{subject}/sub-{subject}_run-{run}_mean_Lflip.nii.gz'
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_prepdwi']
#     shell:
#         "c3d {input} -flip x -o  {output}"

# rule merge_run_means:
#     input: expand('results/perfusion_preprocessing/sub-{{subject}}/sub-{{subject}}_run-{run}_mean_{{H}}.nii.gz', run=[str(r).zfill(2) for r in range(1,9)])
#     output: 'results/maps/sub-{subject}/sub-{subject}_PWI_means_{H}.nii.gz'
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_prepdwi']
#     shell:
#         "fslmerge -t {output} {input}"    

# rule sample_run_means_hippocampus:
#     input:
#         nii = 'results/maps/sub-{subject}/sub-{subject}_PWI_means_{H}.nii.gz',
#         ribbon = rules.extract_gm_ribbon.output,
#         inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
#         midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
#         outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
#     output: 'results/surface_maps/sub-{subject}/sub-{subject}_PWI_means_{H}.native.shape.gii'
#     group: 'map_perfusion_hpc'
#     singularity: config['singularity_connectomewb']   
#     shell:
#         "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"

# Transform ASL data to lores MP2RAGE for cortical analyses
rule warp_perf_to_anatomy:
    input:
        unpack(collect_input),
        nii = 'results/perfusion_clean/sub-{subject}/sub-{subject}_acq-ASL_run-all_moCorr_diCorr_perfusion_calib.nii.gz',
        xfm = 'results/perfusion_clean/sub-{subject}/sub-{subject}_acq-MZeroScan_distCorr_template_reg2Anat_greedy2ants_apply1st.txt'
    output: 'results/maps/sub-{subject}/sub-{subject}_CBF_space-MP2RAGE.nii.gz'
    group: 'map_perfusion_ctx'
    singularity: config['singularity_prepdwi']   
    shell:
        """
        ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}
        antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.nii} -t {input.xfm}
        """

rule sample_perf_cortex:
    input:
        nii = rules.warp_perf_to_anatomy.output,
        inner = 'results/hcp_mmp/sub-{subject}/{hemi}.white.32k_fs_LR.surf.gii',
        midthickness = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.32k_fs_LR.surf.gii',
        outer = 'results/hcp_mmp/sub-{subject}/{hemi}.pial.32k_fs_LR.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_hemi-{hemi}_CBF.32k_fs_LR.shape.gii'
    group: 'map_perfusion_ctx'
    singularity: config['singularity_connectomewb']
    shell:
        """           
        wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -trilinear        
        """

