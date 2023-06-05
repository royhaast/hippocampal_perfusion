# Common rules to merge and reduce (avg, std and cov) surface
# maps across subjects and hemispheres. Averages across subjects
# first
rule extract_gm_ribbon:
    input: 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/coords-IO.nii.gz'
    output: 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/ribbon.nii.gz'
    #group: 'map_general'    
    singularity: config['singularity_prepdwi']
    shell:
        "fslmaths {input} -bin {output}"

def construct_merge_cmd (wildcards,input,output):
    cmd = "wb_command -metric-merge {}".format(output)
    input = input[1:] if 'vol0000' in input[0] else input
    column = ""
    if ("parameter" in wildcards) or ("asl_parameter" in wildcards):
        column = "-column 1"
    for in_metric in input:
        cmd = "{0} -metric {1} {2}".format(cmd, in_metric, column)
    return cmd

rule merge_group_metric_gii:
    input: expand('results/surface_maps/sub-{subject}/sub-{subject}_{{tpl_res}}{{parameter}}_{{H}}.native.shape.gii', subject=subjects)
    output: 'results/surface_maps/group/sub-group_{{tpl_res}}{parameter}_{H}.native.shape.gii'
    params:
        merge_cmd = construct_merge_cmd
    group: 'map_general'    
    singularity: config['singularity_connectomewb']  
    shell:
        "{params.merge_cmd}"    

rule reduce_group_metric_gii:
    input: rules.merge_group_metric_gii.output 
    output:
        avg = 'results/surface_maps/group/sub-group_{tpl_res}{parameter}_{H}_avg.native.shape.gii',
        std = 'results/surface_maps/group/sub-group_{tpl_res}{parameter}_{H}_std.native.shape.gii',
        cov = 'results/surface_maps/group/sub-group_{tpl_res}{parameter}_{H}_cov.native.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'map_general'      
    shell:
        "wb_command -metric-reduce {input} MEAN {output.avg} -only-numeric &&"
        "wb_command -metric-reduce {input} STDEV {output.std} -only-numeric &&"
        "wb_command -metric-reduce {input} COV {output.cov} -only-numeric "   

# Average across subjects AND hemispheres
rule reduce_group_lr_metric_gii:
    input: expand('results/surface_maps/sub-{subject}/sub-{subject}_{{parameter}}_{H}.native.shape.gii', subject=subjects, H=['R','Lflip'])
    output:
        merged = 'results/surface_maps/group/sub-group_{parameter}_LR.native.shape.gii',
        avg = 'results/surface_maps/group/sub-group_{parameter}_LR_avg.native.shape.gii',
        std = 'results/surface_maps/group/sub-group_{parameter}_LR_std.native.shape.gii',
        cov = 'results/surface_maps/group/sub-group_{parameter}_LR_cov.native.shape.gii'
    group: 'map_general'        
    run:
        import numpy as np
        import nibabel as nib
        from scipy.stats import variation 

        metric = nib.load(input[0])
        data = np.zeros((len(metric.darrays[0].data),len(input)))

        gii_merged = nib.gifti.GiftiImage()

        for i,in_metric in enumerate(input):
            data[:,i] = nib.load(in_metric).darrays[0].data
            gii_merged.add_gifti_data_array(
                nib.gifti.GiftiDataArray(data=data[:,i].astype(np.float32)
                )
            )       
        nib.save(gii_merged, output.merged)
        shell("wb_command -set-structure {output.merged} CORTEX_RIGHT")

        metric.darrays[0].data = np.nanmean(data,1)
        nib.save(metric, output.avg)

        metric.darrays[0].data = np.nanstd(data,1)
        nib.save(metric, output.std)

        metric.darrays[0].data = np.nanstd(data,1)/np.nanmean(data,1)
        nib.save(metric, output.cov)

# rule smooth_group_lr_metric_gii:
#     input:
#         metric = 'results/surface_maps/group/sub-group_{parameter}_LR_avg.native.shape.gii',
#         #midthickness = 'results/average_hippocampus/midthickness_hemi-R.nativecrop.surf.gii'
#     output: 'results/surface_maps/group/sub-group_{parameter}_LR_avg_smoothed.native.shape.gii'
#     params:
#         sigma = 0.6
#     group: 'map_general'    
#     singularity: config['singularity_connectomewb']  
#     shell:
#         "wb_command -metric-smoothing {input.metric} {params.sigma} {output} -fix-zeros"       

rule calculate_gradients:
    input:
        metric = 'results/surface_maps/group/sub-group_{parameter}_LR_avg.native.shape.gii',
        midthickness = 'results/average_hippocampus/midthickness_hemi-R.native.surf.gii'
    output: 'results/surface_maps/group/sub-group_{parameter}_LR_avg_gradient.native.shape.gii'
    params:
        sigma = 1    
    group: 'map_general'    
    singularity: config['singularity_connectomewb']  
    shell:
        "wb_command -metric-gradient {input.midthickness} {input.metric} {output} -presmooth {params.sigma} -average-normals"
