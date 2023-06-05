# lookup tables for setting wb structure:
hemi_to_structure = {'L': 'CORTEX_LEFT', 'Lflip': 'CORTEX_LEFT', 'R': 'CORTEX_RIGHT'}
surf_to_secondary_type = {'midthickness': 'MIDTHICKNESS', 'inner': 'PIAL', 'outer': 'GRAY_WHITE'}

# Generate subject-specific warps to go from a template unfolded
# space to the subject's native space (and vice versa). Does this 
# for inner, midthickness and outer levels
def construct_warps_cmd (wildcards, input, output):
    singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

    set_autotop_dir = f"SINGULARITYENV_AUTOTOP_DIR=/project/6050199/rhaast/99_tmp_scratch_backup/HPC_hippunfold/hippocampal_autotop_ali_dev"
    set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
    set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

    cmd =   f"{set_autotop_dir} {set_matlab_lic} {set_java_home} {singularity_cmd} "\
            f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
            f"create_warps_giftis('{input.io_folder}','{output.io_folder}');\"" 
    return cmd  

rule create_unfold_to_native_warps:
    input:
        dseg = 'results/autotop-dev/results/sub-{subject}/anat/sub-{subject}_space-T1w_desc-subfields_modality-segT2w_template-CITI168_dseg.nii.gz',
        io_folder = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop',
    output:
        io_folder = directory('results/surface_warps/sub-{subject}/{H}'),
        flatmaps = expand('results/surface_warps/sub-{{subject}}/{{H}}/{surface}.unfoldedtemplate.surf.gii', surface=['inner','midthickness','outer']),
        native2unfold = 'results/surface_warps/sub-{subject}/{H}/Warp_native2unfold.nii',
        unfold2native = 'results/surface_warps/sub-{subject}/{H}/Warp_unfold2native_extrapolateNearest.nii',
        unfoldtemplate2unfold = 'results/surface_warps/sub-{subject}/{H}/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',
        unfold2unfoldtemplate = 'results/surface_warps/sub-{subject}/{H}/WarpITK_unfold2unfoldtemplate_0InverseWarp.nii.gz',        
    params:
        matlab_cmd = construct_warps_cmd
    group: 'surfaces'          
    shell:
        "SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "{params.matlab_cmd}"

rule convert_warp_to_world:
    input:
        unfoldtemplate2unfold = rules.create_unfold_to_native_warps.output.unfoldtemplate2unfold,
        unfold2unfoldtemplate = rules.create_unfold_to_native_warps.output.unfold2unfoldtemplate
    output:
        unfoldtemplate2unfold = 'results/surface_warps/sub-{subject}/{H}/Warp_unfoldtemplate2unfold.nii.gz',
        unfold2unfoldtemplate = 'results/surface_warps/sub-{subject}/{H}/Warp_unfold2unfoldtemplate.nii.gz'
    singularity: config['singularity_connectomewb']        
    group: 'surfaces'            
    shell:
        "wb_command -convert-warpfield -from-itk {input.unfoldtemplate2unfold} -to-world {output.unfoldtemplate2unfold} && "
        "wb_command -convert-warpfield -from-itk {input.unfold2unfoldtemplate} -to-world {output.unfold2unfoldtemplate}"

rule create_unfolded_from_template:
    input:
        flatmap = 'results/surface_warps/sub-{subject}/{H}/{surface}.unfoldedtemplate.surf.gii',
        unfoldtemplate2unfold = rules.convert_warp_to_world.output.unfoldtemplate2unfold,        
    output:
        flatmap = 'results/surface_warps/sub-{subject}/{H}/{surface}.unfolded.surf.gii',
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.H],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surface]                
    singularity: config['singularity_connectomewb']        
    group: 'surfaces'      
    shell:
        "wb_command -surface-apply-warpfield {input.flatmap} {input.unfoldtemplate2unfold} {output.flatmap} && "
        "wb_command -set-structure {output.flatmap} {params.structure_type} -surface-type FLAT -surface-secondary-type {params.secondary_type}"

def construct_correction_cmd (wildcards, input, output):
    singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

    set_autotop_dir = f"SINGULARITYENV_AUTOTOP_DIR=/project/6050199/rhaast/99_tmp_scratch_backup/HPC_hippunfold/hippocampal_autotop_ali_dev"
    set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
    set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

    cmd =   f"{set_autotop_dir} {set_matlab_lic} {set_java_home} {singularity_cmd} "\
            f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
            f"constrain_surf_to_warp('{input.flatmap}','{input.unfold2native}','{output}');\"" 
    return cmd  

rule correct_vertices_for_bbox:
    input:
        flatmap = 'results/surface_warps/sub-{subject}/{H}/{surface}.unfolded.surf.gii',
        unfold2native = 'results/surface_warps/sub-{subject}/{H}/Warp_unfold2native_extrapolateNearest.nii'
    output: 'results/surface_warps/sub-{subject}/{H}/{surface}_corrected.unfolded.surf.gii'  
    params:
        matlab_cmd = construct_correction_cmd
    group: 'surfaces'         
    shell:
        "SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "{params.matlab_cmd}" 

rule transform_unfolded_to_native:
    input:
        flatmap = rules.correct_vertices_for_bbox.output,
        unfold2native = rules.create_unfold_to_native_warps.output.unfold2native,        
    output:
        native = 'results/surface_warps/sub-{subject}/{H}/{surface}.native.surf.gii'
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.H],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surface]                
    singularity: config['singularity_connectomewb']        
    group: 'surfaces'       
    shell:
        "wb_command -surface-apply-warpfield {input.flatmap} {input.unfold2native} {output.native} && "
        "wb_command -set-structure {output.native} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type {params.secondary_type}"

rule downsample_surface:
    input: 'results/surface_warps/sub-{subject}/{H}/{surface}.native.surf.gii'
    output: 'results/surface_warps/sub-{subject}/{H}/{surface}.downsampled.native.surf.gii'
    singularity: config['singularity_freesurfer']
    shell:
        "mris_decimate -d 0.05 {input} {output}"

rule extract_surface_coordinates:
    input: 'results/surface_warps/sub-{subject}/{H}/{surface}.native.surf.gii'
    output: 'results/surface_warps/sub-{subject}/{H}/coords_{surface}.native.shape.gii'
    singularity: config['singularity_connectomewb']        
    group: 'surfaces'       
    shell:
        "wb_command -surface-coordinates-to-metric {input} {output}"

rule extract_surface_normals:
    input: 'results/surface_warps/sub-{subject}/{H}/{surface}.native.surf.gii'
    output: 'results/surface_warps/sub-{subject}/{H}/norms_{surface}.native.shape.gii'
    singularity: config['singularity_connectomewb']        
    group: 'surfaces'       
    shell:
        "wb_command -surface-normals {input} {output}"

rule surface_info_to_text:
    input:
        rules.extract_surface_coordinates.output,
        rules.extract_surface_normals.output
    output:
        'results/surface_warps/sub-{subject}/{H}/coords_{surface}.native.txt',
        'results/surface_warps/sub-{subject}/{H}/norms_{surface}.native.txt',
    group: 'surface'
    run:
        import numpy as np
        import nibabel as nib

        for in_metric, out_metric in zip(input,output):
            data = nib.load(in_metric)
            x = data.darrays[0].data
            y = data.darrays[1].data
            z = data.darrays[2].data
            xyz = np.vstack((x,y,z))

            np.savetxt(out_metric,np.transpose(xyz))

# # Compute average surfaces across all subjects, per hemisphere 
def construct_avg_surface_cmd (wildcards,input,output):
    cmd = "wb_command -surface-average {} -stddev {}".format(output.surface, output.stddev)
    for in_surface in input:
        cmd = "{0} -surf {1}".format(cmd, in_surface)
    return cmd 

rule create_average_surface:
    input: expand('results/surface_warps/sub-{subject}/{{H}}/{{surface}}.native.surf.gii',subject=subjects)
    output:
        surface = 'results/average_hippocampus/{surface}_hemi-{H}.native.surf.gii',
        stddev = 'results/average_hippocampus/{surface}_hemi-{H}_stdev.native.shape.gii'
    params:
        average_cmd = construct_avg_surface_cmd
        # structure_type = lambda wildcards: hemi_to_structure[wildcards.H],
        # secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surface]
    group: 'surfaces'     
    singularity: config['singularity_connectomewb']        
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "{params.average_cmd}"    



# rule apply_unfold_to_native_warps:
#     input:
#         autotops_folder = rules.create_unfold_to_native_warps.input.io_folder,
#         warps_folder = rules.create_unfold_to_native_warps.output.io_folder,
#         unfoldtemplate2unfold = rules.create_unfold_to_native_warps.unfoldtemplate2unfold, #'results/surface_warps/sub-{subject}/{H}/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',
#         unfold2native = rules.create_unfold_to_native_warps.unfold2native #'results/surface_warps/sub-{subject}/{H}/Warp_unfold2native_extrapolateNearest.nii',
#     output:
#         midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
#     params:
#         autotop_matlab = config['autotop_matlab']
#     singularity: config['singularity_autotop']        
#     group: 'surfaces'     
#     threads: 8
#     resources:
#         mem_mb = 32000         
#     shell:
#         "bash {params.autotop_matlab}/tools/warps_gifti/mapUnfoldToFullGrid_part1.sh {input.autotops_folder} {input.warps_folder}"
   
# rule convert_warp_to_world:
#     input:
#         unfoldtemplate2unfold = rules.create_unfold_to_native_warps.output.unfoldtemplate2unfold,
#         unfold2unfoldtemplate = rules.create_unfold_to_native_warps.output.unfold2unfoldtemplate
#     output:
#         unfoldtemplate2unfold = 'results/surface_warps/sub-{subject}/{H}/Warp_unfoldtemplate2unfold.nii.gz',
#         unfold2unfoldtemplate = 'results/surface_warps/sub-{subject}/{H}/Warp_unfold2unfoldtemplate.nii.gz'
#     singularity: config['singularity_connectomewb']        
#     group: 'surfaces'            
#     shell:
#         "wb_command -convert-warpfield -from-itk {input.unfoldtemplate2unfold} -to-world {output.unfoldtemplate2unfold} && "
#         "wb_command -convert-warpfield -from-itk {input.unfold2unfoldtemplate} -to-world {output.unfold2unfoldtemplate}"

# rule transform_unfolded_to_native:
#     input:
#         flatmap = 'results/surface_warps/sub-{subject}/{H}/{surface}.unfoldedtemplate.surf.gii',
#         unfoldtemplate2unfold = rules.convert_warp_to_world.output.unfoldtemplate2unfold,        
#         unfold2native = rules.create_unfold_to_native_warps.output.unfold2native,        
#     output:
#         flatmap = 'results/surface_warps/sub-{subject}/{H}/{surface}.unfolded.surf.gii',
#         native = 'results/surface_warps/sub-{subject}/{H}/{surface}.native.surf.gii'
#     params:
#         structure_type = lambda wildcards: hemi_to_structure[wildcards.H],
#         secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surface]                
#     singularity: config['singularity_connectomewb']        
#     group: 'surfaces'     
#     threads: 8
#     resources:
#         mem_mb = 32000    
#     shell:
#         "wb_command -surface-apply-warpfield {input.flatmap} {input.unfoldtemplate2unfold} {output.flatmap} && "
#         "wb_command -set-structure {output.flatmap} {params.structure_type} -surface-type FLAT -surface-secondary-type {params.secondary_type} && "
#         "wb_command -surface-apply-warpfield {output.flatmap} {input.unfold2native} {output.native} && "
#         "wb_command -set-structure {output.native} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type {params.secondary_type}"

# def construct_warp_cmd (wildcards, input):
#     singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

#     set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
#     set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

#     cmd = f"{set_matlab_lic} {set_java_home} {singularity_cmd} "\
#             f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
#             f"create_warps_old('{input.io_folder}','{input.io_folder}'); "\
#             f"extrapolate_warp_unfold2native('{input.io_folder}','{input.io_folder}');\"" 
#     return cmd   

# rule create_unfold_to_native_warps:
#     input:
#         dseg = rules.autotop.output,
#         io_folder = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop',
#     output:
#         native2unfold = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/Warp_native2unfold.nii',
#         unfold2native = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/Warp_unfold2native.nii'
#     params:
#         matlab_cmd = construct_warp_cmd
#     group: 'surfaces'        
#     threads: 8
#     resources:
#         mem_mb = 32000    
#     shell:
#         "SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "

# rule transform_unfolded_to_native:
#     input:
#         warp = rules.create_unfold_to_native_warps.output.unfold2native,
#         unfolded = expand('resources/{surface}.unfolded.surf.gii',surface=['inner','midthickness','outer'])
#     output:
#         native = expand('results/autotop-dev/work/autotop/sub-{{subject}}/sub-{{subject}}_hemi-{{H}}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/{surface}.nativecrop.surf.gii',surface=['inner','midthickness','outer'])
#     #singularity: config['singularity_connectomewb']        
#     group: 'surfaces'     
#     threads: 8
#     resources:
#         mem_mb = 32000    
#     run:
#         for in_unfolded, out_native in zip(input.unfolded, output.native):
#             shell("wb_command -surface-apply-warpfield {in_unfolded} {input.warp} {out_native}")
#             shell("wb_command -set-structure {out_native} CORTEX_RIGHT -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS")

# # This uses the original surfaces from Jordan's code and crops them
# # to the same resolution as the surfaces generated using the newer
# # (Ali's) code
# def construct_crop_cmd (wildcards, input):
#     singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

#     set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
#     set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

#     cmd = f"{set_matlab_lic} {set_java_home} {singularity_cmd} "\
#             f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
#             f"crop_SurfMap('{input.io_folder}','{input.io_folder}');\"" 
#     return cmd   

# rule crop_original_surfaces:
#     input:
#         dseg = rules.autotop.output,
#         io_folder = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop'
#     output:
#         expand('results/autotop-dev/work/autotop/sub-{{subject}}/sub-{{subject}}_hemi-{{H}}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-{depth}.surf.gii', depth=['1','2','3','4'])
#     params:
#         matlab_cmd = construct_crop_cmd
#     group: 'surfaces' 
#     threads: 8
#     resources:
#         mem_mb = 32000    
#     shell:
#         "{params.matlab_cmd}"

# # Generate midthickness surface using connectome workbench instead
# # of using the laplace coordinates
# rule generate_midthickness_surf:
#     input:
#         inner = rules.crop_original_surfaces.output[0], #'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-1.surf.gii',
#         outer = rules.crop_original_surfaces.output[3] #'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-4.surf.gii'
#     output: 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-mid.surf.gii',
#     params:
#         structure_type = lambda wildcards: hemi_to_structure[wildcards.H]
#     group: 'surfaces'     
#     singularity: config['singularity_connectomewb'] 
#     shell:
#         "wb_command -set-structure {input.inner} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type PIAL && "
#         "wb_command -set-structure {input.outer} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE && "
#         "wb_command -surface-cortex-layer {input.outer} {input.inner} 0.5 {output}  && "
#         "wb_command -set-structure {output} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS"
