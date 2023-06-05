rule freesurfer_to_gifti:
    input: 'results/freesurfer/sub-{subject}/surf/{hemi}.{surfname}'
    output: 'results/hcp_mmp/sub-{subject}/{hemi}.{surfname}.fs.surf.gii'
    group: 'hcp_mmp'
    singularity: config['singularity_freesurfer']        
    shell: "mris_convert {input} {output}"

rule freesurfer_to_nifti:
    input: 'results/freesurfer/sub-{subject}/mri/T1.mgz'
    output: 'results/hcp_mmp/sub-{subject}/T1.nii.gz'
    group: 'hcp_mmp'
    singularity: config['singularity_freesurfer']        
    shell: "mri_convert {input} {output}"

rule get_tkr2scanner:
    input: rules.freesurfer_to_nifti.output
    output: 'results/hcp_mmp/sub-{subject}/tkr2scanner.xfm'
    group: 'hcp_mmp'
    singularity: config['singularity_freesurfer']
    shell: 'mri_info {input} --tkr2scanner > {output}'
     
rule apply_surf_tkr2scanner:
    input: 
        surf = rules.freesurfer_to_gifti.output,
        tkr2scanner = rules.get_tkr2scanner.output
    output: 'results/hcp_mmp/sub-{subject}/{hemi}.{surfname}.native.surf.gii'
    group: 'hcp_mmp'
    singularity: config['singularity_connectomewb']
    shell: 'wb_command -surface-apply-affine {input.surf} {input.tkr2scanner} {output}'

rule gen_midthickness:
    input:
        white = 'results/hcp_mmp/sub-{subject}/{hemi}.white.native.surf.gii',
        pial = 'results/hcp_mmp/sub-{subject}/{hemi}.pial.native.surf.gii'
    output: 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.native.surf.gii'
    group: 'hcp_mmp'
    singularity: config['singularity_connectomewb']
    shell: 'wb_command -surface-cortex-layer {input.white} {input.pial} 0.5 {output}'
   
rule resample_subj_to_fsaverage_sphere:
    input: 
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}.{surfname}.native.surf.gii',
        current_sphere = 'results/hcp_mmp/sub-{subject}/{hemi}.sphere.reg.fs.surf.gii',
        new_sphere = 'resources/standard_mesh_atlases/resample_fsaverage/'
                        'fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii'
    params:
        method = 'BARYCENTRIC'
    output: 'results/hcp_mmp/sub-{subject}/{hemi}.{surfname}.32k_fs_LR.surf.gii'
    group: 'hcp_mmp'
    singularity: config['singularity_connectomewb']
    shell: 'wb_command -surface-resample {input.surf} {input.current_sphere} {input.new_sphere} {params.method} {output}'

rule gen_inflated:
    input: 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.native.surf.gii'
    output:
        inflated = 'results/hcp_mmp/sub-{subject}/{hemi}.inflated.native.surf.gii',
        very_inflated = 'results/hcp_mmp/sub-{subject}/{hemi}.very_inflated.native.surf.gii'
    group: 'hcp_mmp'
    singularity: config['singularity_connectomewb']      
    shell:
        "wb_command -surface-generate-inflated {input} {output.inflated} {output.very_inflated} -iterations-scale 0.75"

rule resample_labels_to_subj_sphere:
    input:
        label = 'resources/standard_mesh_atlases/{hemi}.hcp-mmp.32k_fs_LR.label.gii',
        current_sphere = 'resources/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii',
        new_sphere = 'results/hcp_mmp/sub-{subject}/{hemi}.sphere.reg.fs.surf.gii',
        current_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.32k_fs_LR.surf.gii',
        new_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.native.surf.gii'
    params:
        method = 'ADAP_BARY_AREA'
    output: 'results/hcp_mmp/sub-{subject}/{hemi}.hcp-mmp.native.label.gii'
    group: 'hcp_mmp'
    singularity: config['singularity_connectomewb']
    shell: 
        'wb_command -label-resample {input.label} {input.current_sphere} {input.new_sphere}'
        ' {params.method} {output} -area-surfs {input.current_surf} {input.new_surf}'

rule map_labels_to_volume_ribbon:
    input: 
        label = rules.resample_labels_to_subj_sphere.output,
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.native.surf.gii',
        vol_ref = rules.freesurfer_to_nifti.output,
        white_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.white.native.surf.gii',
        pial_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.pial.native.surf.gii',
    output: 'results/hcp_mmp/sub-{subject}/{hemi}.hcp-mmp.nii.gz'
    group: 'hcp_mmp'
    singularity: config['singularity_connectomewb']
    shell:
        'wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output}'
        ' -ribbon-constrained {input.white_surf} {input.pial_surf} -greedy'