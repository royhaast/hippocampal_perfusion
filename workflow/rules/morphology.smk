# Gyrification
rule calculate_vertex_area_from_native_surface:
    input: 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_vertexarea_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb'] 
    group: 'morphology' 
    shell:
        "wb_command -surface-vertex-areas {input} {output}"    

rule calculate_vertex_area_from_unfolded_surface:
    input: 'results/surface_warps/sub-{subject}/{H}/midthickness.unfoldedtemplate.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_vertexarea_{H}.unfolded.shape.gii'
    singularity: config['singularity_connectomewb'] 
    group: 'morphology' 
    shell:
        "wb_command -surface-vertex-areas {input} {output}"    

rule calculate_gyrification_from_area:
    input:
        native = 'results/surface_maps/sub-{subject}/sub-{subject}_vertexarea_{H}.native.shape.gii',
        unfolded = 'results/surface_maps/sub-{subject}/sub-{subject}_vertexarea_{H}.unfolded.shape.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_gyrification_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb'] 
    group: 'morphology' 
    shell:    
        """
        wb_command -metric-math x/y {output} -var x {input.native} -var y {input.unfolded}
        """

# Curvature
rule smooth_midthickness:
    input: 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
    output: 'results/surface_warps/sub-{subject}/{H}/midthickness.smoothed.native.surf.gii'
    params:
        strength=0.6,
        iterations=100
    singularity: config['singularity_connectomewb']
    group: 'morphology'       
    shell:
        "wb_command -surface-smoothing {input} {params.strength} {params.iterations} {output}"

rule calculate_curvature_from_surface:
    input: 'results/surface_warps/sub-{subject}/{H}/midthickness.smoothed.native.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_curvature_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'       
    shell:
        "wb_command -surface-curvature {input} -mean {output}"     

# Thickness
rule calculate_thickness_from_surface:
    input:
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'       
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"  

# Calculate distance from subject-specific surfaces to group
# template. Does this for midthickness only.
rule calculate_distance_to_avg:
    input:
        avg = 'results/average_hippocampus/midthickness_hemi-{H}.native.surf.gii',
        subject = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
    output:
        distance = 'results/surface_maps/sub-{subject}/sub-{subject}_surfdist_{H}.native.shape.gii',
        displacement = 'results/surface_maps/sub-{subject}/sub-{subject}_surfdisp_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'              
    shell:    
        "wb_command -surface-to-surface-3d-distance {input.avg} {input.subject} {output.distance} -vectors {output.displacement} && "
        "wb_command -metric-math 'metric' {output.distance} -fixnan 0 -var metric {output.distance} && "
        "wb_command -metric-math 'metric' {output.displacement} -fixnan 0 -var metric {output.displacement}"
        
        
# Generate segmentation result notebooks
rule segmentation_notebooks:
    input:
        seg = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/manual_lbl.nii.gz',
        labels = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/subfields-BigBrain.nii.gz',
        t2w = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/img.nii.gz',
        thickness = 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.native.shape.gii',
        surface = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
    output: 'visualization/segmentation/sub-{subject}/sub-{subject}_hemi-{H}_thickness_folded.png'
    log:
        notebook="logs/notebooks/sub-{subject}/sub-{subject}_hemi-{H}_segmentation.ipynb"
    notebook:
        "../notebooks/segmentation.py.ipynb"
        
rule segmentation_notebooks_html:
    input: 'logs/notebooks/sub-{subject}/sub-{subject}_hemi-{H}_segmentation.ipynb'
    output: 'logs/notebooks/sub-{subject}/sub-{subject}_hemi-{H}_segmentation.html'
    shell:
        """
        jupyter nbconvert \
            --TagRemovePreprocessor.enabled=True \
            --TagRemovePreprocessor.remove_cell_tags snakemake-job-properties \
            --to html {input}
        """