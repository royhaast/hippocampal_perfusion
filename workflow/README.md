## Data processing workflow

This folder contains the data analysis pipeline built using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.

- Snakemake takes the `Snakefile` to automatically construct an analysis pipeline depending on the desired output file(s).
    - `DAG.pdf` contains the automatically constructed directed acyclic graph of the whole workflow
    - The `rules` folder contains sets of analysis steps organized per MRI image modality:
        - `rules/common.smk`: common rules that are called at different steps throughout the pipeline.
        - `rules/mp2rage.smk`: processing steps related to the high-resolution (0.5 mm isotropic) and standard resolution (1 mm isotropic) MP2RAGE data. Steps include removal of non-brain tissue, FreeSurfer-based brain segmentation, parcellation and surface reconstructions. 
        - `rules/hcp-mmp.smk`: converts FreeSurfer output to GiFTI format and projects the HCP-MMP atlas onto the subject's native as well as resampled (i.e., to `fs_LR_32k` vertices) surfaces.
        - `rules/sa2rage.smk`: samples the B1+ map onto the hippocampal surfaces.
        - `rules/autotop.smk`: prepares the manual segmentations, T2w volumes and then runs the actual hippocampal unfolding using the HippUnfold legacy code. 
        - `rules/surfaces.smk`: does some postprocessing of the hippocampal surfaces to correct vertex positioning, calculate native to unfolded warps and averaging across subjects.
        - `rules/morphology.smk`: calculates several morphological metrics using the hippocampal surfaces and generates notebooks with and overview of the T2w data, manual segmentations and HippUnfold output.
        - `rules/t2w.smk`: basic steps to project T2w intensity onto hippocampal surfaces.
        - `rules/tof.smk`: performs TOF data filtering and preparation for semi-automatic vascular segmentation using MeViSlab. Generates notebooks and calculates several metrics for analyses of hippocampal vasculature. 
        - `rules/perfusion.smk`: last but not least, full fit- and run-wise preprocessed ASL data (e.g., perfusion, tSNR, ... maps) are finally projected onto the hippocampal surfaces. 

To complete ...