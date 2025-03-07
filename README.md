[![DOI](https://zenodo.org/badge/649706083.svg)](https://zenodo.org/doi/10.5281/zenodo.10644516)
## Hippocampal perfusion

This repository contains the data analysis code and (visualization) notebooks related to the paper:

> "Insights into hippocampal perfusion using high-resolution, multi-modal 7T MRI"

that can be found [here](https://doi.org/10.1073/pnas.2310044121).

- The analysis pipeline (see `workflow` for more details) is build using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.
- Several `notebooks` are uploaded to support the inspection of the intermediate analysis steps (i.e., hippocampal segmentation and vasculature), as well as the different visualization steps (in `visualization`).
- `results/surface_maps` contains group-averaged surface maps as displayed in the manuscript figures. 
- A more interactive visualization tool of the results is possible using our online data explorer that can be accessed [here](https://tinyurl.com/3z8czuy9/).

![Hippocampal Data Explorer](HippocampalDataExplorer.gif?raw=true "Hippocampal Data Explorer")

#### Dependencies
- Singularity and containers (see `config/config.yml`)
- A virtual environment is recommended with the following packages installed:
    - Python (3.0 or higher)
    - Snakemake 
    - [HippUnfold](https://hippunfold.readthedocs.io/)

**Note**
Raw imaging data, manual segmentations as well as the related code for the image preprocessing steps will be uploaded to a separate data repository after the data manuscript has been accepted for publication.
