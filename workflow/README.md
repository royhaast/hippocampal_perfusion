## Data processing workflow

This folder contains the data analysis pipeline built using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.

- Snakemake takes the `Snakefile` to automatically construct an analysis pipeline depending on the desired output file(s).

	- `DAG.pdf` contains the automatically constructed directed acyclic graph of the whole workflow
	- The `rules` folder contains sets of analysis steps organized per MRI image modality:

		- `rules/mp2rage.smk`: processing steps related to the high-resolution (0.5 mm isotropic) and standard resolution (1 mm isotropic) MP2RAGE data. Steps include ...
		- `rules/sa2rage.smk`:

To complete ...