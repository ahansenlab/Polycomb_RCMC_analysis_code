# Polycomb RCMC Analysis Code
This repository contains source code for the preprint [Dissection of Polycomb complexes reveals canonical PRC1-independent 3D chromatin interactions](https://www.biorxiv.org/content/10.64898/2026.10.13.758654v1), used in the analysis of data presented in the manuscript. 

Code is provided either in the form of Python scripts or Jupyter notebooks to be run in conda environments containing the required packages, or as standalone R scripts requiring installation of the packages listed in the code.

## Code summary
### bigWig generation (makebigwig_mac3.py)
Required packages:
-	MACS3
-	samtools
-	sambamba
-	bedtools
-	wigToBigWig

Python script used to generate bigWig format files from bam files for both ChIP-seq and TTchem-seq experiments.

Example usage:

```
python /path/to/script/makebigwig_macs3.py list_of_bams.tsv paired mm39 --outdir exampleoutdir
```

### TTchem-seq alignment (makebigwig_mac3.py)
Required packages:
-	Bowtie2
-	samtools
-	sambamba
-	STAR

Python script used to aligned paired-end reads from TTchem-seq data with spike-in (will also work for other RNA-seq-type datasets).

Example usage:

```
python /path/to/script/RNA_PE_alignment.py -f list_to_align.tsv -g mm39 -s dm6 -t 12 -o RNAalignmentcounts --outdir alignments
```

### Calculating individual loop strengths (Interaction_Strength.ipynb)
Required packages:
-	seaborn
-	numpy
-	matplotlib
-	pandas
-	coolpuppy
-	cooler
-	bioframe
-	cooltools

Jupyter notebook for running on RCMC mcools to calculate the strengths of individual loops used for filtering.

### Filtering interactions based on individual strengths (interaction_filter.R)
Required packages:
-	ggsignif
-	tidyr
-	plyr
-	dplyr
-	ggplot2
-	reshape2
-	ggExtra
-	purrr
-	stringr
-	InteractionSet

R script used to filter pairwise interactions between Polycomb sites after calculation individual interaction strengths using the above Jupyter notebook.

### Generating pileups of interactions for aggregate analysis (PcG_Pileups.ipynb)
Required packages:
-	seaborn
-	numpy
-	matplotlib
-	pandas
-	coolpuppy
-	cooler
-	bioframe
-	cooltools

Jupyter notebook for running on RCMC mcools to generate aggregate pileup analyses, both for 2D and 1D plots.

### Generating plots of genomics data (PcG_Region_Plotting_Figs.ipynb)
Required packages:
-	seaborn
-	numpy
-	matplotlib
-	pandas
-	coolpuppy
-	cooler
-	bioframe
-	cooltools
-	coolbox
-	pygenometracks

Jupyter notebook for generating snapshot plots of genomics data at given regions for figures.

### Analysing Polycomb bodies (PcG_bodies_analysis.R)
Required packages:
-	readtext
-	ggplot2
-	ggsignif
-	plyr
-	dplyr
-	reshape2
-	stringr
-	ggbeeswarm

R script used to process outputs from Polycomb body segmentation and produce plots for figures.

### Plotting distance distributions for Polycomb locus imaging (Plotting_distances.R)
Required packages:
-	readtext
-	ggplot2
-	ggsignif
-	plyr
-	dplyr
-	reshape2
-	stringr
-	ggbeeswarm

R script used to process outputs from single locus microscopy 3D distance measurements and produce plots for figures.

### Using the modified [Fyrtarn](https://github.com/ahansenlab/Fyrtarn/tree/main) pipeline for snapshots (fyrtarn_for_snapshots/)
To use, follow the below instructions:
1. Create a conda environment
```
conda create --name fyrtarn_for_snapshots python=3.13
conda activate fyrtarn_for_snapshots
```
2. Create a folder to clone the repo to
```
cd directory/for/repo/download
git clone https://github.com/ahansenlab/Polycomb_RCMC_analysis_code
cd Polycomb_RCMC_analysis_code/fyrtarn_for_snapshot
```
3. Follow the below instructions to install necessary packages. Note from the original Fyrtarn: Support for (NVIDIA) GPUs is enabled through CuPy. In order to work properly, this dependency must be configured to the CUDA Toolkit version. For example, an environment with CUDA 12.9 would use the following (i.e. pip install -e .[12])
```
mkdir .install
pip install -e .[CUDA]
pip install -e .[12]
```

## How to cite
This work is shared under an MIT license. If you make use of analysis scripts or data from this work, please cite as follows:

Huseyin, M. K., Kim, S., Toppen, J., Blackledge, N. P. & Hansen, A. S. Dissection of Polycomb complexes reveals canonical PRC1-independent 3D chromatin interactions. 2026.10.13.758654 Preprint at https://doi.org/10.64898/2026.10.13.758654 (2026).
