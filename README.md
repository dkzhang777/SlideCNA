# SlideCNA 

### Introduction
SlideCNA is a method to call copy number alterations (CNA) from spatial transcriptomics data (adapted for Slide-seq data). SlideCNA uses expression smoothing across the genome to extract changes in copy number and implements a spatio-molecular binning process to boost signal and consolidate reads. Based on the CNA profiles, SlideCNA can identify clusters across space. 

Example Jupyter notebooks of SlideCNA applied to Slide-seq, snRNA-seq, and Slide-seq with TACCO bead splitting data are available here: https://github.com/dkzhang777/SlideCNA_Analysis.

### Installation

#### New Conda Environment

Create a new conda environment using the SlideCNA_env.yml file from the SlideCNA repository:
```
conda env create -f "https://github.com/dkzhang777/SlideCNA/blob/main/SlideCNA_env.yml"
```

Install SlideCNA through R from Github:
```
library(devtools)
devtools::install_github("dkzhang777/SlideCNA@main", force=TRUE)
library(SlideCNA)
```
### Preparation
Preparation of Slide-seq data raw counts matrix and meta data with cell type annotations. Metadata should contain the following columns in the provided format:

bc (chr): bead labels \
cluster_type (chr): annotation of the bead as 'Normal' (Non-malignant) or 'Malignant' 

and, if using spatially-aware binning:

pos_x (dbl): x-coordinate bead position\
pos_y (dbl): y-coordinate bead position
    
### Running SlideCNA
```
run_slide_cna(counts, 
              beads_df, 
              gene_pos, 
              output_directory, 
              plot_directory,
              spatial=TRUE,
              roll_mean_window=101,
              avg_bead_per_bin=12,
              pos=TRUE, 
              pos_k=55, 
              ex_k=1)
```

### Parameter Descriptions

`counts` (data.frame): raw counts (genes x beads) \
`beads_df` (data.frame): annotations of each bead (beads x annotations); contains columns 'bc' for bead names, 'cluster_type' for annotations of 'Normal' or 'Malignant', 'pos_x' for x-coordinate bead positions, and 'pos_y' for y-coordinate bead positions \
`gene_pos` (data.frame): table with columns for GENE, chr, start, end, rel_gene_pos (1 : # of genes on chromosome)\
`output_directory` (char): output directory path\
`plot_directory` (char): output plot directory path\
`spatial` (bool): TRUE if using spatial information FALSE if not\
`roll_mean_window `(int): integer number of adjacent genes for which to average over in pyramidal weighting scheme\
`avg_bead_per_bin` (int): integer of average number of beads there should be per bin\
`pos` (bool): TRUE if doing spatial and expressional binning, FALSE if just expressional binning\
`pos_k` (numeric): positional weight\
`ex_k` (numeric): expressional weight

### Results

Results will appear in output_directory and plot_directory. Key output files are described below: \ 
`so.rds` Seurat object of Slide-seq data \
`md.txt` metadata of Slide-seq data with Seurat annotations \
`md_bin.txt` metadata of binned Slide-seq data \
`dat_bin_scaled.txt` expression values of binned Slide-seq data after applying pyramidal weighting scheme and normalizing for UMI per bin\
`best_k_malig.rds` value of optimal number of malignant clusters \
`cluster_labels_all.txt` cluster assignments when performing cluster designation on all binned beads \
`cluster_labels_malig.txt` cluster assignments when performing cluster determination on only malignant binned beads \
`cluster_markers_all.txt` DEGs per cluster when performing cluster designation on all binned beads \
`cluster_markers_malig.txt` DEGs per cluster when performing cluster determination on only malignant binned beads \
`go_terms_all.txt` GO terms per cluster when performing cluster designation on all binned beads \
`go_terms_malig.txt` GO terms per cluster when performing cluster determination on only malignant binned beads \
