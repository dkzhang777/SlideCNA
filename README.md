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
Prepare a Seurat Data object of the Slide-seq data that contains counts matrix and metadata with cell type annotations. Metadata should contain the following columns in the provided format:

bc (chr): bead labels \
seurat_clusters (fct): Seurat-defined clusters<sup>*</sup>\
pos_x (dbl): x-coordinate bead position\
pos_y (dbl): y-coordinate bead position\
cluster_type (chr): annotation of the bead as 'Normal' (Non-malignant) or 'Malignant'
    
<sup>*</sup> Seurat clusters are included to guide interpretation of cell type annotations. RCTD-defined cell types and malignant/non-malignant status should generally correspond to distinct Seurat clusters.

### Running SlideCNA
```
run_slide_cna(so, 
              md, 
              gene_pos,
              plotDir,
              OUPUT_DIRECTORY,
              spatial=TRUE,
              roll_mean_window=101,
              avg_bead_per_bin=12, 
              pos=TRUE, pos_k=55, 
              ex_k=1, 
              hc_function_bin='ward.D2', 
              spatial_vars_to_plot=c("seurat_clusters", "bin_all", "N_bin", 
                                     "umi_bin", "cluster_type"),
              scale_bin_thresh_hard=TRUE, 
              lower_bound_cnv=0.6, 
              upper_bound_cnv=1.4, 
              hc_function_cnv='ward.D2', 
              hc_function_cnv_heatmap = 'ward.D2',
              quantile_plot_cluster_label="seurat_clusters", 
              hc_function_silhouette='ward.D2',
              max_k_silhouette=10, 
              plot_silhouette=TRUE, 
              hc_function_plot_clones="ward.D2", 
              chrom_ord = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 
                            'chr8', 'chr9','chr10', 'chr11', 'chr12', 'chr13', 
                            'chr14', 'chr15', 'chr16', 'chr17','chr18', 'chr19', 
                            'chr20', 'chr21', 'chr22', 'chr23', 'chrX', 'chrY', 'chrM'),
              # Heat map legend chromosome colors
              chrom_colors = c(chr1='#8DD3C7', chr2='#FFFFB3', chr3='#BEBADA', 
                               chr4='#FB8072', chr5='#80B1D3', chr6='#FDB462', 
                               chr7='#B3DE69', chr8='#FCCDE5', chr9='#D9D9D9' , 
                               chr10='#BC80BD', chr11='#CCEBC5', chr12='#FFED6F', 
                               chr13='#1B9E77', chr14='#D95F02', chr15='#7570B3', 
                               chr16='#E7298A', chr17='#66A61E', chr18='#E6AB02', 
                               chr19='#A6761D', chr20='#666666', chr21='#A6CEE3', 
                               chr22='#1F78B4', chrX='#B2DF8A'),
              # Ggplot parameters
              text_size = 16,
              title_size = 18,
              legend_size_pt = 4,
              legend_height_bar = 1.5)
```

### Parameter Descriptions

`so` (Seurat Object): object that contains a counts matrix and metadata with cell type annotations \
`md` (data.table): metadata of each bead (beads x annotations); contains columns 'bc' with values of bead names, 'cluster_type' with values of 'Normal' or 'Malignant', 'pos_x' with x-coordinate bead positions, 'pos_y' with y-coordinate bead positions, and 'nCount_RNA' with RNA counts per bead \
`gene_pos` (data.table): table with columns for GENE, chr, start, end, rel_gene_pos (1 : # of genes on chromosome)\
`plotDir` (char): output plot directory path\
`OUTPUT_DIRECTORY` (char): output directory path\
`spatial` (bool): TRUE if using spatial information FALSE if not\
`roll_mean_window` (int): integer number of adjacent genes for which to average over in pyramidal weighting scheme\
`avg_bead_per_bin` (int): integer of average number of beads there should be per bin \
`pos` (bool): TRUE if doing spatial and expressional binning, FALSE if just expressional binning\
`pos_k` (numeric): positional weight\
`ex_k` (numeric): expressional weight\
`hc_function_bin` (char): hierarchical clustering function for binning. Choose from "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)\
`spatial_vars_to_plot` (char): vector of features to plot/columns of metadata\
`scale_bin_thresh_hard` (bool): TRUE if using strict, absolute thresholds for expression thresholds and FALSE if adjusting depending on the sample\
`lower_bound_cnv` (numeric): lower threshold for CNV scores\
`upper_bound_cnv` (numeric): upper threshold for CNV scores \
`hc_function_cnv` (char): hierarchical clustering function to use for CNV-calling. Choose from "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)\
`hc_function` (char): hierarchical clustering function to use for visualzing CNV heat map. This should be consistent with `hc_function_cnv`. Choose from "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)\
`quantile_plot_cluster_label` (char): column name to keep in quantile plot\
`hc_function_silhouette` (char): hierarchical clustering function to use for the Silhouette method. Choose from "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)\
`max_k_silhouette` (int): max number of clusters to evaluate (2 to max_k_silhouette) in Silhouette method\
`hc_function_plot_clones` (char): hierarchical clustering function to use in plotting clones. Choose from "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)\
`chrom_ord` (char):  vector of order and names of chromosomes\
`chrom_colors` (char): vector of which colors each chromosome should be in heat map\
`text_size` (int): size of text in some ggplots\
`title_size` (int): size of title in some ggplots\
`legend_size_pt` (int): size of legend text size in some ggplots\
`legend_height_bar` (int): height of legend bar in some ggplots

### Results
Results will appear in OUTPUT_DIRECTORY and plotDir
