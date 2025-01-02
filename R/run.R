#' Run SlideCNA workflow 
#'
#' Take a raw expression counts, cell type annotations, and positional cooridnates to identify CNA patterns across space and CNA-based clustering patterns
#'
#' @param counts data.frame of raw counts (genes x beads)
#' @param beads_df data.frame of annotation of each bead (beads x annotations); contains columns 'bc' for bead names, 'cluster_type' for annotations of 'Normal' or 'Malignant', 'pos_x' for x-coordinate bead positions, and 'pos_y' for y-coordinate bead positions
#' @param gene_pos data.frame with columns for GENE, chr, start, end, rel_gene_pos (1 : # of genes on chromosome)
#' @param plot_directory output plot directory path
#' @param output_directory output directory path
#' @param spatial TRUE if using spatial information FALSE if not
#' @param roll_mean_window integer number of adjacent genes for which to average over in pyramidal weighting scheme
#' @param avg_bead_per_bin integer of average number of beads there should be per bin 
#' @param pos TRUE if doing spatial and expressional binning, FALSE if just expressional binning
#' @param pos_k positional weight
#' @param ex_k expressional weight
#' @param hc_function_bin hierarchical clustering function for binning; to feed hclust's method argument, one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
#' @param spatial_vars_to_plot character vector of features to plot/columns of metadata
#' @param scale_bin_thresh_hard TRUE if using strict thresholds for expression thresholds and FALSE if adjusting 
#' thresholds based on 1 + or - the mean of absolute min and max vlaues
#' @param lower_bound_cnv numeric float to represent the lower cap for CNV scores
#' @param upper_bound_cnv numeric float to represent the upper cap for CNV scores 
#' @param hc_function_cnv character for which hierarchical clustering function to use for CNV-calling; to feed hclust's method argument, one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
#' @param hc_function_cnv_heatmap character for which hierarchical clustering function to use for visualzing CNV heat map; to feed hclust's method argument, one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
#' @param quantile_plot_cluster_label character string of which column name to keep in quantile plot
#' @param hc_function_silhouette character string for which hierarchical clustering function to use for 
#'        the Silhouette method; to feed hclust's method argument, one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
#' @param max_k_silhouette integer of number max number of clusters to evaluate (2:max_k_silhouette) 
#'.       in Silhouette method
#' @param plot_silhouette TRUE if plotting silhouette scores for clustering 
#' @param hc_function_plot_clones character string for which hierarchical clustering function to use 
#'        in plotting clones
#' @param use_GO_terms TRUE if using enrichR to get Gene Ontology terms for SlideCNA-defined clusters
#' @param chrom_ord character vector of order and names of chromosomes
#' @param chrom_colors character vector of which colors each chromosome should be in heat map
#' @param text_size integer of size of text in some ggplots
#' @param title_size integer of size of title in some ggplots
#' @param legend_size_pt integer of size of legend text size in some ggplots
#' @param legend_height_bar integer of height of legend bar in some ggplots
#' 
#' @export

run_slide_cna <- function(counts, 
                          beads_df, 
                          gene_pos, 
                          output_directory, 
                          plot_directory,
                          spatial=TRUE,
                          roll_mean_window=101,
                          avg_bead_per_bin=12, 
                          pos=TRUE, 
                          pos_k=55, 
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
                          use_GO_terms=TRUE,
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
                          legend_height_bar = 1.5) {
    
    setwd(output_directory)
    
    # Set up logger
    log_file <- paste0(output_directory, "/SlideCNA.log")
    if (file.exists(log_file)) {file.remove(log_file)}
    futile.logger::flog.appender(futile.logger::appender.file(log_file))
    
    # Set a global error handler
    options(error = function(e) {
        futile.logger::flog.error(paste("Global error handler: Error encountered:", conditionMessage(e)))
    })
    
    # Check if enrichR is installed if obtaining GO terms for SlideCNA-based clusters
    if (isTRUE(use_GO_terms)) {
      if (!requireNamespace("packageName", quietly = TRUE)) {
          stop("The 'enrichR' package is required if use_GO_terms is true.",
               "Please install it using install_github('wjawaid/enrichR').")
      }
    }

    # Make initial Seurat object and meta data object
    
    # Reformat counts and create a sparse matrix for input into Seurat
    counts_mat <- counts %>% 
              dplyr::select(beads_df$bc) %>%
              data.table::as.data.table() %>%
              mltools::sparsify()
    row.names(counts_mat) <- row.names(counts)
    
    # Ensure beads_df has correct columns 
    if (isTRUE(spatial)) {
        if (!("bc" %in% colnames(beads_df)) | 
            !("cluster_type" %in% colnames(beads_df)) | 
            !("pos_x" %in% colnames(beads_df)) | 
            !("pos_y" %in% colnames(beads_df))) {
            error_msg <- "beads_df does not contain all necessary columns (bc, cluster_type, pos_x, pos_y)"
            futile.logger::flog.error(error_msg)
            stop(error_msg)
        }
    }
    else {
        if (!("bc" %in% colnames(beads_df)) | 
            !("cluster_type" %in% colnames(beads_df))) {
            error_msg <- "beads_df does not contain all necessary columns (bc, cluster_type"
            futile.logger::flog.error(error_msg)
            stop(error_msg)
        }
    }
        
    # Make Seurat object
    futile.logger::flog.info("making initial Seurat object")
    
    so <- SlideCNA::make_seurat_annot(counts_mat, 
                                      beads_df, 
                                      seed_FindClusters = 0, 
                                      seed_RunTSNE = 1, 
                                      seed_RunUMAP = 42)
    
    # Capture meta data as a separate object
    md <- data.table::as.data.table(so@meta.data)
    rownames(md) <- md$bc
    
    md[md$cluster_type == 'Normal',]$cluster_type <- "Non-malignant" # Rename normal beads as "Non-malignant"

    # save Seurat object and meta data
    futile.logger::flog.info("saving initial Seurat object and meta data")
    saveRDS(so, paste0(output_directory, "/so.rds"))
    saveRDS(md, paste0(output_directory, "/md.rds"))
    utils::write.table(md, file="md.txt")
        
    normal_beads <- md[md$cluster_type == 'Non-malignant',]$bc
    
    gene_pos <- data.table::as.data.table(gene_pos)

    # Normalize counts, cap extreme values, and adjust for reference bead expression    
    futile.logger::flog.info("preprocessing counts")
    prep_dat <- SlideCNA::prep(so, 
                              normal_beads, 
                              gene_pos=gene_pos, 
                              chrom_ord=chrom_ord,
                              logTPM=FALSE)
    
    # Apply pyramidal average weighting scheme to expression values
    futile.logger::flog.info("getting the pyramidal weighted rolling mean of expression values")
    rm <- SlideCNA::weight_rollmean(prep_dat, 
                                   roll_mean_window)
    
    # Center Data
    futile.logger::flog.info("centering expression values")
    centered_rm <- SlideCNA::center_rm(rm)

    # Adjust for reference beads
    futile.logger::flog.info("adjusting expression values for reference beads")
    rm_adj <- SlideCNA::ref_adj(centered_rm, 
                               normal_beads)

    # Reverse log transformation
    futile.logger::flog.info("reversing log transformation of expression values")
    dat <- cbind(rm_adj[,c("GENE","chr","start","end","rel_gene_pos","length")], 
                 2 ** rm_adj[,!c("GENE","chr","start","end","rel_gene_pos","length")])
    
    # With spatial information
    if (isTRUE(spatial)) {
        # Expressional/positional binning
        futile.logger::flog.info("binning beads")
        md <- SlideCNA::bin_metadata(md, 
                                    dat, 
                                    avg_bead_per_bin, 
                                    pos, 
                                    pos_k, 
                                    ex_k, 
                                    hc_function_bin, 
                                    plot_directory)
        saveRDS(md, file="md_bin.rds")
        utils::write.table(md, file="md_bin.txt")

        # Convert data to long format and combine with metadata
        dat_long <- SlideCNA::dat_to_long(dat, md)

        # Plot spatial information
        futile.logger::flog.info("making spatial plots of binned beads")
        SlideCNA::SpatialPlot(dat_long, 
                             spatial_vars_to_plot, 
                             text_size,
                             title_size,
                             legend_size_pt, 
                             legend_height_bar, 
                             plot_directory)

        # Convert data to be by bins
        dat_bin <- SlideCNA::long_to_bin(dat_long, 
                                        plot_directory)

        # Scale data by nUMIs
        futile.logger::flog.info("scaling data by nUMIs")
        dat_bin <- SlideCNA::scale_nUMI(dat_bin, 
                                       scale_bin_thresh_hard)
        saveRDS(dat_bin, file="dat_bin_scaled.rds")
        utils::write.table(dat_bin, file="dat_bin_scaled.txt")

        # Identify CNVs 
        futile.logger::flog.info("obtaining CNA scores")
        cnv_data <- SlideCNA::prep_cnv_dat(dat_bin, 
                                          lower_bound_cnv, 
                                          upper_bound_cnv, 
                                          hc_function_cnv, 
                                          plot_directory)

        # Display CNVs across heatmap of beads x genes
        futile.logger::flog.info("making CNA plots")
        SlideCNA::cnv_heatmap(cnv_data, 
                             md, 
                             chrom_colors=chrom_colors,
                             hc_function_cnv_heatmap, 
                             plot_directory)

        # Display CNV Score Plots
        SlideCNA::quantile_plot(cnv_data, 
                               quantile_plot_cluster_label,
                               text_size,
                               title_size,
                               legend_height_bar,
                               plot_directory)
        
        SlideCNA::mean_cnv_plot(cnv_data, 
                               text_size,
                               title_size,
                               legend_height_bar,
                               plot_directory)
        
        # With just malignant beads
        futile.logger::flog.info("determining number of malignant clusters")
        best_k_malig <- SlideCNA::get_num_clust(cnv_data, 
                                               hc_function_silhouette, 
                                               max_k_silhouette, 
                                               plot_silhouette, 
                                               malig=TRUE, 
                                               k=NA, 
                                               plot_directory) 
        best_k_all <- best_k_malig + 1
        futile.logger::flog.info("Optimal number of malignant clusters: %s", best_k_malig)
        futile.logger::flog.info("Optimal number of total clusters (# of malignant clusters + 1): %s", best_k_all)
        saveRDS(best_k_malig, file="best_k_malig.rds")
        
        cnv_data2 <- cnv_data
        
        # Get clones over all beads from their CNVs
        hcl_sub_all <- SlideCNA::plot_clones(cnv_data, 
                                            md, 
                                            k=best_k_malig+1, 
                                            type='all', 
                                            chrom_colors,
                                            text_size,
                                            title_size,
                                            legend_size_pt, 
                                            legend_height_bar, 
                                            hc_function_plot_clones, 
                                            plot_directory,
                                            spatial=spatial)
        utils::write.table(hcl_sub_all, "cluster_labels_all.txt")
        
        cnv_data2$all <- merge(cnv_data$all, 
                               data.table::as.data.table(hcl_sub_all), 
                               by='variable')

        # Make binned seurat object
        futile.logger::flog.info("making Seurat object of all binned beads")
        so_bin_all <- SlideCNA::make_so_bin(so,
                                            md,
                                            hcl_sub_all)
        
        # Find DEGs and GO markers per clone over all beads
        futile.logger::flog.info("trying to find DEGs and GO markers of each cluster")
        so_clone_all <- SlideCNA::clone_so(so, 
                                          hcl_sub_all,
                                          md)

        cluster_markers_all_obj <- try(SlideCNA::find_cluster_markers(so_clone=so_bin_all, 
                                                                     type="all",
                                                                      n_markers=5,
                                                                      value="log2_expr",
                                                                     text_size=text_size,
                                                                     title_size=title_size,
                                                                     legend_size_pt=legend_size_pt,
                                                                      bin=TRUE,
                                                                     plot_directory=plot_directory))
        try(utils::write.table(cluster_markers_all_obj$markers_clone, "cluster_markers_all.txt"))
        
        if(isTRUE(use_GO_terms)) {
            go_terms_all_obj <- try(SlideCNA::find_go_terms(cluster_markers_obj=cluster_markers_all_obj,
                                                       type='all',
                                                       text_size=text_size,
                                                       title_size=title_size,
                                                       plot_directory=plot_directory))
            try(utils::write.table(go_terms_all_obj$en_clone, "go_terms_all.txt"))
        }
        else {
            go_terms_all_obj <- NULL
        }

        # Get clones over malignant beads from their CNVs
        hcl_sub_malig <- SlideCNA::plot_clones(cnv_data, 
                                              md, 
                                              k=best_k_malig, 
                                              type='malig', 
                                              chrom_colors,
                                              text_size,
                                              title_size,
                                              legend_size_pt, 
                                              legend_height_bar, 
                                              hc_function_plot_clones, 
                                              plot_directory, 
                                              spatial=spatial)
        utils::write.table(hcl_sub_malig, "cluster_labels_malig.txt")
        
        cnv_data2$malig <- merge(cnv_data$malig, 
                                 data.table::as.data.table(hcl_sub_malig), 
                                 by='variable')
        
        cnv_data <- cnv_data2 
        
        # Make binned Seurat object with malignant binned beads
        futile.logger::flog.info("making Seurat object of malignant binned beads")
        so_bin_malig <- SlideCNA::make_so_bin(so,
                                              md,
                                              hcl_sub_malig,
                                              mal=TRUE)
        
        # Find DEGs and GO markers per clone over malignant beads
futile.logger::flog.info("trying to find DEGs and GO markers of each malignant cluster")
        so_clone_malig <- SlideCNA::clone_so(so, 
                                            hcl_sub_malig, 
                                            md, 
                                            mal=TRUE)
        cluster_markers_malig_obj <- try(SlideCNA::find_cluster_markers(so_clone=so_bin_malig, 
                                                                       type="malig", 
                                                                        n_markers=5,
                                                                        value="log2_expr",
                                                                       text_size=text_size,
                                                                       title_size=title_size,
                                                                       legend_size_pt=legend_size_pt,
                                                                        bin=TRUE,
                                                                       plot_directory=plot_directory))
        try(utils::write.table(cluster_markers_malig_obj$markers_clone, "cluster_markers_malig.txt"))
        
        if(isTRUE(use_GO_terms)) {
            go_terms_malig_obj <- try(SlideCNA::find_go_terms(cluster_markers_obj=cluster_markers_malig_obj, 
                                                         type="malig", 
                                                         text_size=text_size,
                                                         title_size=title_size,
                                                         plot_directory=plot_directory))
            try(utils::write.table(go_terms_malig_obj$en_clone, "go_terms_malig.txt"))
        }
        else {
            go_terms_malig_obj <- NULL
        }
    
        cnv_data$hc_sub_all <- hcl_sub_all
        try(cnv_data$cluster_markers_all <- cluster_markers_all_obj)
        try(cnv_data$go_terms_all <- go_terms_all_obj)
        cnv_data$hc_sub_malig <- hcl_sub_malig
        try(cnv_data$cluster_markers_malig <- cluster_markers_malig_obj)
        try(cnv_data$go_terms_malig<- go_terms_malig_obj)
        
        saveRDS(cnv_data, file="cnv_data.rds")

        futile.logger::flog.info("Done!")

    }
    # Without spatial information
    else {
        # Expressional binning
        futile.logger::flog.info("binning beads")
        md <- SlideCNA::bin_metadata(md, 
                                    dat, 
                                    avg_bead_per_bin, 
                                    pos=FALSE, 
                                    pos_k, 
                                    ex_k, 
                                    hc_function_bin, 
                                    plot_directory)
        saveRDS(md, file="md_bin.rds")
        utils::write.table(md, file="md_bin.txt")

        # Convert data to long format and combine with metadata
        dat_long <- SlideCNA::dat_to_long(dat, md)

        # Convert data to be by bins
        dat_bin <- SlideCNA::long_to_bin(dat_long, 
                                        plot_directory, 
                                        spatial=spatial)

        # Scale data by nUMIs
        futile.logger::flog.info("scaling data by nUMIs")
        dat_bin <- SlideCNA::scale_nUMI(dat_bin, 
                                       scale_bin_thresh_hard)
        saveRDS(dat_bin, file="dat_bin_scaled.rds")
        utils::write.table(dat_bin, file="dat_bin_scaled.txt")

        # Identify CNVs
        futile.logger::flog.info("obtaining CNA scores")
        cnv_data <- SlideCNA::prep_cnv_dat(dat_bin, 
                                          lower_bound_cnv, 
                                          upper_bound_cnv, 
                                          hc_function_cnv, 
                                          plot_directory)

        # Display CNVs across heatmap of beads x genes
        futile.logger::flog.info("making CNA plots")
        SlideCNA::cnv_heatmap(cnv_data, 
                             md, 
                             chrom_colors,
                             hc_function_cnv_heatmap, 
                             plot_directory)

        # Get heatmap of silhouette scores across all k for all methods
        # Silhouette method to get ideal number of clusters
        futile.logger::flog.info("determining number of malignant clusters")
        best_k_malig <- SlideCNA::get_num_clust(cnv_data, 
                                               hc_function_silhouette, 
                                               max_k_silhouette, 
                                               plot_silhouette,
                                               malig=TRUE, 
                                               k=NA, 
                                               plot_directory) # With just malignant beads
        best_k_all <- best_k_malig + 1
        futile.logger::flog.info("Optimal number of malignant clusters: %s", best_k_malig)
        futile.logger::flog.info("Optimal number of total clusters (# of malignant clusters + 1): %s", best_k_all)
        saveRDS(best_k_malig, file="best_k_malig.rds")
        
        cnv_data2 <- cnv_data

        # Get clones over all beads from their CNVs
        hcl_sub_all <- SlideCNA::plot_clones(cnv_data, 
                                            md, 
                                            k=best_k_malig+1, 
                                            type='all', 
                                            chrom_colors,
                                            text_size,
                                            title_size,
                                            legend_size_pt, 
                                            legend_height_bar, 
                                            hc_function_plot_clones, 
                                            plot_directory,
                                            spatial=spatial)
        utils::write.table(hcl_sub_all, "cluster_labels_all.txt")
        
        cnv_data2$all <- merge(cnv_data$all, 
                               data.table::as.data.table(hcl_sub_all), 
                               by='variable')
        
        # Make binned seurat object
        futile.logger::flog.info("making Seurat object of all binned beads")
        so_bin_all <- SlideCNA::make_so_bin(so, 
                                            md,
                                            hcl_sub_all)
        
        # Find DEGs and GO markers per clone over all beads
        futile.logger::flog.info("trying to find DEGs and GO markers of each cluster")
        so_clone_all <- SlideCNA::clone_so(so, 
                                           hcl_sub_all, 
                                           md)

        cluster_markers_all_obj <- try(SlideCNA::find_cluster_markers(so_clone=so_bin_all, 
                                                                     type="all", 
                                                                      n_markers=5,
                                                                      value="log2_expr",
                                                                     text_size=text_size,
                                                                     title_size=title_size,
                                                                     legend_size_pt=legend_size_pt,
                                                                      bin=TRUE,
                                                                     plot_directory=plot_directory))
        try(utils::write.table(cluster_markers_all_obj$markers_clone, "cluster_markers_all.txt"))
        
        go_terms_all_obj <- try(SlideCNA::find_go_terms(cluster_markers_obj=cluster_markers_all_obj, 
                                                       type='all', 
                                                       text_size=text_size,
                                                       title_size=title_size,
                                                       plot_directory=plot_directory))
        try(utils::write.table(go_terms_all_obj$en_clone, "go_terms_all.txt"))
    
        # Get clones over malignant beads from their CNVs
        hcl_sub_malig <- SlideCNA::plot_clones(cnv_data, 
                                              md, 
                                              k=best_k_malig, 
                                              type='malig', 
                                              chrom_colors,
                                              text_size,
                                              title_size,
                                              legend_size_pt, 
                                              legend_height_bar, 
                                              hc_function_plot_clones, 
                                              plot_directory, 
                                              spatial=spatial)
        utils::write.table(hcl_sub_malig, "cluster_labels_malig.txt")
        
        cnv_data2$malig <- merge(cnv_data$malig, 
                                 data.table::as.data.table(hcl_sub_malig), 
                                 by='variable')
        
        cnv_data <- cnv_data2
        
        # Make binned seurat object with malignant beads
        futile.logger::flog.info("making Seurat object of malignant binned beads")
        so_bin_malig <- SlideCNA::make_so_bin(so,
                                              md,
                                              hcl_sub_malig,
                                              mal=TRUE)
        
        # Find DEGs and GO markers per clone over malignant beads
        futile.logger::flog.info("trying to find DEGs and GO markers of each malignant cluster")
        so_clone_malig <- SlideCNA::clone_so(so, 
                                             hcl_sub_malig, 
                                             md, 
                                             mal=TRUE)
        
        cluster_markers_malig_obj <- try(SlideCNA::find_cluster_markers(so_clone=so_bin_malig, 
                                                                       type="malig", 
                                                                        n_markers=5,
                                                                        value="log2_expr",
                                                                       text_size=text_size,
                                                                       title_size=title_size,
                                                                       legend_size_pt=legend_size_pt,
                                                                        bin=TRUE,
                                                                       plot_directory=plot_directory))
        try(utils::write.table(cluster_markers_malig_obj$markers_clone, "cluster_markers_malig.txt"))

        
        go_terms_malig_obj <- try(SlideCNA::find_go_terms(cluster_markers_obj=cluster_markers_malig_obj, 
                                                         type="malig", 
                                                         text_size=text_size,
                                                         title_size=title_size,
                                                         plot_directory=plot_directory))
        try(utils::write.table(go_terms_malig_obj$en_clone, "go_terms_malig.txt"))
        
        cnv_data$hc_sub_all <- hcl_sub_all
        try(cnv_data$cluster_markers_all <- cluster_markers_all_obj)
        try(cnv_data$go_terms_all <- go_terms_all_obj)
        cnv_data$hc_sub_malig <- hcl_sub_malig
        try(cnv_data$cluster_markers_malig <- cluster_markers_malig_obj)
        try(cnv_data$go_terms_malig<- go_terms_malig_obj)

        saveRDS(cnv_data, file="cnv_data.rds")
        
        futile.logger::flog.info("Done!")
        
    }
 
}
