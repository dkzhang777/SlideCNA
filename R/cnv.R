#' Prepare data for CNV heat map
#'
#' This function caps CNV scores, adds annotation columns for plotting,
#' performs hierarchical clustering of bins based on similar CNV score, and plots nUMI per bin 
#'
#' @param dat_bin data.table of CNV scores per bin
#' @param lower numeric float to represent the lower cap for CNV scores
#' @param upper numeric float to represent the upper cap for CNV scores 
#' @param hc_function character for which hierarchical clustering function to use
#' @param plotDir output plot directory path
#' @return A list object for downstream cnv plotting and analysis
#'         all = data.table of CNV scores of all bins x (metadata + genes)
#'         malig = data.table of CNV scores of just malignant bins x (metadata + genes)
#'         all_wide = data.frame in wide format of CNV scores of all bins x (metadata + genes)
#'         malig_wide = data.frame in wide format of CNV scores of just malignant bins x (metadata + genes)
#'         hcl = hclust object that describes the hierarchical clustering for malignant bins 
#'         hcl_all = hclust object that describes the hierarchical clustering for all bins 

#' @export
prep_cnv_dat <- function(dat_bin, 
                         lower=0.6, 
                         upper=1.4, 
                         hc_function = 'ward.D2', 
                         plotDir) {                                  
    # separate data into all (malignant + non-malignant), malignant, and non-malignant cells
    malig_cells=unique(dat_bin[cluster_type == 'Malignant',]$variable)
    ref_cells=unique(dat_bin[cluster_type == 'Non-malignant',]$variable)
    all_cells=unique(dat_bin$variable)
    
    dat_malig <- dat_bin[variable %in% malig_cells]
    dat_ref <- dat_bin[variable %in% ref_cells]
    
    # cap data
    dat_bin[value>upper]$value <- upper
    dat_malig[value>upper]$value <- upper
    dat_bin[value<lower]$value <- lower
    dat_malig[value<lower]$value <- lower
    
    # factorize cells, chrs, and add order column
    dat_malig[,variable:=factor(variable,levels=malig_cells),]
    dat_malig$chr <- factor(dat_malig$chr, 
                            levels = unique(dat_malig$chr))
    dat_malig[,total_order:=c(1:nrow(dat_malig)),]
    dat_bin[,variable:=factor(variable,levels=all_cells),]
    dat_bin$chr <- factor(dat_bin$chr, 
                          levels = unique(dat_bin$chr))
    dat_bin[,total_order:=c(1:nrow(dat_bin)),]
    
    # transform to wide format for a bins vs genes data frame
    malig_wide=as.data.frame(data.table::dcast(dat_malig,variable~chr+plot_order,
                                               value.var = "value"))
    row.names(malig_wide)=malig_wide$variable
    malig_wide$variable=NULL
    malig_wide <- malig_wide[,stringr::str_sort(colnames(malig_wide), numeric = TRUE)]

    all_wide=as.data.frame(data.table::dcast(dat_bin,variable~chr+plot_order,
                                             value.var = "value"))
    row.names(all_wide)=all_wide$variable
    all_wide$variable=NULL
    all_wide <- all_wide[,stringr::str_sort(colnames(all_wide),
                                   numeric = TRUE)]

    # hierarchical clustering of bins by CNV score
    hcl=stats::hclust(stats::dist(malig_wide), method=hc_function)
    hcl_all=stats::hclust(stats::dist(all_wide), method=hc_function)
    
    # plot dendrogram of hierarchical clustering
    grDevices::pdf(file = paste0(plotDir,"/hcl_cnv.pdf"), width = 10, height = 6)
    graphics::par(mfrow=c(1,2)) 
    plot(hcl,labels=FALSE)
    plot(hcl_all,labels=FALSE)
    grDevices::dev.off()
    
    graphics::par(mfrow=c(1,2)) 
    plot(hcl,labels=FALSE)
    plot(hcl_all,labels=FALSE)
    
    dat_malig[,variable:=factor(variable,levels=hcl$labels[hcl$order])]
    dat_bin[,variable:=factor(variable,levels=hcl_all$labels[hcl_all$order])]

    # plot distribution of values across number of beads/bin and umi/bin
    grDevices::pdf(file = paste0(plotDir,"/beads_umi_per_bin.pdf"), width = 10, height = 6)
    graphics::par(mfrow=c(1,2)) 
    graphics::boxplot(dat_malig$value ~ dat_malig$N_bin)
    graphics::boxplot(dat_malig$value ~ dat_malig$umi_bin)
    grDevices::dev.off()
    
    graphics::par(mfrow=c(1,2)) 
    graphics::boxplot(dat_malig$value ~ dat_malig$N_bin)
    graphics::boxplot(dat_malig$value ~ dat_malig$umi_bin)
    
    cnv_data <- list("all" = dat_bin, 
                     "malig" = dat_malig, 
                     "all_wide" = all_wide, 
                     "malig_wide" = malig_wide, 
                     "hcl"=hcl,
                     "hcl_all" = hcl_all)
    return(cnv_data)
}
utils::globalVariables(c("cluster_type", "variable", "value", "total_order"))

#' Plot CNV scores on a heat map
#'
#' This function prepares data for plotting and makes a heat map of CNV scores per bead across all genes
#'
#' @param cnv_data list object of cnv data from SlideCNA::prep_cnv_dat()
#' @param md data.table of metadata of each bead
#' @param chrom_colors vector of colors labeled by which chromosome they correspond to
#' @param hc_function character for which hierarchical clustering function to use
#' @param plotDir output plot directory path

#' @export
cnv_heatmap <- function(cnv_data, 
                        md, 
                        chrom_colors, 
                        hc_function = 'ward.D2', 
                        plotDir) {
    
    # get chromosome names
    annot_col <- data.frame(chr=matrix(unlist(strsplit(colnames(cnv_data$all_wide), "_")), 
                                       ncol=2, 
                                       byrow=TRUE)[,1], 
                            row.names=colnames(cnv_data$all_wide))
    
    # Select columns to plot
    md_sub <- unique(md[,c('bin_all', 'cluster_type')])
    annot_row <- data.frame(md_sub$cluster_type, row.names=md_sub$bin_all)
    colnames(annot_row)[1] = 'Tissue Type'
    
    # Order by plot order (genomic position)
    cnv_data$all_wide <- cnv_data$all_wide[order(as.numeric(row.names(cnv_data$all_wide))),]
    cnv_data$all_wide <- rbind(cnv_data$all_wide[1:nrow(cnv_data$malig_wide),]
                          [match(rownames(cnv_data$malig_wide),
                                 rownames(cnv_data$all_wide[1:nrow(cnv_data$malig_wide),])),],
                               cnv_data$all_wide[(nrow(cnv_data$malig_wide)+1):nrow(cnv_data$all_wide),])
    
    # Choose colors to represnt tissue types
    ann_colors = list("Tissue Type" = c("Malignant"="#F8776D", "Non-malignant"="#00BFC4"), chr = chrom_colors)
    
    # Plot heat map
    grDevices::png(file = paste0(plotDir,"/cnv_heatmap.png"), width = 1500, height = 1000) # png version
    print(pheatmap::pheatmap(cnv_data$all_wide, 
                             color = grDevices::colorRampPalette(c("navy", "white","firebrick3"))(50), 
                             cluster_rows = TRUE, 
                             cluster_cols = FALSE, 
                             clustering_distance_rows = "euclidean",
                             clustering_method = hc_function, 
                             annotation_row = annot_row, 
                             annotation_col = annot_col,
                             annotation_colors=ann_colors, 
                             show_rownames=FALSE, 
                             show_colnames=FALSE, 
                             filename=paste0(plotDir,"/cnv_heatmap.png")))
    grDevices::dev.off()
 
    # Plot heat map (pdf)
    pheatmap::pheatmap(cnv_data$all_wide, 
                             color = grDevices::colorRampPalette(c("navy", "white","firebrick3"))(50), 
                             cluster_rows = TRUE, 
                             cluster_cols = FALSE, 
                             clustering_distance_rows = "euclidean",
                             clustering_method = hc_function, 
                             annotation_row = annot_row, 
                             annotation_col = annot_col,
                             annotation_colors=ann_colors, 
                             show_rownames=FALSE, 
                             show_colnames=FALSE, 
                             filename=paste0(plotDir,"/cnv_heatmap.pdf"),
                              width=18,
                              height=15)
}

#' Plot CNV score quantiles per bin and per chromosome
#'
#' This function colors and plots each bin by its CNV score quantiles (min, 1st quartile, median, 3rd quartile, max)
#' on spatial coordinates for each chromosome
#'
#' @param cnv_data list object of cnv data from SlideCNA::prep_cnv_dat()
#' @param cluster_label character string of which column name to keep 
#' @param text_size integer of text size for ggplot
#' @param title_size integer of title size for ggplot
#' @param legend_height_bar integer of bar height of legend for ggplot
#' @param plotDir output plot directory path
#’ @Import ggplot2

#' @export
quantile_plot <- function(cnv_data, 
                          cluster_label="seurat_clusters", 
                          text_size, 
                          title_size, 
                          legend_height_bar, 
                          plotDir) {
    plot_data <- cnv_data$all
    plot_data$chr <- factor(plot_data$chr, levels = unique(plot_data$chr))
    
    # get quantile values
    min_plot <- plot_data[,new_value:=min(value),by=c("variable","chr")]
    min_plot$level <- '1'
    first_plot <- plot_data[,new_value:=stats::quantile(value, .25),by=c("variable","chr")]
    first_plot$level <- '2'
    med_plot <- plot_data[,new_value:=stats::median(value),by=c("variable","chr")]
    med_plot$level <- '3'
    third_plot <- plot_data[,new_value:=stats::quantile(value, .75),by=c("variable","chr")]
    third_plot$level <- '4'
    max_plot <- plot_data[,new_value:=max(value),by=c("variable","chr")]
    max_plot$level <- '5'

    # combined data.frame of quamtiles
    quant_combined <- rbind(min_plot, first_plot, med_plot, third_plot, max_plot)
    quant_combined=dplyr::distinct(dplyr::select(quant_combined,
                                                 c(variable, chr, pos_x, pos_y, cluster_label, new_value, level)))
    # Plot CNV scores across space
    legend_title='CNV Score'
    gg <- ggplot2::ggplot(quant_combined) +
    geom_point(aes(pos_x, pos_y, color = new_value)) + 
    facet_wrap(chr~level,scales="free_x") +
    scale_color_gradient2(midpoint=1, low="blue", mid="white", high="red") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "black"),
          axis.title=element_text(size=title_size, face='bold'), 
          strip.text.x = element_text(size=text_size),
          legend.title=element_text(size=title_size, face='bold'),
          legend.text=element_text(size=text_size),
          legend.key.height = unit(legend_height_bar,"cm"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Position X") + 
    ylab("Position Y") + 
    labs(color = legend_title)
    
    grDevices::png(file = paste0(plotDir,"/cnv_score_quantiles.png"), width = 1000, height = 1200) # png version
    print(gg)
    grDevices::dev.off()
    print(gg)
}                              
utils::globalVariables(c("new_value", "value", "variable", "chr", "pos_x", "pos_y", "level"))

#' Plot mean CNV scores per bin and per chromosome
#'
#' This function colors and plots each bin by its mean CNV score
#' on spatial coordinates for each chromosome
#'
#' @param cnv_data list object of cnv data from SlideCNA::prep_cnv_dat()
#' @param text_size integer of text size for ggplot
#' @param title_size integer of title size for ggplot
#' @param legend_height_bar integer of bar height of legend for ggplot
#' @param plotDir output plot directory path
#’ @Import ggplot2

#' @export
mean_cnv_plot <- function(cnv_data, 
                          text_size, 
                          title_size, 
                          legend_height_bar, 
                          plotDir) {
    
#     plot_data <- data.table()
#     for (col in colnames(cnv_data$all)) {
#         plot_data <- cbind(plot_data, cnv_data$all[[col]])
#     }
#     colnames(plot_data) <- colnames(cnv_data$all)
    
    plot_data <- data.table::copy(cnv_data$all)
    # add new column that takes the mean cnv score across all genes
    mean_plot <- plot_data[,cnv_score:=mean(value),by=c("variable","chr")]
    
    # plot CNV scores
    legend_title='CNV Score'
    gg <- ggplot2::ggplot(mean_plot) +
    geom_point(aes(pos_x, pos_y, color = cnv_score)) + 
    facet_wrap(~chr,scales="free_x") +
    scale_color_gradient2(midpoint=1, low="blue", mid="white", high="red") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "black"), 
          axis.text=element_text(size=text_size), 
          axis.title=element_text(size=title_size, face='bold'), 
          strip.text.x = element_text(size=text_size),
          legend.title=element_text(size=title_size, face='bold'), 
          legend.text=element_text(size=text_size),
          legend.key.height = unit(legend_height_bar,"cm"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("Position X") + ylab("Position Y") + labs(color = legend_title) 
    
    grDevices::png(file = paste0(plotDir,"/mean_cnv_score.png"), width = 1200, height = 1200) # png version
    print(gg)
    grDevices::dev.off()
    print(gg)
}                              
utils::globalVariables(c("cnv_score", "value", "pos_x", "pos_y"))