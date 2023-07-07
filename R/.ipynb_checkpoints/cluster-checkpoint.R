#' Find optimal number of clusters
#'
#' This function uses the Silhouette Method applied to CNV scores to determine the 
#' best number of clusters to divide the binned beads into
#'
#' @param data cnv_data list object of cnv data from SlideCNA::prep_cnv_dat()
#' @param hc_func character string for which hierarchical clustering function to use
#' @param max_k integer of number max number of clusters to evaluate (2:max_k)
#' @param plot TRUE if plotting silhoutte scores per cluster
#' @param malig TRUE if only using malignant bins and FALSE if using all bins
#' @param k integer of optimal number of clusters, if known, and NA if not known
#' @param plotDir output plot directory path
#' @return An integer representing the number of clusters that optimizes the silhouette score

#' @export
get_num_clust <- function(data, 
                          hc_func='ward.D2', 
                          max_k=10, 
                          plot=TRUE, 
                          malig=FALSE, 
                          k=NA, 
                          plotDir) {
    
    if (malig == TRUE) {
        type = 'malig'
    }
    else {
        type = 'all'
    }
    
    # Malignant cells only
    if (malig == TRUE) {
        hcl <- data$hcl
        expr_data <- data$malig_wide
    }
    # All cells
    else {
        hcl <- data$hcl_all
        expr_data <- data$all_wide
    }
    
    # Sometimes thresholding results in columns with all the same value (e.g. 0.6 or 1.4)
    # The following cannot work w/ 0 variance, so add 0.01 to top row of these 0 variance columns
    expr_data[1,which(apply(expr_data, 2, var)==0)] <- expr_data[1,which(apply(expr_data, 2, var)==0)] + 0.01
    
    # if no input k, automatically find optimal k
    if (is.na(k)) {
        # Get the k (number of clusters) that maximizes silhouette score
        best_score <- -Inf
        for (k in 2:max_k) {
            sil_cl <- cluster::silhouette(cutree(hcl, k=k), 
                                          dist(expr_data))
            avg_sil_score <- mean(sil_cl[,3])
            if (avg_sil_score > best_score) { #in case of ties, takes lower k
                best_score <- avg_sil_score
                best_k <- k
            }
        }
    }
    # if input k given, set best k to be that value
    else {
        best_k <- k
    }

    
    # Optional plotting
    if (isTRUE(plot)) {
        
        # Cluster embedding plot
        fviz <- factoextra::fviz_cluster(list(data=expr_data, 
                                              cluster=cutree(hcl, k=best_k)))
        pdf(file = paste0(plotDir, "/", type, "_clone_dim_red.pdf"), width = 8, height = 6)
        print(fviz)
        dev.off()
        print(fviz)
        
        # Silhouette plot
        sil_cl <- cluster::silhouette(cutree(hcl, k=best_k), 
                                      dist(expr_data))
        pdf(file = paste0(plotDir,"/", type, "_clone_silhouette_plot.pdf"), width = 8, height = 6)
        plot(sil_cl, border=NA)
        dev.off()
        plot(sil_cl, border=NA)

    }
    return(best_k)
}            

#' Plot cluster/clone information
#'
#' This function plots cluster dendrograms, spatial assignment, and the CNV heat map
#'
#' @param data cnv_data list object of cnv data from SlideCNA::prep_cnv_dat()
#' @param md data.table of metadata of each bead
#' @param k integwer of number of clusters/clones
#' @param type character string, being "all" if using all binned beads, or "malig" if just malignant binned beads
#' @param hc_function character string for which hierarchical clustering function to use
#' @param plotDir output plot directory path
#' @param spatial TRUE if using spatial information
#' @return A hierarchical clustering object of the clusters
#' @export
#' @import ggplot2

#' @export
plot_clones = function(cnv_data, 
                       md, 
                       k, 
                       type, 
                       chrom_colors, 
                       text_size, 
                       title_size, 
                       legend_size_pt, 
                       legend_height_bar, 
                       hc_function = 'ward.D2', 
                       plotDir, 
                       spatial=TRUE) {
        
    # Initialize data based on CNV object
    if (type == 'all') {
        sub = cnv_data$all
        sub_wide = cnv_data$all_wide
        hcl = cnv_data$hcl_all
    }
    else if (type == 'malig') {
        sub = cnv_data$malig
        sub_wide = cnv_data$malig_wide
        hcl = cnv_data$hcl
    }
    
    # Plot hierarchical clustering dendrogram of how clones are separated    
    gg_dend <- ggplot2::ggplot(dendextend::as.ggdend(as.dendrogram(hcl) %>%
      dendextend::set("branches_k_color", k = k)))

    pdf(file = paste0(plotDir,"/", type, "_", k, "_clones_dend.pdf"), width = 6, height = 8)
    print(gg_dend)
    dev.off()
    
    # Cut tree into k clones
    hcl_sub <- as.data.frame(cutree(hcl, k = k))
    colnames(hcl_sub) <- c('clone')
    hcl_sub$variable <- as.factor(rownames(hcl_sub))
    
    sub_clone <- merge(sub, hcl_sub, by='variable')

    # Only malignant binned beads
    if (type=='malig') {
        normal_bins <- setdiff(cnv_data$all$variable, cnv_data$malig$variable)
        
        sub_normal <- cnv_data$all[cnv_data$all$variable %in% normal_bins, ]
        sub_normal[, clone:=(k+1)]
        # sub_normal$clone <- k+1 # Make normal binned beads another clone
        sub_all <- rbind(sub_clone, sub_normal)
        
        if (spatial==TRUE) {
            sub_all_unique <- unique(sub_all[,c('variable', 'pos_x', 'pos_y', 'clone')])
        }
        else {
            sub_all_unique <- unique(sub_all[,c('variable', 'clone')])
        }
    }    
    
    # Plot and color clones spatially
    if (spatial==TRUE) {
        sub_clone_unique <- unique(sub_clone[,c('variable', 'pos_x', 'pos_y', 'clone')])
        legend_title = 'Cluster'
        
        if (type=='all') {
            df_to_plot <- sub_clone_unique
            spatial_colors <- scales::hue_pal()(k)
        }
        else if (type=='malig') {
            df_to_plot <- sub_all_unique
            spatial_colors <- c(scales::hue_pal()(k), "black")
        }

        gg_spatial <- ggplot2::ggplot(df_to_plot) +
        geom_point(aes(pos_x, pos_y, color = as.factor(clone)), size=2) +
        theme_bw() + 
        theme(axis.text=element_text(size=text_size),
              axis.title=element_text(size=title_size,face='bold'),
              legend.title=element_text(size=title_size, face='bold'),
              legend.text=element_text(size=text_size),
              plot.background = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        xlab("Position X") + 
        ylab("Position Y") + labs(color = legend_title) + 
        coord_fixed(ratio=1) +
        guides(color = guide_legend(override.aes = list(size = legend_size_pt))) +
        scale_color_manual(values=spatial_colors)

        pdf(file = paste0(plotDir,"/", type, "_", k, "_clones_spatial.pdf"), width = 6, height = 8)
        print(gg_spatial)
        dev.off()
        print(gg_spatial)    
    }
    else {
        sub_clone_unique <- unique(sub_clone[,c('variable', 'clone')])
    }
    
    # Add annotations for column information
    # Chromsomes
    annot_col <- data.frame(chr=matrix(unlist(strsplit(colnames(sub_wide), "_")), 
                                       ncol=2,
                                       byrow=TRUE)[,1], 
                            row.names=colnames(sub_wide)) 
    md_sub <- unique(md[,c('bin_all', 'cluster_type')])
    annot_row <- data.frame(cluster_type=md_sub$cluster_type, # tissue type
                            variable=md_sub$bin_all) # bin 
    annot_row <- merge(annot_row, 
                       as.data.frame(sub_clone_unique[,c('variable', 'clone')]),
                       by="variable")
    rownames(annot_row) <- annot_row$variable
    annot_row <- annot_row[,c('cluster_type', 'clone')]
    annot_row$clone <- factor(annot_row$clone)
    colnames(annot_row) = c('Tissue Type', 'Cluster')
    
    # Make each clone its own color
    hues <- scales::hue_pal()(k)
    clone_colors = c()
    for (i in 1:k) {
        clone_colors <- append(clone_colors, hues[i])
    }
    names(clone_colors)<- as.character(c(1:k)) 
    
    # Set colors
    ann_colors = list(
    'Tissue Type' = c("Malignant"="#F8776D", "Non-malignant"="#00BFC4"),                       
    'Cluster' = clone_colors,
    chr = chrom_colors)

    # Plot CNV heatmap with clone labelling    
    png(file = paste0(plotDir,"/", type, "_", k, "_clones_cnv_heatmap.png"), width = 1500, height = 1000) #png version
    print(pheatmap::pheatmap(sub_wide, 
                             color = colorRampPalette(c("navy", "white","firebrick3"))(50), 
                             cluster_rows = TRUE, 
                             cluster_cols = FALSE, 
                             clustering_distance_rows = "euclidean",
                             clustering_method = hc_function, 
                             annotation_row = annot_row, 
                             annotation_col = annot_col,
                             annotation_colors=ann_colors, 
                             show_rownames=FALSE, 
                             show_colnames=FALSE, 
                             filename= paste0(plotDir,"/", type, "_", k, "_clones_cnv_heatmap.png")))
    dev.off()
    
     # Plot CNV heatmap with clone labelling (pdf)
     pheatmap::pheatmap(sub_wide, 
                             color = colorRampPalette(c("navy", "white","firebrick3"))(50), 
                             cluster_rows = TRUE, 
                             cluster_cols = FALSE, 
                             clustering_distance_rows = "euclidean",
                             clustering_method = hc_function, 
                             annotation_row = annot_row, 
                             annotation_col = annot_col,
                             annotation_colors=ann_colors, 
                             show_rownames=FALSE, 
                             show_colnames=FALSE, 
                             filename= paste0(plotDir,"/", type, "_", k, "_clones_cnv_heatmap.pdf"),
                       width=18,
                       height=15)
    
    return(hcl_sub)
}
            
#' Add clone information to meta data of seurat object    
#'
#' This function adds another column for cluster designation to a seurat object's meta data
#'
#' @param data so seurat object of beads and their meta data
#' @param hcl_sub hierarchical clustering object of cluster assignemnt as outputted from SlideCNA::plot_clones()
#' @param md data.table of metadata of each bead
#' @param mal TRUE if only using malignant beads
#' @return A seurat object updated with clone information

#' @export
clone_so <- function(so, 
                     hcl_sub, 
                     md, 
                     mal=FALSE) {
    
    # Only malignant beads
    if (isTRUE(mal)) {
        md <- md[cluster_type=='Malignant']
    }
    
    so_clone <- so
    hcl_sub_int <- hcl_sub
    hcl_sub_int$variable <- as.integer(hcl_sub_int$variable)

    # Add cluster designation to meta data
    md_clone <- merge(md, hcl_sub_int, by.x='bin_all',by.y='variable')
    md_merge <- merge(so_clone@meta.data, md_clone)
    rownames(md_merge) <- md_merge$bc
    # Update seurat object
    so_clone=Seurat::AddMetaData(so_clone, metadata = md_merge)
    
    # Subset to only malignant beads
    if (mal==TRUE) {
        so_clone <- subset(so_clone, cells =(data.table::as.data.table(so_clone@meta.data)[cluster_type=='Malignant'])$bc)
    }
    
    so_clone=Seurat::SetIdent(so_clone, value = so_clone@meta.data$clone)
    
    return(so_clone)
}       

#' Find and plot top n DEGs per cluster  
#'
#' This function uses Seurat's marker finding capability to find DEGs of each cluster
#'
#' @param data so seurat object with clone information
#' @param type character string that is 'all' if using malignant and normal clusters and 'malig'
#'        if just using malignant clusters
#' @param md logfc.threshold numeric float that is seurat parameter,
#'        representing the minimum log2 fold change for DEGs to be significant
#' @param min.pct numeric Seurat function parameter 
#' @param only.pos TRUE if only using DEGs with positive log2 fold change
#' @param n_markers integer of number of top DEGs to plot/use
#' @param plotDir output plot directory path
#' @return A list object with cluster marker information
#'         markers_clone = data.table of all cluster markers
#'         top_markers_clone = data.table of just top cluster markers
#'         top_clone_vis = data.frame formatted for plot visualization of top cluster markers
#' @import ggplot2

#' @export
find_cluster_markers <- function(so_clone, 
                                 type, 
                                 logfc.threshold=0.2, 
                                 min.pct=0, 
                                 only.pos=TRUE, 
                                 n_markers=5, 
                                 text_size, 
                                 title_size, 
                                 legend_size_pt, 
                                 plotDir) {
    
    # get all marker genes per cluster
    markers_clone=data.table::as.data.table(Seurat::FindAllMarkers(so_clone,
                                                                   logfc.threshold=logfc.threshold, 
                                                                   min.pct=min.pct,
                                                                   only.pos=only.pos))
    # select top n markers 
    top_markers_clone <- markers_clone[order(avg_log2FC,decreasing=TRUE),
                                       .SD[1:n_markers,], by="cluster"][p_val_adj<.05,]

    # recreate counts table
    counts_re <- as.data.frame(Seurat::GetAssayData(object = so_clone, slot = "counts"))
    
    # Put the top DEGs in a form that can be easily visualized, by getting the average expression of each DEG per cluster and the percent of beads in that cluster it is expressed in (count>0)
    genes <- c()
    gene_pct <- c()
    avg_expr <- c()
    cluster <- c()
    for (gene in top_markers_clone$gene) {
        for (clust in unique(top_markers_clone$cluster)) {
            counts_clust <- dplyr::select(counts_re, 
                                          tidyselect::all_of(so_clone@meta.data[so_clone@meta.data$clone==clust,]$bc))
            genes <- c(genes,gene)
            cluster <- c(cluster, clust)
            gene_pct <- c(gene_pct, length(which(counts_clust[gene,]>0))/(length(counts_clust[gene,])))
            avg_expr <- c(avg_expr, mean(as.numeric(counts_clust[gene,])))
        }
    }
    top_clone_vis <- data.frame(gene=genes,
                                cluster=cluster,
                                gene_pct=gene_pct,
                                avg_expr=avg_expr)
    
    # Plot of top DEGs per cluster and their average expressions + percentage of beads they are expressed in per cluster
    size_title = "Percent Expressed"
    color_title = "Average Expression"

    gg <- ggplot2::ggplot(top_clone_vis, 
                          aes(x = gene, 
                              y = cluster, 
                              size = gene_pct, 
                              colour = avg_expr)) +
    geom_point() + 
    scale_colour_viridis_c(option='C') + 
    scale_size_area(max_size = 15) + 
    theme_bw() +
    theme(axis.text=element_text(size=text_size), 
          axis.text.x=element_text(angle = 45, vjust=.9, hjust=.9), 
          axis.title=element_text(size=title_size, face='bold'),
          legend.title=element_text(size=title_size, face='bold'),
          legend.text=element_text(size=text_size),
          plot.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Gene") + 
    ylab("Cluster") + 
    labs(color = color_title, 
         size = size_title) + 
    guides(color = guide_legend(override.aes = list(size = legend_size_pt)))
    
    pdf(file = paste0(plotDir,"/", "top_", n_markers, "_markers_", type, ".pdf"), 
        width = 10, height = 6)
    print(gg)
    dev.off()
    print(gg)

    # Create combined cluster object
    if (type=='all') {
        cluster_markers_all_obj <- list("markers_clone"=markers_clone, 
                                        "top_markers_clone"=top_markers_clone,
                                        "top_clone_vis"=top_clone_vis) 
        return(cluster_markers_all_obj)
    }
    else if (type=='malig') {
        cluster_markers_malig_obj <- list("markers_clone"=markers_clone, 
                                          "top_markers_clone"=top_markers_clone,
                                          "top_clone_vis"=top_clone_vis) 
        return(cluster_markers_malig_obj)
    }                             
}                            

#' Find and plot top n GO-enriched terms per cluster  
#'
#' This function utilizes cluster-specific DEGs to identify cluster-specifc GO biological processes
#' and plots these if they occur
#'
#' @param data cluster_markers_obj list object with cluster marker information
#' @param type character string that is 'all' if using malignant and normal clusters and 'malig'
#'        if just using malignant clusters
#' @param n_terms integer of number of top DEGs to plot/use
#' @param plotDir output plot directory path
#' @return A list object with cluster GO term information
#'         en_clone = data.table of cluster GO terms
#'         top_en_clone = data.table of just top cluster GO terms
#' @import ggplot2
                             
#' @export
find_go_terms <- function(cluster_markers_obj, 
                          type, 
                          n_terms=5, 
                          text_size, 
                          title_size, 
                          plotDir) {
    
    # Get enriched GO terms from marker genes for each cluster
    en_clone=cluster_markers_obj$markers_clone[order(avg_log2FC, 
                                                     decreasing=TRUE), 
                                               SlideCNA::run_enrichr(gene,50),
                                               by="cluster"]   
    
    # Select top n_terms GO terms
    top_en_clone <-en_clone[GO_Biological_Process_2018.Adjusted.P.value<0.05][order(GO_Biological_Process_2018.Adjusted.P.value),.SD[1:n_terms],by="cluster"] 

    top_en_clone$neglog10p <- -log10(top_en_clone$GO_Biological_Process_2018.Adjusted.P.value)
    top_en_clone$GO_Biological_Process_2018.Term <- gsub("\\s*\\([^\\)]+\\)","",top_en_clone$GO_Biological_Process_2018.Term) # Rename to look nicer                            
    top_en_clone$cluster <- factor(top_en_clone$cluster, 
                                   levels = sort(as.numeric(levels(top_en_clone$cluster))))

    # Plot the top GO enriched terms per cluster with bars representing -log 10 p-value
    gg <- ggplot2::ggplot(top_en_clone, 
                          aes(x=stringr::str_wrap(GO_Biological_Process_2018.Term,20),
                              y=neglog10p)) + 
    geom_bar(stat="identity") + 
    facet_wrap(~cluster,scale="free_y") + 
    theme_bw() + 
    theme(axis.text=element_text(size=text_size), 
          axis.title=element_text(size=title_size,face='bold'), 
          strip.text.x = element_text(size = text_size),
          plot.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("GO Biological Process 2018 Term") + 
    ylab((("-Log_10(p)"))) +
    coord_flip()
    
    pdf(file = paste0(plotDir,"/", "top_", n_terms, "_go_terms_", type, ".pdf"), width = 15, height = 6)
    print(gg)
    dev.off()
    print(gg)   
    
    # Combine data structures in list object format
    if (type=='all') {
        go_terms_all_obj <- list("en_clone"=en_clone, 
                                        "top_en_clone"=top_en_clone) 
        return(go_terms_all_obj)
    }
    else if (type=='malig') {
        go_terms_malig_obj <- list("en_clone"=en_clone, 
                                        "top_en_clone"=top_en_clone) 
        return(go_terms_malig_obj)
    }
}                             

### Subfunction to get significantly enriched GO terms given a set of signfiicant beads and genes
#' @export
run_enrichr=function(genes,
                     ngenes) {
    res=data.table::as.data.table(enrichR::enrichr(genes[1:ngenes]
                                                   ,databases = c("GO_Biological_Process_2018")))
    return(res)
}

