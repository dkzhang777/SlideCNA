#' Spatio-molecular binning of relative expression intensities
#'
#' This function combines metadata with binned relative expression intensities 
#'
#' @param md data.table of metadata of each bead
#' @param dat data.table of smoothed relative expression intensities 
#' @param avg_bead_per_bin integer of average number of beads there should be per bin 
#' @param pos TRUE if doing spatial and expressional binning, FALSE if just expressional binning
#' @param pos_k positional weight
#' @param ex_k expressional weight
#' @param hc_function hierarchical clustering function
#' @param plotDir output plot directory path
#' @return A data.table of bead metadata combined with binned expression intensities for all genes for all beads

bin_metadata <- function(md, 
                         dat, 
                         avg_bead_per_bin=12, 
                         pos=TRUE, 
                         pos_k=55, 
                         ex_k=1, 
                         hc_function='ward.D2',
                         plotDir) {
    
    # bin malignant beads
    md_mal <- md[cluster_type=='Malignant']
    dat_mal <- cbind(dat[,c("GENE", "chr", "start", "end", "rel_gene_pos", "length")], 
                     dat[,which(dat[,c(colnames(dat)%in%md_mal$bc)]), with=FALSE])
    
    n_bin_mal <- round(nrow(md_mal)/avg_bead_per_bin)
    md_mal <- package::bin(dat_mal, md_mal, n_bin_mal, pos, pos_k, ex_k, hc_function, plotDir)
    
    # bin reference beads
    md_ref <- md[cluster_type=='Normal']
    dat_ref <- cbind(dat[,c("GENE", "chr", "start", "end", "rel_gene_pos", "length")], 
                     dat[,which(dat[,c(colnames(dat)%in%md_ref$bc)]), with=FALSE])
    
    n_bin_ref <- round(nrow(md_ref)/avg_bead_per_bin)
    md_ref <- package::bin(dat_ref, md_ref, n_bin_ref, pos, pos_k, ex_k, hc_function, plotDir)
    md_ref$bin_all <- md_ref$bin_all + max(md_mal$bin_all) 
    
    md <- rbind(md_mal, md_ref)
    md[,pct_mal:=sum(icluster_type)/.N,by=bin_all]
    return(md)
    
}   

### Subfunction of bin_metadata() for Expression/Positional Binning 
                              
bin <- function(dat, 
                md, 
                k, 
                pos=TRUE, 
                pos_k=55, 
                ex_k=1, 
                hc_function = 'ward.D2', 
                plotDir) {
    
    dat_var <- t(dat[,!c("chr", "start", "end", "rel_gene_pos", "length")]) 
    colnames(dat_var) <- dat_var[1,]
    dat_var <- as.data.frame(dat_var[-1,])
    dat_var$bc <- rownames(dat_var)
    dat_var <- as.data.frame(merge(dplyr::select(md, bc), dat_var, by="bc"))
    row.names(dat_var) <- dat_var$bc
    dat_var <- dat_var[,-1]  
    dat_var <- scale(sapply(dat_var, as.numeric))
    expr_distance <- dist(dat_var) # get distances from expression matrix
    
    md[,icluster_type:=ifelse(cluster_type=='Normal',0,1)]
    icluster_type <- as.data.frame(dplyr::select(md, bc, icluster_type))
    row.names(icluster_type) <- icluster_type$bc
    icluster_type <- icluster_type[,-1]
    
    # If using spatial information
    if (pos==TRUE) {
        pos <- as.data.frame(dplyr::select(md, bc, pos_x, pos_y))
        row.names(pos) <- pos$bc
        pos <- pos[,-1]
        pos_distance <- dist(pos) # get distances from position matrix
        
        # linearly combine expression and position matrices
        distance <- pos_k * pos_distance + ex_k * expr_distance
    }
    else {
        distance <- expr_distance
    }
    
    # hierarchical clustering on combined distances
    hcl <- hclust(distance, method=hc_function)
    pdf(file = paste0(plotDir, "/bin_hclust_dend.pdf"), width = 10, height = 6)
    plot(hcl)
    dev.off()
    
    # get bins 
    hcl_sub <- cutree(hcl, k = k)
    
    # update bead metadata with bin designations
    new_md <- dplyr::mutate(md, bin_all = hcl_sub)
    new_md[,N_bin:=.N,by=bin_all] #.N is number of instances that bin combo exists in dataset
    new_md[,umi_bin:=sum(nCount_RNA),by=bin_all]
    
    return(new_md)
}

#' Convert data to long format and add in metadata
#'
#' This function will create rows for each bead and gene combination, adding in new metadata with bin designations
#'
#' @param dat data.table of smoothed relative expression intensities 
#' @param md data.table of metadata per bead 
#' @return A data.table of bead expression intensities per gene with metadata in long format 

dat_to_long <- function(dat, 
                        md) {
    # change to long format
    dat_long=data.table::as.data.table(reshape2::melt(dat,id.vars = c("GENE","chr","start","end","rel_gene_pos", "length")))
    # add metadata
    dat_long=merge(dat_long,md,by.x="variable",by.y="bc")
    return(dat_long)
}  

#' Spatial plots of meta data
#'
#' This function will plot information about beads and bins on x and y coordinates
#'
#' @param dat_long data.table of bead expression intensities per gene with metadata in long format
#' @param vars character vector of features to plot/columns of metadata
#' @param plotDir output plot directory path

SpatialPlot <- function(dat_long, 
                        vars=NULL, 
                        text_size, 
                        title_size, 
                        legend_size_pt, 
                        legend_height_bar, 
                        plotDir) {
    dat_distinct <- dplyr::distinct(dplyr::select(dat_long, 
                                                  c("variable", "pos_x", "pos_y", vars)))
    
    for (var in vars) {
        
        if (var == "seurat_clusters") {
            legend_title = "Seurat Cluster"
        }
        else if (var == "nmf_ct") {
            legend_title = "NMF Cell Types"
        }
        else if (var == "Cross_Section") {
            legend_title = "Tumor Subsample"
        }
        else if (var == "bin_all") {
            legend_title = "Bin"
        }
        else if (var == "N_bin") {
            legend_title = "Number of Beads per Bin"
        }
        else if (var == "umi_bin") {
            legend_title = "Number of UMIs per Bin"
        }
        else if (var == "cluster_type") {
             legend_title = "Tissue Type"
        }
        
        if (var %in% c('seurat_clusters', 'nmf_ct', 'cluster_type', 'Cross_Section')) {
            gg <- ggplot2::ggplot(dat_distinct) +
            geom_point(aes(pos_x, pos_y, color = eval(as.name(var))), size = 0.5) + 
            theme_bw() +
            theme(axis.text=element_text(size=text_size), 
                  axis.title=element_text(size=title_size, face='bold'),
                  legend.title=element_text(size=title_size, face='bold'), 
                  legend.text=element_text(size=text_size),
                  plot.background = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            xlab("Position X") + 
            ylab("Position Y") + 
            labs(color = legend_title) +
            coord_fixed() +
            guides(color = guide_legend(override.aes = list(size = legend_size_pt)))
            
            pdf(file = paste0(plotDir,"/", var, "_spatial.pdf"), width = 6, height = 8)
            print(gg)
            dev.off()
            print(gg)
        }
        
        else if (var == 'bin_all') {
            gg <- ggplot2::ggplot(dat_distinct) +
            geom_point(aes(pos_x, pos_y, color = eval(as.name(var))), size=0.5) + 
            theme_bw() + 
            theme(axis.text=element_text(size=text_size), 
                  axis.title=element_text(size=title_size, face='bold'),
                  legend.title=element_text(size=title_size, face='bold'),
                  legend.text=element_text(size=text_size),
                  legend.key.height = unit(legend_height_bar,"cm"),
                  plot.background = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            xlab("Position X") + 
            ylab("Position Y") + 
            labs(color = legend_title) + 
            coord_fixed() +
            scale_color_gradientn(colours = sample(rainbow(500), 500))
            
            pdf(file = paste0(plotDir,"/", var, "_spatial.pdf"), width = 6, height = 8)
            print(gg)
            dev.off()
            print(gg)
        }
        
        else {
            gg <- ggplot2::ggplot(dat_distinct) + 
            geom_point(aes(pos_x, pos_y, color = eval(as.name(var))), size = 0.6) + 
            theme_bw() + 
            theme(axis.text=element_text(size=text_size), 
                  axis.title=element_text(size=title_size, face='bold'),
                  legend.title=element_text(size=title_size, face='bold'), 
                  legend.text=element_text(size=text_size),
                  legend.key.height = unit(legend_height_bar,"cm"),
                  plot.background = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            xlab("Position X") + 
            ylab("Position Y") + 
            labs(color = legend_title) +
            coord_fixed()
            
            pdf(file = paste0(plotDir,"/", var, "_spatial.pdf"), width = 6, height = 8)
            print(gg)
            dev.off()
            print(gg)

        }        
             
    }
    
}

#' Convert to wide bin x genes + metadata format
#'
#' This function will combine beads into bins, taking the average expression intensities, average positions,  
#' most common cluster seurat cluster, and most common cluster/tissue type of constituent beads
#'
#' @param dat_long data.table of bead expression intensities per gene with metadata in long format
#' @param plotDir output plot directory path
#' @param spatial True if using spatial information
#' @return data.table of expression intensities at aggregated bin level
#' @export
#' @import ggplot2

long_to_bin <- function(dat_long, 
                        plotDir, 
                        spatial=TRUE) {
    
    # Using spatial information
    if (isTRUE(spatial)) {
        # combine constituent beads to bins
        dat_bin=dat_long[,.(value=mean(value),
                            pos_x=mean(pos_x),
                            pos_y=mean(pos_y),
                            seurat_clusters=mode(seurat_clusters), 
                            cluster_type=mode(cluster_type)),
                            by=c("bin_all","GENE","chr","start","end","rel_gene_pos", 
                                 "length", "N_bin", "umi_bin")]
        data.table::setnames(dat_bin,"bin_all","variable")
        dat_bin[,plot_order:=c(1:length(value)),by=c("variable","chr")] # plot order is order of genes on chrom.

        # Plot max seurat cluster of a given bin
        gg <- ggplot2::ggplot(unique(dat_bin[,c('pos_x','pos_y','seurat_clusters')])) +
          geom_jitter(aes(pos_x, pos_y, color = seurat_clusters), width = 0.1)

        pdf(file = paste0(plotDir,"/seurat_clusters_bin_spatial.pdf"), width = 10, height = 6)
        print(gg)
        dev.off()
        print(gg)
    }
    
    # Not using spatial information
    else {
        # combine constituent beads to bins
        dat_bin=dat_long[,.(value=mean(value),
                            labels_cl_unif=mode(labels_cl_unif), 
                            cluster_type=mode(cluster_type)), 
                            by=c("bin_all", "GENE", "chr",
                                 "start","end","rel_gene_pos",
                                 "length", "N_bin", "umi_bin")]
        
        data.table::setnames(dat_bin,"bin_all","variable")
        dat_bin[,plot_order:=c(1:length(value)),by=c("variable","chr")] # plot order is order of genes on chrom.
    }
    
    return(dat_bin)
}                               

### Subfunction for long_to_bin() that finds mode of vector/column

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Scale for nUMI (UMI Count) to generate CNV scores
#'
#' This function re-scales expression intensities to be in a smaller range, normalizes for nUMI per bin,
#' and subtracts reference bead signal
#'
#' @param dat_bin data.table of relative expression intensities per bin
#' @param thresh_hard TRUE if using strict thresholds for expression thresholds and FALSE if adjusting 
#' thresholds based on 1 + or - the mean of absolute min and max vlaues
#' @return data.table of CNV scores per bin

scale_nUMI <- function(dat_bin, 
                       thresh_hard=FALSE) {
    
    # Setting thresholds for scale range
    if (!isTRUE(thresh_hard)) {
        bot <- quantile(dat_bin$value, na.rm=TRUE)[[1]] # min
        top <- quantile(dat_bin$value, na.rm=TRUE)[[5]] # max

        thresh=mean(abs(c(bot, top)))
        bot_thresh <- 1-thresh
        top_thresh <- 1+thresh
    }
    else {
        bot_thresh=0.5
        top_thresh=1.5
    }

    for (nbin in unique(dat_bin$umi_bin)) {

            start = bot_thresh
            end = top_thresh
        
            # Capping extreme values and normalizing for nUMIs for each bin
            dat_bin[dat_bin$umi_bin==nbin,]$value <- package::scalefit(dat_bin, nbin, start, end)
        
            # Adjust for reference beads
            normal_mean <- mean(dat_bin[dat_bin$cluster_type=='Normal',]$value) # mean of reference beads
            dat_bin$value <- dat_bin$value - normal_mean + 1 
        }
    
    return(dat_bin)
}                  

### Subfunction for scale_nUMI that normalizes a given bin for UMI count and centers the  mean CNV score at 1

scalefit <- function(obj, 
                     nbin, 
                     start, 
                     end) {
    scaled <- scales::rescale(obj[obj$umi_bin==nbin,]$value, to=c(start,end)) # scale for nUMIs 
    scaled <- scale(scaled, scale=FALSE) + 1 # re-center at 1
    return(scaled)
}