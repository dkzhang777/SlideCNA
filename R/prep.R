#' Infercnv-based preparation of relative gene expression intensities
#'
#' This function takes in a data table of raw counts and a vector of reference/normal beads
#' to normalize counts and adjust for reference expression. 
#'
#' @param dat data.table of raw counts
#' @param normal_beads vector of names of normal beads 
#' @param gene_pos data.table with columns for GENE, chr, start, end, rel_gene_pos (1 : # of genes on chromosome)
#' @param chrom_ord vector of the names of chromosomes in order
#' @param logTPM TRUE if performing adjustment with logTPM
#' @return A data.table of normalized, capped, and ref-adjusted counts with genomic psoition info 

#' @export
prep <- function(dat, 
                 normal_beads, 
                 gene_pos, 
                 chrom_ord, 
                 logTPM=FALSE) {
    
    dat <- as.data.frame(dat)
    normal_i <- which(colnames(dat[,-1]) %in% normal_beads)
    
    dat=merge(gene_pos,dat,by="GENE")
    
    # Normalize for gene length then sequencing depth using TPM -> log-2(TPM + 1)
    gene_pos$length <- gene_pos$end - gene_pos$start
    dat_gl <- merge(gene_pos[,c(1,6)], dat, by="GENE")
    
    if (logTPM==FALSE) {
        dat_inter <- dat_gl[,7:ncol(dat_gl)]/dat_gl$length
        dat_inter <- sweep(dat_inter,2,colSums(dat_inter),`/`)
        dat_inter <- dat_inter * 1000000
        dat_inter <- log((dat_inter + 1), 2)
        dat_inter <- sweep(dat_inter, MARGIN=1, STATS= rowMeans(dat_inter[,normal_i]))
    }
    else {
        dat_inter <- dat_gl[,7:ncol(dat_gl)]
        dat_inter <- sweep(dat_inter, MARGIN=1, STATS= rowMeans(dat_inter[,..normal_i]))
    }

    lower_bound <- mean(apply(dat_inter, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(dat_inter, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))

    thresh=mean(abs(c(lower_bound, upper_bound)))
    
    # Cap extreme values
    dat_inter[dat_inter >= 3] <- 3
    dat_inter[dat_inter <= -3] <- -3
    
    # Center by subtracting average expression                              
    dat <- cbind(GENE=dat_gl$GENE, dat_inter)
    dat=merge(gene_pos,dat,by="GENE")
    dat <- dat[dat$chr!='chrM' & dat$chr!= 'chrY',]
    
    # Order by chromosome and start position
    dat=dat[order(ordered(dat$chr, levels = chrom_ord), dat$start),]

    dat <- data.table::as.data.table(dat)
    return(dat)
}

#' Expressional smoothing along a chromosome using a weighted pyramidal moving average
#'
#' Take in a data.table of genomic positions and bead normalized/modified counts and apply pyramidal weighting 
#' with a window size k to create smoothed expression intensities
#'
#' @param dat data.table of normalized/adjusted counts
#' @param k size of window for weighting 
#' @return A data.table of expression intensities
                              
#' @export
weight_rollmean <- function(dat, 
                            k=101) {
    
    rm_mat <- matrix(nrow = 1, ncol = ncol(dat)-6)
    
    #get rolling means per chromosome
    for (chrom in unique(dat$chr)) {
        # no. of genes on chromosome larger than window size
        if (nrow(dat[dat$chr==chrom,]) > k) {
            rm_mat <- rbind(rm_mat, 
                            SlideCNA::weight_rollmean_sub(as.matrix(dat[dat$chr==chrom,7:ncol(dat)]), k=k))
        }
        # no. of genes on chromosome is just 1
        else if (nrow(dat[dat$chr==chrom,]) == 1) {
            rm_mat <- rbind(rm_mat, as.matrix(dat[chr==chrom,7:ncol(dat)]))
        }
        # no. of genes on chromosome is more than one but less than window size -> set window size that length
        else {
            rm_mat <- rbind(rm_mat, SlideCNA::weight_rollmean_sub(as.matrix(dat[dat$chr==chrom,7:ncol(dat)]), k= 
                                                        (nrow(dat[dat$chr==chrom,]))))
        }
    }
        
    rm_mat <- rm_mat[-1,]
    rm <- as.data.frame(dat)
    rm[,7:ncol(rm)] <- as.data.frame(rm_mat)
    rm <- data.table::as.data.table(rm)
    return(rm)
    
}
                   
#' Subfunction of weight_rollmean
#'
#' Take in a counts matrix and apply pyramidal weighting 
#' with a window size k to create smoothed expression intensities
#'
#' @param mat matrix of normalized/adjusted counts
#' @param k size of window for weighting 
#' @return A matrix of smoothed counts
                              
#' @export
weight_rollmean_sub <- function(mat, 
                                k) {
    new_mat <- mat
    
    # k is even
    if ((k %% 2) == 0) {
        k <- k - 1
    }
        
    side <- (k - 1) / 2
    weights <- c(1:side, side+1, side:1) / sum(c(1:side, side+1, side:1))
    
    # Smooth per bead
    for (c in 1:ncol(mat)) {
        
        # Smooth each gene using weighted average intensities of its symmetrical k neighboring genes
        for (r in 1:nrow(mat)) {
            # first row
            if (r == 1) {
                first_weights <- c((side+1):1) / sum(c((side+1):1))
                new_mat[r,c] <- sum(mat[r:(side+1),c] * first_weights)
            }
            
            # top edge
            else if (r <= side) {
                top_weights <- c((side+2-r):side, (side+1):1) / sum(c((side+2-r):side, (side+1):1))
                new_mat[r,c] <- sum(mat[1:(r+side),c] * top_weights)
            }
            
            # last row
            else if (r == nrow(mat)) {
                last_weights <- c(1:(side+1)) / sum(c(1:(side+1)))
                new_mat[r,c] <- sum(mat[(r-side):r,c] * last_weights)
            }
            
            #bottom edge
            else if (r > (nrow(mat) - side)) {
                bottom_weights <- c(1:(side+1),side:((side+1)-(nrow(mat)-r))) / 
                sum(c(1:(side+1),side:((side+1)-(nrow(mat)-r))))
                new_mat[r,c] <- sum(mat[(r-side):nrow(mat),c] * bottom_weights)
            }
            
            else{
                new_mat[r,c] <- sum(mat[(r-side):(r+side),c] * weights)
            }
        }
    }
    
    return(new_mat)
}
                              
#' Center expression intensities
#'
#' Take in a data.table of genomic positions and smoothed expression intensities counts and center by 
#' subtracting average intensity across all beads for each gene
#'
#' @param rm data.table of smoothed expression intensities counts
#' @return centered_rm data.table of smoothed, centered expression intensities
                            
#' @export
center_rm <- function(rm) {
    meta <- rm[,c(1:6)]
    rm <- rm[,c(-1,-3,-4,-5,-6)]

    # subtract average intensities of all beads to center at 0
    centered_rm <- sweep(rm[,2:ncol(rm)], 1, rowMeans(as.matrix(rm[,2:ncol(rm)])), '-')
    
    centered_rm <- cbind(meta, centered_rm)
    return(centered_rm)

}
                              
#' Adjust for Reference (Normal) Beads
#'
#' Take in a data.table of genomic positions and smoothed, centered expression intensities counts and
#' adjust for reference beads by subtracting average intensities of reference beads for each gene.
#' This is the second reference adjustment.                             
#'
#' @param centered_rm data.table of smoothed, centered expression intensities counts
#' @param normal_beads vector of names of normal beads                               
#' @return rm_adj data.table of smoothed relative expression intensities
                              
#' @export
ref_adj <- function(centered_rm, 
                    normal_beads) {
    rm_adj <- as.data.frame(centered_rm[,-c(1:6)])
    normal_mean <- rowMeans(rm_adj[,normal_beads]) # mean of reference beads
    
    # Subtract mean of normal beads from all beads for each gene
    rm_adj <- rm_adj - normal_mean
    
    rm_adj <- cbind(centered_rm[,c(1:6)], rm_adj)
    rm_adj <- data.table::as.data.table(rm_adj)
    return(rm_adj)
}                              
