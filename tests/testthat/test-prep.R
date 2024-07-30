library(data.table)
library(dplyr)


### weight_rollmean() 

### Testing strategy
# dat
#    number of genes on chr = 1, < k, > k
#    genes on chr have same expression value, different expression values 
#    1, >1 chromosome
#.   1, >1 bead
# k
#    = 1, 101

# Single gene dat
# covers
# dat
#    number of genes on chr = 1
#    1 chromosome
#    1 bead
# k
#    = 1
test_that("smoothing works with one gene", {
    n_genes <- 1
    dat_one_gene <- data.frame(GENE = 'GENE1',
                                 chr = 'chr1',
                                 start = 1,
                                 end = 2,
                                 rel_gene_pos = 1,
                                 length = 1,
                                 bead1 = c(5)) %>% data.table()
    
    rm_one_gene <- SlideCNA::weight_rollmean(dat=dat_one_gene, k=1)

    # should stay same
    expect_equal(rm_one_gene, dat_one_gene)
})

# Small mock data frame w/ different expression values
n_genes <- 5
dat_less_genes <- data.frame(GENE = paste0("GENE", 1:n_genes),
                             chr = 'chr1',
                             start = c(1:n_genes),
                             end = c(2:(n_genes+1)),
                             rel_gene_pos = c(1:n_genes),
                             length = 1,
                             bead1 = c(5:1),
                             bead2 = c(1:5)) %>% data.table()

# Small dat, # of genes on chr < k
# covers
# dat
#    number of genes on chr < k
#    genes on chr have different expression values
#    > 1 bead
# k
#    = 101
test_that("smoothing works on diff expr values and # genes < k", {
    n_genes <- 5
    
    rm_less_genes <- SlideCNA::weight_rollmean(dat=dat_less_genes, k=101)
    
    rm_less_genes_exp <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = 'chr1',
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 bead1 = c(13/3, 15/4, 3, 9/4, 5/3),
                                 bead2 = c(5/3, 9/4, 3, 15/4, 13/3)) %>% data.table()

    expect_equal(rm_less_genes, rm_less_genes_exp)
})

# Small dat, # of genes on chr > k
# covers
# dat
#    number of genes on chr > k
test_that("smoothing works # genes > k", {
    n_genes <- 5
    
    rm_less_genes <- SlideCNA::weight_rollmean(dat=dat_less_genes, k=3)
    
    rm_less_genes_exp <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = 'chr1',
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 bead1 = c(14/3, 4, 3, 2, 4/3),
                                 bead2 = c(4/3, 2, 3, 4, 14/3)) %>% data.table()

    expect_equal(rm_less_genes, rm_less_genes_exp)
})

# Same expr dat
# covers
# dat
#    number of genes on chr > k
#    genes on chr have same expr value
#    > 1 chr
# k
#    = 101
test_that("smoothing works on genes w/ same expr value on chr and > 1 chr", {
    
    # 1st chromosome has genes w/ same expr value
    # 2nd chromosome has # of genes < k
    n_genes_chr1 <- 102
    n_genes_chr2 <- 5
    n_genes <- n_genes_chr1 + n_genes_chr2

    dat_same_expr <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = c(replicate(n_genes_chr1, 'chr1'),
                                         replicate(n_genes_chr2, 'chr2')),
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 # all the same expr value on chrs
                                 bead1 = c(replicate(n_genes_chr1, 2),
                                          replicate(n_genes_chr2, 0))) %>% data.table()
    
    rm_same_expr <- SlideCNA::weight_rollmean(dat=dat_same_expr, k=101)

    expect_equal(rm_same_expr, dat_same_expr)
})

### weight_rollmean() 

### Testing strategy
# normal_beads
#    number of normal beads = 1, > 1
#    normal beads have mean expr value == 0, != 0

# number of normal beads have mean expr == 0
# covers
# normal_beads
#    normal beads have mean expr value == 0
test_that("ref adj works for normal mean expr value == 0", {
    
    n_genes <- 5
    rm_zero <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = 'chr1',
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 bead_malig = c(5, 4, 3, 2, 1),
                                 bead_norm = replicate(n_genes, 0)) %>% data.table()
        
    rm_adj_zero <- SlideCNA::ref_adj(centered_rm=rm_zero, 
                                     normal_beads=c('bead_norm'))

    # should stay same
    expect_equal(rm_adj_zero, rm_zero)
})

# number of normal beads = 1, normal beads have mean expr != 0
# covers
# normal_beads
#    number of normal beads = 1
#    normal beads have mean expr value != 0
test_that("ref adj works for # of normal beads = 1 and normal expr value != 0", {
    
    n_genes <- 5
    rm_one_bead <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = 'chr1',
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 bead_malig = c(5, 4, 3, 2, 1),
                                 bead_norm = replicate(n_genes, 1)) %>% data.table()
        
    rm_adj_one_bead <- SlideCNA::ref_adj(centered_rm=rm_one_bead, 
                                         normal_beads=c('bead_norm'))
    
    rm_adj_one_bead_exp <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = 'chr1',
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 bead_malig = c(4, 3, 2, 1, 0),
                                 bead_norm = replicate(n_genes, 0)) %>% data.table()

    expect_equal(rm_adj_one_bead, rm_adj_one_bead_exp)
})

# number of normal beads > 1
# covers
# normal_beads
#    number of normal beads > 1
test_that("ref adj works for # of normal beads > 1", {
    
    n_genes <- 5
    rm_multi_beads <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = 'chr1',
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 bead_malig = c(5, 4, 3, 2, 1),
                                 bead_norm1 = replicate(n_genes, 1),
                                 bead_norm2 = replicate(n_genes, 3)) %>% data.table()
        
    rm_adj_multi_beads <- SlideCNA::ref_adj(centered_rm=rm_multi_beads, 
                                         normal_beads=c('bead_norm1', 'bead_norm2'))
    
    rm_adj_multi_beads_exp <- data.frame(GENE = paste0("GENE", 1:n_genes),
                                 chr = 'chr1',
                                 start = c(1:n_genes),
                                 end = c(2:(n_genes+1)),
                                 rel_gene_pos = c(1:n_genes),
                                 length = 1,
                                 bead_malig = c(3, 2, 1, 0, -1),
                                 bead_norm1 = replicate(n_genes, -1),
                                 bead_norm2 = replicate(n_genes, 1)) %>% data.table()

    expect_equal(rm_adj_multi_beads, rm_adj_multi_beads_exp)
})
