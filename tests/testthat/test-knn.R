
source("./testing_helpers.R")

sc <- data.frame(row.names = c("CD19", "CD3", "CD4", "CD8", "ACTB", "CD45"),
                 "B Cell"     = c(4, 0, 0, 0, 0, 3),
                 "CD8 T Cell" = c(0, 5, 0, 8, 6, 1),
                 "CD4 T Cell" = c(0, 6, 4, 0, 7, 2),
                 "HSC"        = c(1, 2, 3, 4, 8, 1))

sc <- as.matrix(sc)

ca <- suppressWarnings(APL::cacomp(sc,
                                   princ_coords = 3,
                                   dims = 4, 
                                   ntop = nrow(sc)))

ca_dists <- calc_distances(caobj = ca)


############## Tests ############

test_that("bigraph without loops, sparse Matrix",{
  
  
  bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                            gene_dists = ca_dists[["gg"]],
                            cell_gene_assr = ca_dists[["cg"]],
                            gene_cell_assr = ca_dists[["gc"]],
                            k_c = 2,
                            k_g = 2,
                            k_cg = 2,
                            k_gc = 2,
                            loops = FALSE,
                            select_genes = TRUE,
                            prune_overlap = TRUE,
                            overlap = 0.2,
                            calc_gene_cell_kNN = FALSE,
                            marker_genes = c("CD8", "CD4"))
  
  bigraph2 = create_bigraph(caobj = ca,
                            k_c = 2,
                            k_g = 2,
                            k_cg = 2,
                            k_gc = 2,
                            loops = FALSE,
                            select_genes = TRUE,
                            prune_overlap = TRUE,
                            overlap = 0.2,
                            calc_gene_cell_kNN = FALSE,
                            marker_genes = c("CD8", "CD4"))
  
  expect_equal(bigraph1, bigraph2)
  
})


test_that("bigraph with loops, sparse Matrix",{
    
    
    bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                                     gene_dists = ca_dists[["gg"]],
                                     cell_gene_assr = ca_dists[["cg"]],
                                     gene_cell_assr = ca_dists[["gc"]],
                                     k_c = 2,
                                     k_g = 2,
                                     k_cg = 2,
                                     k_gc = 2,
                                     loops = TRUE,
                                     select_genes = TRUE,
                                     prune_overlap = TRUE,
                                     overlap = 0.2,
                                     calc_gene_cell_kNN = FALSE,
                                     marker_genes = c("CD8", "CD4"))
    
    bigraph2 = create_bigraph(caobj = ca,
                              k_c = 2,
                              k_g = 2,
                              k_cg = 2,
                              k_gc = 2,
                              loops = TRUE,
                              select_genes = TRUE,
                              prune_overlap = TRUE,
                              overlap = 0.2,
                              calc_gene_cell_kNN = FALSE,
                              marker_genes = c("CD8", "CD4"))
    
    expect_equal(bigraph1, bigraph2)
    
})


test_that("no gene selection, no loop",{
    
    
    bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                                     gene_dists = ca_dists[["gg"]],
                                     cell_gene_assr = ca_dists[["cg"]],
                                     gene_cell_assr = ca_dists[["gc"]],
                                     k_c = 2,
                                     k_g = 2,
                                     k_cg = 2,
                                     k_gc = 2,
                                     loops = TRUE,
                                     select_genes = FALSE,
                                     prune_overlap = FALSE,
                                     overlap = 0.2,
                                     calc_gene_cell_kNN = FALSE,
                                     marker_genes = c("CD8", "CD4"))
    
    bigraph2 = create_bigraph(caobj = ca,
                              k_c = 2,
                              k_g = 2,
                              k_cg = 2,
                              k_gc = 2,
                              loops = TRUE,
                              select_genes = FALSE,
                              prune_overlap = FALSE,
                              overlap = 0.2,
                              calc_gene_cell_kNN = FALSE,
                              marker_genes = c("CD8", "CD4"))
    
    expect_equal(bigraph1, bigraph2)
})



test_that("calculate gc-graph",{
    
    
    bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                                     gene_dists = ca_dists[["gg"]],
                                     cell_gene_assr = ca_dists[["cg"]],
                                     gene_cell_assr = ca_dists[["gc"]],
                                     k_c = 2,
                                     k_g = 2,
                                     k_cg = 2,
                                     k_gc = 2,
                                     loops = TRUE,
                                     select_genes = FALSE,
                                     prune_overlap = FALSE,
                                     overlap = 0.2,
                                     calc_gene_cell_kNN = TRUE,
                                     marker_genes = c("CD8", "CD4"))
    
    bigraph2 = create_bigraph(caobj = ca,
                              k_c = 2,
                              k_g = 2,
                              k_cg = 2,
                              k_gc = 2,
                              loops = TRUE,
                              select_genes = FALSE,
                              prune_overlap = FALSE,
                              overlap = 0.2,
                              calc_gene_cell_kNN = TRUE,
                              marker_genes = c("CD8", "CD4"))
    
    expect_equal(bigraph1, bigraph2)
})


test_that("no marker genes, differing number of k",{
    
    
    bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                                     gene_dists = ca_dists[["gg"]],
                                     cell_gene_assr = ca_dists[["cg"]],
                                     gene_cell_assr = ca_dists[["gc"]],
                                     k_c = 3,
                                     k_g = 3,
                                     k_cg = 2,
                                     k_gc = 1,
                                     loops = TRUE,
                                     select_genes = FALSE,
                                     prune_overlap = FALSE,
                                     overlap = 0.2,
                                     calc_gene_cell_kNN = TRUE)
    
    bigraph2 = create_bigraph(caobj = ca,
                              k_c = 3,
                              k_g = 3,
                              k_cg = 2,
                              k_gc = 1,
                              loops = TRUE,
                              select_genes = FALSE,
                              prune_overlap = FALSE,
                              overlap = 0.2,
                              calc_gene_cell_kNN = TRUE)
    
    expect_equal(bigraph1, bigraph2)
})


test_that("marker genes, differing number of k",{
    
    
    bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                                     gene_dists = ca_dists[["gg"]],
                                     cell_gene_assr = ca_dists[["cg"]],
                                     gene_cell_assr = ca_dists[["gc"]],
                                     k_c = 3,
                                     k_g = 3,
                                     k_cg = 2,
                                     k_gc = 1,
                                     loops = TRUE,
                                     select_genes = FALSE,
                                     prune_overlap = FALSE,
                                     overlap = 0.2,
                                     calc_gene_cell_kNN = FALSE,
                                     marker_genes = c("CD45"))
    
    bigraph2 = create_bigraph(caobj = ca,
                              k_c = 3,
                              k_g = 3,
                              k_cg = 2,
                              k_gc = 1,
                              loops = TRUE,
                              select_genes = FALSE,
                              prune_overlap = FALSE,
                              overlap = 0.2,
                              calc_gene_cell_kNN = FALSE,
                              marker_genes = c("CD45"))
    
    expect_equal(bigraph1, bigraph2)
})


test_that("very high overlap",{
    
    
    bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                                     gene_dists = ca_dists[["gg"]],
                                     cell_gene_assr = ca_dists[["cg"]],
                                     gene_cell_assr = ca_dists[["gc"]],
                                     k_c = 2,
                                     k_g = 2,
                                     k_cg = 2,
                                     k_gc = 2,
                                     loops = FALSE,
                                     select_genes = TRUE,
                                     prune_overlap = TRUE,
                                     overlap = 0.9,
                                     calc_gene_cell_kNN = FALSE,
                                     marker_genes = c("CD8", "CD4"))
    
    bigraph2 = create_bigraph(caobj = ca,
                              k_c = 2,
                              k_g = 2,
                              k_cg = 2,
                              k_gc = 2,
                              loops = FALSE,
                              select_genes = TRUE,
                              prune_overlap = TRUE,
                              overlap = 0.9,
                              calc_gene_cell_kNN = FALSE,
                              marker_genes = c("CD8", "CD4"))
    
    expect_equal(bigraph1, bigraph2)
    
})

test_that("very low overlap",{
    
    
    bigraph1 = create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                                     gene_dists = ca_dists[["gg"]],
                                     cell_gene_assr = ca_dists[["cg"]],
                                     gene_cell_assr = ca_dists[["gc"]],
                                     k_c = 2,
                                     k_g = 2,
                                     k_cg = 2,
                                     k_gc = 2,
                                     loops = FALSE,
                                     select_genes = TRUE,
                                     prune_overlap = TRUE,
                                     overlap = 0.05,
                                     calc_gene_cell_kNN = FALSE,
                                     marker_genes = c("CD8", "CD4"))
    
    bigraph2 = create_bigraph(caobj = ca,
                              k_c = 2,
                              k_g = 2,
                              k_cg = 2,
                              k_gc = 2,
                              loops = FALSE,
                              select_genes = TRUE,
                              prune_overlap = TRUE,
                              overlap = 0.05,
                              calc_gene_cell_kNN = FALSE,
                              marker_genes = c("CD8", "CD4"))
    
    expect_equal(bigraph1, bigraph2)
    
})