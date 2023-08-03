
sc <- readRDS("./testdata/mini_lympho_example.rds")

sc["ACTB",] <- c(0, 6, 7, 8)
sc[,"HSC"] <- c(1, 2, 3, 8, 0)

sc <- rbind(sc, "CD45" = c(3, 1, 2,1))

ca <- suppressWarnings(APL::cacomp(sc,
                                   princ_coords = 3,
                                   dims = 4, 
                                   ntop = nrow(sc)))

ca_dists <- calc_distances(caobj = ca)

test_that("bigraph without loops, sparse Matrix",{
  
  
  bigraph1 = create_bigraph(cell_dists = ca_dists[["cc"]],
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
  
  bigraph2 = create_bigraph_biocneighbors_spmat(caobj = ca,
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