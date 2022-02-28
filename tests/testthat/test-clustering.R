
library(APL)
df <- data.frame(row.names = c("CD19", "CD3", "CD4", "CD8", "ACTB"),
                 "B Cell"     = c(4, 0, 0, 0, 0),
                 "CD8 T Cell" = c(0, 5, 0, 8, 8),
                 "CD4 T Cell" = c(0, 6, 4, 0, 8),
                 # "Platelet"   = c(0, 0, 0, 0, 5, 3),
                 "HSC"       = c(2, 2, 2, 2, 8))
df <- as.matrix(df)


jaccard <- function(a, b, common_edges) {
  union = sum(a) + sum(b) - common_edges
  return (common_edges/union)
}

ca <- suppressWarnings(APL::cacomp(df, princ_coords = 3, ntop = nrow(df)))

ca_dists <- calc_distances(caobj = ca)

adj <- create_bigraph(cell_dists = ca_dists[["cc"]],
                      gene_dists = ca_dists[["gg"]],
                      cell_gene_assr = ca_dists[["cg"]],
                      gene_cell_assr = ca_dists[["gc"]],
                      k_c = 2,
                      k_g = 2,
                      k_cg = 2,
                      k_gc = 2,
                      overlap = 0,
                      prune_overlap = FALSE,
                      select_genes = FALSE,
                      calc_gene_cell_kNN = TRUE)
saveRDS(adj, "./tests/testthat/testdata/handmade_gcKNN_adj.rds")

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "out"
)

saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_withLoops_gcKNN.rds")

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "out"
)

saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_noLoops_gcKNN.rds")



stopifnot(length(unique(rowSums(adj)))==1)
k.param <- sum(adj[1,])

GSG.idx <- matrix(NA, nrow = nrow(adj), ncol = k.param)
rownames(GSG.idx) <- rownames(adj)

for (r in 1:nrow(adj)){
  GSG.idx[r,] <- which(adj[r,] == 1)
}

snn_matrix_seu <- Seurat:::ComputeSNN(
  nn_ranked = GSG.idx,
  prune = 1/15)

saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_noLoops_transpose_gcKNN.rds")



stopifnot(length(unique(rowSums(adj)))==1)
k.param <- sum(adj[1,])

GSG.idx <- matrix(NA, nrow = nrow(adj), ncol = k.param)
rownames(GSG.idx) <- rownames(adj)

for (r in 1:nrow(adj)){
  GSG.idx[r,] <- which(adj[r,] == 1)
}

snn_matrix_seu <- Seurat:::ComputeSNN(
  nn_ranked = GSG.idx,
  prune = 1/15)

saveRDS(snn_matrix_seu, "./tests/testthat/testdata/SNN_seurat_outgoing_nodiag_gcKNN.rds")

####################

adj <- create_bigraph(cell_dists = ca_dists[["cc"]],
                      gene_dists = ca_dists[["gg"]],
                      cell_gene_assr = ca_dists[["cg"]],
                      gene_cell_assr = ca_dists[["gc"]],
                      k_c = 2,
                      k_g = 2,
                      k_cg = 2,
                      k_gc = 2,
                      overlap = 0,
                      prune_overlap = FALSE,
                      select_genes = FALSE,
                      calc_gene_cell_kNN = FALSE)

saveRDS(adj, "./tests/testthat/testdata/handmade_gcKNN_transpose_adj.rds")


snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "out"
)

saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_withLoops_transpose_gcKNN.rds")

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "out"
)


saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_noLoops_transpose_gcKNN.rds")

######################
######################

gcKNN_adj <- readRDS("./tests/testthat/testdata/handmade_gcKNN_adj.rds") %>%
gcKNN_transpose_adj <- readRDS("./tests/testthat/testdata/handmade_gcKNN_transpose_adj.rds")

adj <- as.matrix(adj)
ncommon <- adj %*% t(adj)
stopifnot(isSymmetric(ncommon))


stopifnot(nrow(adj) == ncol(adj))
snn_to_test <- matrix(NA, nrow = nrow(adj), ncol = nrow(adj))

for(i in seq_len(nrow(adj))){
  for (j in seq_len(nrow(adj))){
    snn_to_test[i,j] <- jaccard(adj[i,], adj[j,], ncommon[i,j])
  }
}
##################
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj,diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "out"
)
rownames(snn_igraph) <- gcKNN_adj
readr::write_csv(as.data.frame(snn_igraph), "./tests/testthat/testdata/handmade_gcKNN_transpose_adj.csv")

snn.matrix <- ComputeSNNjaccard( as(adj, "dgCMatrix"), 0, mode = "out")
stopifnot(sum(snn_igraph != snn.matrix) == 0)

#################
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj,diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "in"
)

snn.matrix <- ComputeSNNjaccard( as(adj, "dgCMatrix"), 0, mode = "in")
stopifnot(sum(snn_igraph != snn.matrix) == 0)
#################
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj,diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "all"
)
snn.testall <- ComputeSNNjaccard( as(adj, "dgCMatrix"), 0, mode = "all")
snn.testall = as.matrix(snn.testall)
# identical(snn_igraph, snn.testall)
stopifnot(sum(snn_igraph != snn.testall) == 0)
#################
stopifnot(length(unique(rowSums(adj)))==1)
k.param <- sum(adj[1,])

GSG.idx <- matrix(NA, nrow = nrow(adj), ncol = k.param)
rownames(GSG.idx) <- rownames(adj)

for (r in 1:nrow(adj)){
  GSG.idx[r,] <- which(adj[r,] == 1)
}


snn.matrix.seu <- Seurat:::ComputeSNN(
  nn_ranked = GSG.idx,
  prune = 0)

identical(snn.matrix, snn.matrix.seu)
# if(!is(adj, "dgCMatrix")){
#   adj <- as(adj, "dgCMatrix")
# }

View(as.matrix(adj))

isSymmetric(as.matrix(adj))



# if(!is(adj, "dgCMatrix")){
#   adj <- as(adj, "dgCMatrix")
# }

View(as.matrix(adj))

isSymmetric(as.matrix(adj))



rownames(snn.matrix) <- rownames(adj)
colnames(snn.matrix) <- rownames(adj)

SNN <- create_SNN(caobj = ca,
                  distances = ca_dists,
                  k_c = 2,
                  k_g = 2,
                  k_cg = 2,
                  k_gc = 2,
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = FALSE)

View(as.matrix(SNN))



adj_handmade <-
test_that("Correct symmetric SNN graph", {
  expect_equal(2 * 2, 4)
})


test_that("Correct assymetric SNN graph", {
  expect_equal(2 * 2, 4)
})
