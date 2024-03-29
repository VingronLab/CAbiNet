source("./tests/testthat/testing_helpers.R")

library(APL)
library(Matrix)
# df <- data.frame(row.names = c("CD19", "CD3", "CD4", "CD8", "ACTB"),
#                  "B Cell"     = c(4, 0, 0, 0, 0),
#                  "CD8 T Cell" = c(0, 5, 0, 8, 8),
#                  "CD4 T Cell" = c(0, 6, 4, 0, 8),
#                  "HSC"       = c(2, 2, 2, 2, 8))


df <- data.frame(row.names = c("CD19", "CD3", "CD4", "CD8", "ACTB", "CD45"),
                   "B Cell"     = c(4, 0, 0, 0, 0, 3),
                   "CD8 T Cell" = c(0, 5, 0, 8, 6, 1),
                   "CD4 T Cell" = c(0, 6, 4, 0, 7, 2),
                   "HSC"        = c(1, 2, 3, 4, 8, 1))
df <- as.matrix(df)

dists = as.matrix(dist(df))
saveRDS(dists, file.path("./tests/testthat/testdata/testdata_Euclidean_dist.rds"))

saveRDS(df, "./tests/testthat/testdata/mini_lympho_example.rds")

ca <- suppressWarnings(APL::cacomp(df, princ_coords = 3, ntop = nrow(df)))

ca_dists <- calc_distances(caobj = ca)

#####################
#####################
# no loops 

adj <- create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                      gene_dists = ca_dists[["gg"]],
                      cell_gene_assr = ca_dists[["cg"]],
                      gene_cell_assr = ca_dists[["gc"]],
                      k_c = 2,
                      k_g = 2,
                      k_cg = 2,
                      k_gc = 2,
                      loops = FALSE,
                      overlap = 0,
                      prune_overlap = FALSE,
                      select_genes = FALSE,
                      calc_gene_cell_kNN = TRUE)

saveRDS(adj, "./tests/testthat/testdata/testdata_gcKNN_noLoops_adj.rds")
adj <- readRDS("./tests/testthat/testdata/testdata_gcKNN_noLoops_adj.rds")

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "out"
)

rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)

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

rownames(snn_matrix_seu) <- rownames(adj)
colnames(snn_matrix_seu) <- colnames(adj)
saveRDS(snn_matrix_seu,
        "./tests/testthat/testdata/SNN_Seurat_prune1_15_outgoing_noLoops_gcKNN.rds")

## no loops in & outgoing

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "all"
)
rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_all_noLoops_gcKNN.rds")


## no loops incoming
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "in"
)
rownames(snn_igraph) = rownames(adj)
colnames(snn_igraph) = colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_incoming_noLoops_gcKNN.rds")

#####################
#####################
# no loops transpose

adj <- create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                      gene_dists = ca_dists[["gg"]],
                      cell_gene_assr = ca_dists[["cg"]],
                      gene_cell_assr = ca_dists[["gc"]],
                      k_c = 2,
                      k_g = 2,
                      k_cg = 2,
                      k_gc = 2,
                      loops = FALSE,
                      overlap = 0,
                      prune_overlap = FALSE,
                      select_genes = FALSE,
                      calc_gene_cell_kNN = FALSE)

saveRDS(adj, "./tests/testthat/testdata/testdata_trans_gcKNN_noLoops_adj.rds")
adj <- readRDS("./tests/testthat/testdata/testdata_trans_gcKNN_noLoops_adj.rds")

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "out"
)
rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)

saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_noLoops_trans_gcKNN.rds")


## no loops in & outgoing transpose
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "all"
)
rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_all_noLoops_trans_gcKNN.rds")

## no loops incoming transpose
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = FALSE),
  method = c("jaccard"),
  loops = FALSE,
  mode = "in"
)
rownames(snn_igraph) = rownames(adj)
colnames(snn_igraph) = colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_incoming_noLoops_trans_gcKNN.rds")

#####################
#####################
# with loops


adj <- create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                      gene_dists = ca_dists[["gg"]],
                      cell_gene_assr = ca_dists[["cg"]],
                      gene_cell_assr = ca_dists[["gc"]],
                      k_c = 2,
                      k_g = 2,
                      k_cg = 2,
                      k_gc = 2,
                      loops = TRUE,
                      overlap = 0,
                      prune_overlap = FALSE,
                      select_genes = FALSE,
                      calc_gene_cell_kNN = TRUE)

saveRDS(adj, "./tests/testthat/testdata/testdata_gcKNN_withLoops_adj.rds")
adj <- readRDS("./tests/testthat/testdata/testdata_gcKNN_withLoops_adj.rds")

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "out"
)
rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_withLoops_gcKNN.rds")


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
rownames(snn_matrix_seu) <- rownames(adj)
colnames(snn_matrix_seu) <- colnames(adj)
saveRDS(snn_matrix_seu,
        "./tests/testthat/testdata/SNN_Seurat_prune1_15_outgoing_withLoops_gcKNN.rds")

## with loops in & outgoing

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "all"
)
rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_all_withLoops_gcKNN.rds")

## with loops incoming
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "in"
)
rownames(snn_igraph) = rownames(adj)
colnames(snn_igraph) = colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_incoming_withLoops_gcKNN.rds")


#####################
#####################
# with loops transpose 


adj <- create_bigraph_manual(cell_dists = ca_dists[["cc"]],
                      gene_dists = ca_dists[["gg"]],
                      cell_gene_assr = ca_dists[["cg"]],
                      gene_cell_assr = ca_dists[["gc"]],
                      k_c = 2,
                      k_g = 2,
                      k_cg = 2,
                      k_gc = 2,
                      loops = TRUE,
                      overlap = 0,
                      prune_overlap = FALSE,
                      select_genes = FALSE,
                      calc_gene_cell_kNN = FALSE)

saveRDS(adj, "./tests/testthat/testdata/testdata_trans_gcKNN_withLoops_adj.rds")
adj <- readRDS("./tests/testthat/testdata/testdata_trans_gcKNN_withLoops_adj.rds")

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "out"
)
rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_outgoing_withLoops_trans_gcKNN.rds")


## with loops in & outgoing transpose 

snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "all"
)
rownames(snn_igraph) <- rownames(adj)
colnames(snn_igraph) <- colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_all_withLoops_trans_gcKNN.rds")

## with loops incoming transpose
snn_igraph <- igraph::similarity(
  igraph::graph_from_adjacency_matrix(adj, diag = TRUE),
  method = c("jaccard"),
  loops = TRUE,
  mode = "in"
)
rownames(snn_igraph) = rownames(adj)
colnames(snn_igraph) = colnames(adj)
saveRDS(snn_igraph, "./tests/testthat/testdata/SNN_igraph_incoming_withLoops_trans_gcKNN.rds")

####################
