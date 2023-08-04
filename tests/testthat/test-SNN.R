# data_dir <- file.path("./tests/testthat/testdata")
# data <- readRDS("./tests/testthat/testdata/mini_lympho_example.rds")
data_dir <- file.path("./testdata")
# data <- readRDS("./testdata/mini_lympho_example.rds")

data <- data.frame(row.names = c("CD19", "CD3", "CD4", "CD8", "ACTB", "CD45"),
           "B Cell"     = c(4, 0, 0, 0, 0, 3),
           "CD8 T Cell" = c(0, 5, 0, 8, 6, 1),
           "CD4 T Cell" = c(0, 6, 4, 0, 7, 2),
           "HSC"        = c(1, 2, 3, 4, 8, 1))

data <- as.matrix(data)

ca <- suppressWarnings(APL::cacomp(data, princ_coords = 3,dims = 4, ntop = nrow(data)))
ca_dists <- calc_distances(caobj = ca)

#####################
## test SNN graph ###
#####################

test_that("SNN graph with no loops and transposed gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_noLoops_trans_gcKNN.rds"
  ))
  
  
  
  caclust <- make_SNN(caobj = ca, 
                    k = 2,
                    loops = FALSE,
                    SNN_prune = 0,
                    mode = "out",
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = FALSE)
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
})

test_that("SNN graph with no loops and gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_noLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = FALSE,
                  mode = "out",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  SNN_seu <- readRDS(file.path(
    data_dir,
    "SNN_Seurat_prune1_15_outgoing_noLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = FALSE,
                  mode = "out",
                  SNN_prune = 1/15,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  
  expect_equal(caclust@SNN, SNN_seu)
  
})

test_that("SNN graph with loops and transposed gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_withLoops_trans_gcKNN.rds"
  ))

  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = TRUE,
                  mode = "out",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = FALSE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
})

test_that("SNN graph with loops and gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_withLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = TRUE,
                  mode = "out",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  SNN_seu <- readRDS(file.path(
    data_dir,
    "SNN_Seurat_prune1_15_outgoing_withLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = TRUE,
                  mode = "out",
                  SNN_prune = 1/15,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  
  expect_equal(caclust@SNN, SNN_seu)
  
})



test_that("SNN graph with loops and gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_withLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = TRUE,
                  mode = "all",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
 
  
})

test_that("SNN graph no loops and gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_noLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = FALSE,
                  mode = "all",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  
})

test_that("SNN graph no loops and transposed gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_noLoops_trans_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = FALSE,
                  mode = "all",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = FALSE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  
})

test_that("SNN graph with loops and transposed gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_withLoops_trans_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = TRUE,
                  mode = "all",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = FALSE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  
})

test_that("SNN graph with no loops and transposed gcKNN, incoming edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_incoming_noLoops_trans_gcKNN.rds"
  ))
  
  
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = FALSE,
                  SNN_prune = 0,
                  mode = "in",
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = FALSE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
})

test_that("SNN graph with no loops and gcKNN, incoming edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_incoming_noLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = FALSE,
                  mode = "in",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  
})

test_that("SNN graph with loops and transposed gcKNN, incoming edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_incoming_withLoops_trans_gcKNN.rds"
  ))
  
  
  caclust <- make_SNN(caobj = ca, 
                      k = 2,
                      loops = TRUE,
                      mode = "in",
                      SNN_prune = 0,
                      select_genes = FALSE,
                      prune_overlap = FALSE,
                      overlap = 0,
                      calc_gene_cell_kNN = FALSE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
})

test_that("SNN graph with loops and gcKNN, incoming edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_incoming_withLoops_gcKNN.rds"
  ))
  
  caclust <- make_SNN(caobj = ca, 
                  k = 2,
                  loops = TRUE,
                  mode = "in",
                  SNN_prune = 0,
                  select_genes = FALSE,
                  prune_overlap = FALSE,
                  overlap = 0,
                  calc_gene_cell_kNN = TRUE)
  
  SNN <- as.matrix(caclust@SNN)
  
  expect_equal(SNN, SNN_igraph)
  
})



