
data_dir <- file.path("./testdata")
data <- readRDS("./testdata/mini_lympho_example.rds")

ca <- suppressWarnings(APL::cacomp(data, princ_coords = 3, ntop = nrow(df)))
ca_dists <- calc_distances(caobj = ca)

#####################
## test SNN graph ###
#####################

test_that("SNN graph with no loops and transposed gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_noLoops_trans_gcKNN.rds"
  ))
  
  
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = FALSE,
                    SNN_prune = 0,
                    mode = "out",
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = FALSE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
  
})

test_that("SNN graph with no loops and gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_noLoops_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = FALSE,
                    mode = "out",
                    SNN_prune = 0,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = TRUE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  SNN_seu <- readRDS(file.path(
    data_dir,
    "SNN_Seurat_prune1_15_outgoing_noLoops_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = FALSE,
                    mode = "out",
                    SNN_prune = 1/15,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = TRUE)
  
  expect_equal(SNN, SNN_seu)
  
})

test_that("SNN graph with loops and transposed gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_withLoops_trans_gcKNN.rds"
  ))

  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = TRUE,
                    mode = "out",
                    SNN_prune = 0,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = FALSE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
  
})

test_that("SNN graph with loops and gcKNN, outgoing edges only", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_outgoing_withLoops_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = TRUE,
                    mode = "out",
                    SNN_prune = 0,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = TRUE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  SNN_seu <- readRDS(file.path(
    data_dir,
    "SNN_Seurat_prune1_15_outgoing_withLoops_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = TRUE,
                    mode = "out",
                    SNN_prune = 1/15,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = TRUE)
  
  expect_equal(SNN, SNN_seu)
  
})


test_that("SNN graph with loops and gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_withLoops_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = TRUE,
                    mode = "all",
                    SNN_prune = 0,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = TRUE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
 
  
})

test_that("SNN graph no loops and gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_noLoops_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = FALSE,
                    mode = "all",
                    SNN_prune = 0,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = TRUE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  
})

test_that("SNN graph no loops and transposed gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_noLoops_trans_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = FALSE,
                    mode = "all",
                    SNN_prune = 0,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = FALSE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  
})

test_that("SNN graph with loops and transposed gcKNN, all edges", {
  
  
  SNN_igraph <- readRDS(file.path(
    data_dir,
    "SNN_igraph_all_withLoops_trans_gcKNN.rds"
  ))
  
  SNN <- create_SNN(caobj = ca, 
                    distances = ca_dists,
                    k_c = 2,
                    k_g = 2,
                    k_cg = 2,
                    k_gc = 2,
                    loops = TRUE,
                    mode = "all",
                    SNN_prune = 0,
                    select_genes = FALSE,
                    prune_overlap = FALSE,
                    overlap = 0,
                    calc_gene_cell_kNN = FALSE)
  SNN <- as.matrix(SNN)
  
  expect_equal(SNN, SNN_igraph)
  
  
})
