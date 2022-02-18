

#' Calculate euclidean distance
#'
#' @description Calculates euclidean distance.
#'
#' @param mat matrix.
#' @return Matrix with pairwise euclidean distances.
#'
calc_euclidean <- function(mat) {
  n = nrow(mat)
  smat <- apply(mat, 1, crossprod)
  mat1 <- matrix(smat, nrow = n, ncol = n)
  mat3 <- tcrossprod(mat)
  mat4 <- mat1 + t(mat1) - 2 * mat3
  diag(mat4) <- 0
  mat5 <- sqrt(mat4)
  mat5[is.na(mat5)] = 0
  return(mat5)
}


#' Calculate Association Ratio
#'
#' @param caobj A cacomp object with principal and standard coordinates calculated.
#' @param direction Either "cells" or "genes". If "cells" the association ratio
#' of each gene for each cell is calculated. If "genes" the assoc. ratio of each
#' cell for each gene is calculated.
#'
#' @return
#' Matrix with pairwise association ratios.
calc_assR <- function(caobj, direction = "cells"){

  stopifnot(is(caobj, "cacomp"))

  if (direction == "cells") {
    assR <- ca@std_coords_cols %*% t(ca@prin_coords_rows)
  } else if (direction == "genes"){
    assR <- ca@std_coords_rows %*% t(ca@prin_coords_cols)
  } else {
    stop("Please choose a valid direction. Either 'cells' or 'genes'")
  }

  return(assR)

}


#' Calculate distances for bigraph
#' 
#' @param caobj  A cacomp object with principal and standard coordinates calculated.
#' 
#' @return 
#' List with elements:
#' * "cc": cell-cell euclidean distances
#' * "gg": gene-gene euclidean distances
#' * "cg": cell-gene association ratio
#' * "gc": gene-cell association ratio
#' 
#' @md
calc_distances <- function(caobj){
  
  cell_dists <- calc_euclidean(caobj@prin_coords_cols)
  gene_dists <- calc_euclidean(caobj@prin_coords_rows)
  cell_gene_assr <- calc_assR(caobj, direction = "cells")
  gene_cell_assr <- calc_assR(caobj, direction = "genes")
  
  return(list("cc" = cell_dists,
              "gg" = gene_dists,
              "cg" = cell_gene_assr,
              "gc" = gene_cell_assr))
}


#' Make a kNN graph
#' @description Given a distance matrix `make_knn()` outputs the adjacency
#' matrix of the k-nearest-neighbours graph.
#'
#' @param dists distance matrix
#' @param k numeric. Number of nearest neighbours.
#' @param decr boolean. Whether the the values in `dists` should be sorted into
#' decreasing (TRUE) or increasing (FALSE) order.
#'
#' @return
#' adjacency matrix of kNN graph: sparse matrix of type "dgCMatrix".
#'
#' @export
make_knn <- function(dists,
                     k,
                     decr = TRUE) {

  n.row <- nrow(dists)
  n.col <- ncol(dists)

  if (n.col < k) {
    warning(
      "k set larger than number of genes Setting k to number of genes - 1.",
      call. = FALSE
    )
    k <- n.col - 1
  }
  knn.mat <- matrix(data = 0, ncol = k, nrow = n.row)
  # knd.mat <- knn.mat
  for (i in 1:n.row) {
    knn.mat[i, ] <- order(dists[i, ], decreasing = decr)[1:k]
    # knd.mat[i, ] <- dists[i, knn.mat[i, ]]
  }
  nn.ranked <- knn.mat[, 1:k]

  # convert nn.ranked into a Graph
  j <- as.numeric(t(nn.ranked))
  i <- ((1:length(j)) - 1) %/% k + 1

  nn.matrix <- sparseMatrix(i = i,
                            j = j,
                            x = 1,
                            dims = c(nrow(x = dists), ncol(x = dists)))
  rownames(nn.matrix) <- rownames(dists)
  colnames(nn.matrix) <- colnames(dists)
  return(nn.matrix)
}

#' Calculates the overlap of chosen genes between neighbouring cells
#'
#' @description
#' For the cell indicated with `idx` in the rows of `cc_adj` and `cg_adj`,
#' the genes and cells with and edge to the chosen cell are determined.
#' Among the connected it is then determined how many have an edge to the same
#' genes as the chosen cell. The overlap (cells with edge to same gene)/
#' (Nr. of cells) is returned.
#'
#' @param idx indices of cell of interest
#' @param cg_adj cell-gene adjacency matrix
#' @param cc_adj cell-cell adjacency matrix
#'
#' @return
#' numeric vector indicating the overlap between neighbours for each gene that
#' is a neighbour in `cg_adj`. Same length as number of nearest neighbour genes
#' in `cg_adj` for each cell.
#'
get_per_gene_overlap <- function(idx, cg_adj, cc_adj){

  stopifnot(length(idx) == 1)

  genes_picked <- which(cg_adj[idx,] == 1)
  neighbours <- which(cc_adj[idx,] == 1)

  cg_adj_nns <- cg_adj[neighbours, genes_picked]
  overlap <- Matrix::colSums(cg_adj_nns)/length(neighbours)


  return(overlap)
}


#' Determines the overlap of chosen genes between nearest neighbour cells
#' @description
#' Determines how many genes are shared (have an edge) between
#' neighbouring cells for all cells in the dataset.
#'
#' @param cg_adj cell-gene adjacency matrix
#' @param cc_adj cell-cell adjacency matrix
#'
#' @return
#' Matrix (cells x genes) showing the overlap of each gene between the nearest
#' neighbour cells.
#'
#' @export
determine_overlap <- function(cg_adj, cc_adj){

  overlap_mat <- matrix(0, nrow = nrow(cg_adj), ncol = ncol(cg_adj))

  for (cell in seq_len(nrow(cg_adj))){

    genes_picked <- which(cg_adj[cell,] == 1)
    overlap <- get_per_gene_overlap(idx = cell,
                                    cg_adj = cg_adj,
                                    cc_adj = cc_adj)

    overlap_mat[cell,genes_picked] <- overlap
  }

  rownames(overlap_mat) <- rownames(cg_adj)
  colnames(overlap_mat) <- colnames(cg_adj)

  return(overlap_mat)
}



#' Combine kNN graphs to large cell-gene adjecency matrix
#'
#' @description
#' TODO
#' 
#' @param cell_dists cell-cell euclidean distances
#' @param gene_dists gene-gene euclidean distances
#' @param cell_gene_assr cell-gene association ratios
#' @param gene_cell_assr gene-cell association ratios
#' @param k_c k for cell-cell kNN
#' @param k_g k for gene-gene kNN
#' @param k_cg k for cell-gene kNN
#' @param k_gc k for gene-cell kNN
#' @param select_genes TRUE/FALSE. Should genes be selected by wether they have
#' an edge in the cell-gene kNN graph?
#' @param prune_overlap TRUE/FALSE. If TRUE edges to genes that share less
#' than `overlap` of genes with the nearest neighbours of the cell are removed.
#' @param overlap Numeric between 0 and 1. Overlap cutoff if
#' prune_overlap = TRUE.
#' @param calc_cell_gene_kNN TRUE/FALSE. If TRUE a cell-gene graph is calculated
#' by choosing the `k_gc` nearest cells for each gene. If FALSE the cell-gene
#' graph is transposed.
#'
#' @return
#' "Bigraph" of type `dgCMatrix`. The combined adjacency matrix consists of the
#' cell-cell graph, gene-gene graph and cell-gene/gene-cell graph.
#'
#' @export
create_bigraph <- function(cell_dists,
                           gene_dists,
                           cell_gene_assr,
                           gene_cell_assr,
                           k_c,
                           k_g,
                           k_cg,
                           k_gc,
                           select_genes = TRUE,
                           prune_overlap = TRUE,
                           overlap = 0.2,
                           calc_gene_cell_kNN = FALSE) {

  cgg_nn <- make_knn(cell_gene_assr,
                     k = k_cg,
                     decr = TRUE)

  if(isTRUE(select_genes)){
    
    idx <- Matrix::colSums(cgg_nn) > 0
    cgg_nn <- cgg_nn[,idx]
    gene_dists <- gene_dists[idx,idx]
    gene_cell_assr <- gene_cell_assr[idx,]
  }


  ccg_nn <- make_knn(cell_dists,
                     k = k_c,
                     decr = FALSE)


  if(isTRUE(prune_overlap)){

    overlap_mat <- determine_overlap(cg_adj = cgg_nn,
                                     cc_adj = ccg_nn)

    cgg_nn[overlap_mat < overlap] <- 0

    idx <- Matrix::colSums(cgg_nn) > 0
    cgg_nn <- cgg_nn[,idx]
    gene_dists <- gene_dists[idx,idx]
    gene_cell_assr <- gene_cell_assr[idx,]

  }



  ggg_nn <- make_knn(gene_dists,
                     k = k_g,
                     decr = FALSE)

  if(isFALSE(calc_gene_cell_kNN)){
    gcg_nn <- Matrix::t(cgg_nn)

  } else if(isTRUE(calc_gene_cell_kNN)){
    gcg_nn <- make_knn(gene_cell_assr,
                       k = k_gc,
                       decr = TRUE)
  } else {
    stop("calc_cell_gene_kNN has to be either TRUE or FALSE!")
  }



  GSG_1 <- cbind(ccg_nn, cgg_nn)
  GSG_2 <- cbind(gcg_nn, ggg_nn)

  GSG <- rbind(GSG_1, GSG_2)
  return(GSG)
}

#' Create SNN graph from caobj
#'
#' @description 
#' TODO
#' 
#' @param caobj A cacomp object with standard and principal coordinates 
#' calculated.
#' @param distances A list containing the cell-cell, gene-gene, cell-gene and 
#' gene-cell distances/association ratios. Must be named list as follows:
#' * "cc": cell-cell euclidean distances
#' * "gg": gene-gene euclidean distances
#' * "cg": cell-gene association ratio
#' * "gc": gene-cell association ratio
#' @param k_c k for cell-cell kNN
#' @param k_g k for gene-gene kNN
#' @param k_cg k for cell-gene kNN
#' @param k_gc k for gene-cell kNN
#' @param select_genes TRUE/FALSE. Should genes be selected by wether they have
#' an edge in the cell-gene kNN graph?
#' @param prune_overlap TRUE/FALSE. If TRUE edges to genes that share less
#' than `overlap` of genes with the nearest neighbours of the cell are removed.
#' @param overlap Numeric between 0 and 1. Overlap cutoff if
#' prune_overlap = TRUE.
#' @param calc_cell_gene_kNN TRUE/FALSE. If TRUE a cell-gene graph is calculated
#' by choosing the `k_gc` nearest cells for each gene. If FALSE the cell-gene
#' graph is transposed.
#' 
#' @returns 
#' A sparse adjacency Matrix of type "dgCMatrix". The values in the matrix
#' are the Jaccard similarity between nodes in the graph. The range between 0
#' and 1, with 0 meaning that no edges are shared between nodes, wheras 1 means 
#' all edges are shared between nodes.
#' 
#' @md 
#' @export
create_SNN <- function(caobj,
                       distances,
                       k_c,
                       k_g,
                       k_cg,
                       k_gc,
                       select_genes = TRUE,
                       prune_overlap = TRUE,
                       overlap = 0.2,
                       calc_gene_cell_kNN = FALSE) {
  
  
  stopifnot(all(c("cc", "gg", "cg", "gc") %in% names(distances)))
  
  adj <- create_bigraph(cell_dists = distances[["cc"]],
                        gene_dists = distances[["gg"]],
                        cell_gene_assr = distances[["cg"]],
                        gene_cell_assr = distances[["gc"]],
                        k_c = k_c,
                        k_g = k_g,
                        k_cg = k_cg,
                        k_gc = k_gc,
                        overlap = overlap,
                        prune_overlap = prune_overlap,
                        select_genes = select_genes,
                        calc_cell_gene_kNN = calc_cell_gene_kNN)
  
  if(!is(adj, "dgCMatrix")){
    adj <- as(adj, "dgCMatrix")  
  }
  
  snn.matrix <- ComputeSNNasym(adj, prune_cutoff)
  
  rownames(snn.matrix) <- rownames(adj)
  colnames(snn.matrix) <- rownames(adj)
  return(snn.matrix)
}

#' Leiden clustering on bigraph
#' 
#' @description 
#' TODO
#' 
#' @param SNN dense or sparse matrix. A SNN graph.
#' @param resolution resolution for leiden algorithm.
#' @param n.int Number of iterations for leiden algorithm.
#' @param seed Random seed.
#' 
#' @return 
#' vector of type `factor.` Assigned clusters of cells and genes. 
#' The names of cells and genes are saved in the names of the 
#' vector. 
#' 
run_leiden <- function(SNN, 
                       resolution = 1,
                       n.int = 10, 
                       rand_seed = 2358) {
  
  clusters <- leiden(object = SNN,
                     resolution_parameter = resolution,
                     partition_type = "RBConfigurationVertexPartition",
                     initial_membership = NULL,
                     weights = NULL,
                     node_sizes = NULL,
                     n_iterations = n.int,
                     seed = rand_seed)
  
  
  clusters <- as.factor(clusters)
  names(clusters) <- rownames(SNN)
  return(clusters)
}

#' run spectral clustering
run_spectral <- function() {
  #TODO
}

#' Run biclustering
#' 
#' @description 
#' TODO
#' 
#' @param caobj A cacomp object with principal and standard coordinates 
#' calculated.
#' @param k_sym k for cell-cell and gene-gene kNN graph
#' @param k_asym k for cell-gene and gene-cell kNN graph
#' @param algorithm Algorithm for clustering. Options are "leiden" or "spectral".
#' @param return_umap TRUE/FALSE. Whether umap coordinates should be returned.
#' @param k_umap k for UMAP.
#' @inheritParams create_bigraph
#' @inheritParams create_SNN
#' @inheritParams run_leiden
#' 
#' @return
#' Returns list:
#' * 'cell_clusters': factor vector containing assigned clusters for the cell 
#' specified in the elements name.
#' * 'gene_clusters': factor vector containing assigned clusters for the genes
#' specified in the elements name.
#' * 'umap_coords': If `return_umap = TRUE` additionally a data frame with 
#' the umap_coordinates of cells and genes is returned.
#'  
#' 
#' @md
#' @export
run_caclust <- function(caobj,
                        k_sym,
                        k_asym = k_sym,
                        algorithm = "leiden",
                        select_genes = TRUE,
                        prune_overlap = TRUE,
                        overlap = 0.2,
                        calc_gene_cell_kNN = FALSE,
                        resolution = 1,
                        n.int = 10,
                        rand_seed = 2358,
                        return_umap = TRUE,
                        k_umap = k_sym) {
  
  distances <- calc_distances(caobj = caobj)
  
  SNN <- create_SNN(caobj = caobj, 
                    distances = distances,
                    k_c = k_sym,
                    k_g = k_sym,
                    k_cg = k_asym,
                    k_gc = k_asym,
                    select_genes = select_genes,
                    prune_overlap = prune_overlap,
                    overlap = overlap,
                    calc_gene_cell_kNN = calc_gene_cell_kNN)
  
  if (algorithm == "leiden"){
    
    clusters <- run_leiden(SNN = SNN, 
                           resolution = resolution,
                           n.int = n.int, 
                           rand_seed = rand_seed)
    
  } else if (algorithm == "spectral"){
    stop("Spectral clustering not yet implemented. sorry :(")
  }

  
  cell_clusters <- clusters[1:nrow(ca@prin_coords_cols)]
  gene_clusters <- clusters[-(1:nrow(ca@prin_coords_cols))]
  
  if (isTRUE(return_umap)) {
    
    umap_out <- run_biUMAP_leiden(caobj = caobj,
                                  SNN = SNN,
                                  k_umap = ,
                                  rand_seed = rand_seed)
    
    return(list("cell_clusters" = cell_clusters,
                "gene_clusters" = gene_clusters,
                "umap_coords" = umap_out))
  } 
  
  
  return(list("cell_clusters" = cell_clusters,
              "gene_clusters" = gene_clusters))
}


#' Plots UMAP depicting both cells and genes.
#' 
#' @description 
#' TODO
#' 
#' @param caobj A cacomp object with principal and standard coordinates 
#' calculated.
#' @param SNN dense or sparse matrix. A SNN graph.
#' @param k_umap Number of nearest neighbours to use from the SNN graph for
#' UMAP
#' @param rand_seed Random seed for UMAP.
#' 
#' @return 
#' data frame containing the UMAP coordinates of cells and genes.
#' 
#' @export 
run_biUMAP_leiden <- function(caobj, SNN, k_umap, rand_seed = 2358){
  
  k_snn = ncol(SNN)
  SNN_idx <- matrix(data = 0, ncol = k_snn, nrow = nrow(SNN))
  
  for (i in seq_len(nrow(SNN))){
    SNN_idx[i,] <- order(SNN[i,], decreasing = TRUE)
  }
  
  SNN_jacc <- matrix(as.matrix(SNN)[SNN_idx],
                     nrow = nrow(SNN_idx),
                     ncol = ncol(SNN_idx))
  
  
  rownames(SNN_idx) <- rownames(SNN)
  rownames(SNN_jacc) <- rownames(SNN)
  
  custom.config = umap.defaults
  custom.config$random_state = rand_seed
  
  snn_umap_graph = umap::umap.knn(indexes = SNN_idx,
                                  distances = SNN_jac)
  
  assym <- rbind(caobj@std_coords_cols, cobja@prin_coords_rows)
  assym <- assym[rownames(assym) %in% rownames(SNN),]
  
  caclust_umap = umap::umap(assym,
                            config = custom.config,
                            n_neighbors = k_umap, 
                            knn = snn_umap_graph)
  
  umap_coords <- as.data.frame(caclust_umap$layout)
  colnames(umap_coords) <- c("x", "y")
  umap_coords$type <- "gene"
  umap_coords$type[1:nrow(ca@std_coords_cols)] <- "cell"
  umap_coords$cluster <- as.factor(clusters)
  umap_coords$name <- rownames(umap_coords)
  umap_coords <- umap_coords %>% arrange(desc(type))

  
  return(umap_coords)
}

#' Plot biUMAP
#' 
#' @param umap_coords data frame as outputted by `run_biUMAP_*`
#' @param color_by Either "type" or "cluster". "type" colors by the type 
#' (cell or gene) while "cluster" colors by the assigned cluster.
#' 
#' @return 
#' ggplot of UMAP
#' 
#' @export
plot_biUMAP <- function(umap_coords, color_by = "type"){
  
  if(color_by == "type"){
    p <- ggplot(umap_coords, aes(x=x, y=y, color = type,
                                 text = paste0(
                                   "Type: ", type, "\n",
                                   "Name: ", name, "\n",
                                   "Cluster: ", cluster))) +
      geom_point(alpha = 0.4) +
      theme_bw()
  } else if (color_by == "cluster"){
    
    p <- ggplot(umap_coords, aes(x=x, y=y, color = cluster,
                                 text = paste0(
                                   "Type: ", type, "\n",
                                   "Name: ", name, "\n",
                                   "Cluster: ", cluster)))+
      geom_point(alpha = 0.4) +
      theme_bw()
  } else {
    stop("color_by has to be either 'type' or 'cluster'.")
  }

  
  return(p)
  
}


