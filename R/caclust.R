

#' Calculate euclidean distance
#'
#' @description Calculates euclidean distance.
#'
#' @param mat symmetric matrix.
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
    # sgg_nn <- as.matrix(sgg_nn)
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

#' Run SNN on graph
create_SNN <- function(graph) {

}

#' Leiden clustering on bigraph
run_leiden <- function() {
}

#' run spectral clustering
run_spectral <- function() {
}

#' Run biclustering
caclust <- function() {
}
