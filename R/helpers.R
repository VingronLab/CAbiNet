
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
#' @param direction Either "cells" or "genes". If "cells-genes" the association ratio
#' of each gene for each cell is calculated. If "genes-cells" the assoc. ratio of each
#' cell for each gene is calculated.
#'
#' @return
#' Matrix with pairwise association ratios.
calc_assR <- function(caobj, direction = "cells-genes"){
  
  stopifnot(is(caobj, "cacomp"))
  
  if (direction == "cells-genes") {
    assR <- caobj@std_coords_cols %*% t(caobj@prin_coords_rows)
  } else if (direction == "genes-cells"){
    assR <- caobj@std_coords_rows %*% t(caobj@prin_coords_cols)
  } else if (direction == "cells-cells"){
    assR <- caobj@std_coords_cols %*% t(caobj@std_coords_cols)
  } else if (direction == "genes-genes"){
    assR <- caobj@std_coords_rows %*% t(caobj@std_coords_rows)
  } else {
    stop("Please choose a valid direction.")
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
  cell_gene_assr <- calc_assR(caobj, direction = "cells-genes")
  #   gene_cell_assr <- calc_assR(caobj, direction = "genes")
  # no need to calculate again, just take the transpose
  gene_cell_assr <- t(cell_gene_assr) 
  
  return(list("cc" = cell_dists,
              "gg" = gene_dists,
              "cg" = cell_gene_assr,
              "gc" = gene_cell_assr))
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





#' Detecting eigengap.
#' @description 
#' Searching for the largest gap between conjuncted eigenvalues which are sorted 
#' in an ascending order, the most important eigenvalues and eigenvectors are the ones 
#' before the detected eigengap.
#' @param e Eigenvalues of type 'vector'
#' @param v Eigenvector matirx of type 'matrix/array' whose columns are single eigenvectors.
#' @return 
#' The selected most important eigenvectors for downstream clustering, of type 'matrix/array'.
#' 
eigengap = function(e, v){
  # e is eigenvalues, v is eigenvector
  
  idx = order(e, decreasing = FALSE)
  
  # order eigenvalues and eigenvectors
  v = v[, idx]
  e = e[idx]
  
  n = length(e)
  gaps = abs(e[1:(n-1)] - e[2:n])
  idx = which(gaps == max(gaps))[1]
  
  message(paste0('Using eigengap at index ', idx, '\n'))
  
  if ( idx <2){
    selected_eigen = v[, 1:2]
  }else{
    selected_eigen = v[, 1:idx]
  }
  
  return(selected_eigen)
}





#' calculates association ratio between to columns in a matrix.
#' 
#' @param matrix a matrix
#' @param origin origin column (name or index)
#' @param target target column (name or index)
#' 
#' @returns
#' Returns the association ratio between the columns from origin to target.
aR_metric <- function(matrix, origin, target){
  
  assR <- matrix[,origin] %*% matrix[,target]
  assR <- drop(assR)
  return(assR)
  # assR <- caobj@std_coords_rows %*% t(caobj@prin_coords_cols)
}



#' Helper function to check if object is empty.
#' @param x object
#' @return TRUE if x has length 0 and is not NULL. FALSE otherwise
is.empty <- function(x) return(isTRUE(length(x) == 0 & !is.null(x)))

#' Determines majority in a vector
#' @description 
#' Changed version of mclust::majorityVote. Ties are broken randomly.
#' 
#' @param x a vector
#' 
#' @returns 
#' A single element of x that has the highest count.
#' 
get_majority <- function(x){
  
  x <- as.vector(x)
  tally <- table(x)
  max_idx <- seq_along(tally)[tally == max(tally, na.rm = TRUE)]
  
  if(length(max_idx) > 1){
    max_idx <- sample(max_idx, size = 1)
  }
  
  majority <- names(tally)[max_idx]
  
  return(majority)
  
}
