
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




#' 
#' Determines the overlap of chosen genes between nearest neighbour cells.
#' It works faster on sparse matrix than on dense matrix.
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

  overlap_mat <- cc_adj %*% cg_adj
  overlap_mat <- overlap_mat/rowSums(cc_adj)

  overlap_mat <- overlap_mat * cg_adj
  
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


#' Shuffle rows of a data frame for better plotting.
#' @param df data.frame
#' @export
mix <- function(df){
  df <- df[sample(seq_len(nrow(df)), size = nrow(df)),]
}



rowNorm <- function(x){
    
    if (is.matrix(x)){
        norm <- sqrt(rowSums(x^2))
    } else if (is.null(dim(x))) {
        norm <- sqrt(sum(x^2))
    } else {
        stop("Unknown object.")
    }
    
    return(norm)
}


#' Transforms vectors such that the max. inner product (MIP) search is
#' equal to a NN search with euclidean distances.
#' 
#' @param vectors Matrix of row vectors
#' 
#' @returns 
#' A new matrix of row vectors with one additional (first) dimension that is
#' the square root of the difference between the squared max. vector norm and 
#' the squared vector norm. This makes the norm of all vectors effectively
#' the same.
#' 
#' @details 
#' See also:
#' https://github.com/benfred/implicit/blob/42832574f1a29c71b3263e219fc471fc97328552/implicit/utils.py#L60
#' https://towardsdatascience.com/maximum-inner-product-search-using-nearest-neighbor-search-algorithms-c125d24777ef
#' 
#' @references  
#' Yoram Bachrach, Yehuda Finkelstein, Ran Gilad-Bachrach, Liran Katzir,
#'  Noam Koenigstein, Nir Nice, Ulrich Paquet. 
#'  Speeding up the Xbox recommender system using a euclidean transformation 
#'  for inner-product spaces. RecSys 2014
augment_vector <- function(vectors){
    
    rnorms <- rowNorm(vectors)
    max_norm <- max(rnorms)
    
    extra_dim <- sqrt(max_norm^2-rnorms^2)
    
    return(cbind(extra_dim, vectors))
    
}

#' Add an extra 0 dimension to vector
#' @param vectors Matrix of row vectors
#' 
add_zero_dim <- function(vectors){
    return(cbind(0, vectors))
}


#' Convert index matrix to sparse matrix.
#' 
#' @param indx_mat row-wise index matrix. Usually a kNN graph.
#' @param row_names the rownames of the resulting sparse matrix.
#' @param col_names the column names that correspond to the indexes in indx_mat.
#' 
#' @returns 
#' a `Matrix::sparseMatrix` object.
#' 
indx_to_spmat <- function(indx_mat,
                          row_names,
                          col_names){
    
    j <- as.numeric(t(indx_mat))
    i <- ((1:length(j)) - 1) %/% ncol(indx_mat) + 1
    
    nn.matrix <- Matrix::sparseMatrix(i = i,
                                      j = j,
                                      x = 1,
                                      dims = c(length(row_names), length(col_names)))
    
    rownames(nn.matrix) <- row_names
    colnames(nn.matrix) <- col_names
    
    return(nn.matrix)
}



# 
# add2knn <- function(vec, to_add){
#     
#     nvec <- length(vec)
#     
#     vec <- na.omit(vec)
#     vec <- union(vec, to_add)
#     
#     if (length(vec) > nvec) {
#         stop("adding too many elements!")
#     } else {
#         
#         nna <- nvec-length(vec)
#         vec <- c(vec, rep(NA_integer_, nna))
#     }
#     
#     return(vec)
#     
# }