

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
#' @param caobj
calc_assR <- function(caobj){}



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
  #   nn.matrix <- as(object = sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x = object), ncol(x = object))), Class = "Graph")

  rownames(nn.matrix) <- rownames(dists)
  colnames(nn.matrix) <- colnames(dists)
  return(nn.matrix)
  # neighbor.graphs <- list(nn = nn.matrix, nn.ranked = nn.ranked)
}


#' Combine kNN graphs to large cell-gene adjecency matrix
create_bigraph <- function(sample_dists,
                           gene_dists,
                           sample_gene_dists,
                           gene_sample_dists,
                           sk,
                           gk,
                           sgk,
                           gsk,
                           overlap = 0.2,
                           prune_overlap = TRUE) {

  sgg_nn <- knn_sg(sample_gene_dists,
                   k.param = sgk,
                   decr = TRUE)

  sgg_nn <- as.matrix(sgg_nn)
  idx <- colSums(sgg_nn) > 0

  sgg_nn <- sgg_nn[,idx]
  gene_dists <- gene_dists[idx,idx]
  gene_sample_dists <- gene_sample_dists[idx,]

  ssg_nn <- knn_sg(sample_dists,
                   k.param = sk,
                   decr = FALSE)

  overlap_mat <- prune_genes(cg_adj = sgg_nn,
                             cc_adj = ssg_nn)

  sgg_nn[overlap_mat < overlap] <- 0

  idx <- colSums(sgg_nn) > 0
  sgg_nn <- sgg_nn[,idx]
  gene_dists <- gene_dists[idx,idx]
  gene_sample_dists <- gene_sample_dists[idx,]

  ggg_nn <- knn_sg(gene_dists,
                   k.param = gk,
                   decr = FALSE)$nn


  gsg_nn <- t(sgg_nn)
  #   gsg_nn <- knn_sg(gene_sample_dists,
  #                    k.param = gsk,
  #                    decr = TRUE)$nn

  GSG_1 <- cbind(ssg_nn, sgg_nn)
  GSG_2 <- cbind(gsg_nn, ggg_nn)

  GSG <- rbind(GSG_1, GSG_2)
  return(GSG)
}

#' Run SNN on graph
create_SNN <- function() {
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
