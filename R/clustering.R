

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
calc_distances <- function(caobj, appx = TRUE){
  
  cell_dists <- calc_euclidean(caobj@prin_coords_cols)
  gene_dists <- calc_euclidean(caobj@prin_coords_rows)
  if (isTRUE(appx)){
    cell_gene_assr <- calc_assR(caobj, direction = "cells-genes")
  }else{
    cell_gene_assr <- t(caobj@AR)
  }
#   gene_cell_assr <- calc_assR(caobj, direction = "genes")
    # no need to calculate again, just take the transpose
  gene_cell_assr <- t(cell_gene_assr) 

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
#' @param loops TRUE/FALSE. If TRUE self-loops are allowed, otherwise not.
#'
#' @return
#' adjacency matrix of kNN graph: sparse matrix of type "dgCMatrix".
#'
#' @export
make_knn <- function(dists,
                     k,
                     decr = TRUE,
                     loops = FALSE) {

  n.row <- nrow(dists)
  n.col <- ncol(dists)

  if ((n.col-1 < k) |(n.row-1 < k)) {
    warning(
      "k set larger than number of selected genes Setting k to number of genes - 1.",
      call. = FALSE
    )
    k = min(n.col-1, n.row -1)
    # k <- n.col - 1
  }
  knn.mat <- matrix(data = 0, ncol = k, nrow = n.row)
  # knd.mat <- knn.mat
  
  if (isTRUE(loops)){
    nns <- seq_len(k)
    
  } else if (isFALSE(loops)) {
    nns <- seq_len(k)+1
  } else {
    stop()
  }
  
  for (i in 1:n.row) {
    knn.mat[i, ] <- order(dists[i, ], decreasing = decr)[nns]
    # knd.mat[i, ] <- dists[i, knn.mat[i, ]]
  }
  nn.ranked <- knn.mat[,seq_len(k)]

  # convert nn.ranked into a Graph
  j <- as.numeric(t(knn.mat))
  i <- ((1:length(j)) - 1) %/% k + 1
  nn.matrix <- Matrix::sparseMatrix(i = i,
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
#' @param loops TRUE/FALSE. If TRUE self-loops are allowed, otherwise not.
#' @param select_genes TRUE/FALSE. Should genes be selected by whether they have
#' an edge in the cell-gene kNN graph?
#' @param prune_overlap TRUE/FALSE. If TRUE edges to genes that share less
#' than `overlap` of genes with the nearest neighbours of the cell are removed.
#' Pruning is only performed if select_genes = TRUE.
#' @param overlap Numeric between 0 and 1. Overlap cutoff if
#' prune_overlap = TRUE.
#' @param calc_gene_cell_kNN TRUE/FALSE. If TRUE a cell-gene graph is calculated
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
                           loops = FALSE,
                           select_genes = TRUE,
                           prune_overlap = TRUE,
                           overlap = 0.2,
                           calc_gene_cell_kNN = FALSE,
                           marker_genes = NULL) {

  cgg_nn <- make_knn(cell_gene_assr,
                     k = k_cg,
                     decr = TRUE,
                     loops = loops)
  
  if (!is.null(marker_genes)){
    stopifnot(is(marker_genes, "character"))
    
    idx <- which(colnames(cgg_nn) %in% marker_genes)
    marker_knn <- cgg_nn[,idx]
    cgg_nn <- cgg_nn[,-idx]
    
    marker_dists <- gene_dists[idx,]
    marker_assr <- gene_cell_assr[idx,]
    
    gene_dists <- gene_dists[-idx, -idx]
    gene_cell_assr <- gene_cell_assr[-idx,]
  }

  if(isTRUE(select_genes)){
    
    idx <- Matrix::colSums(cgg_nn) > 0
    cgg_nn <- cgg_nn[,idx]
    gene_dists <- gene_dists[idx,idx]
    gene_cell_assr <- gene_cell_assr[idx,]
  }


  ccg_nn <- make_knn(cell_dists,
                     k = k_c,
                     decr = FALSE,
                     loops = loops)


  if(isTRUE(select_genes) & isTRUE(prune_overlap)){

    overlap_mat <- determine_overlap(cg_adj = cgg_nn,
                                     cc_adj = ccg_nn)
    
    cgg_nn[overlap_mat < overlap] <- 0

    idx <- Matrix::colSums(cgg_nn) > 0
    if (sum(idx) <=1){
      cat('The given value of "overlap" is too large, all gene nodes 
           are trimmed. Use a smaller overlap instead! 99% quantile of overlaps is',
          quantile(overlap_mat,0.99), '. Overlap should be smaller than ', max(overlap_mat), '.\n')
      return(ccg_nn)
    }
    cgg_nn <- cgg_nn[,idx]
    gene_dists <- gene_dists[idx,idx]
    gene_cell_assr <- gene_cell_assr[idx,]

  }

  
  if(!is.null(marker_genes)){
    
    cgg_nn <- cbind(cgg_nn, marker_knn)
    
    marker_dists <- marker_dists[,c(colnames(gene_dists), rownames(marker_dists))]
    gene_dists <- cbind(rbind(gene_dists, marker_dists[,colnames(gene_dists)]), t(marker_dists))
    gene_cell_assr <- rbind(gene_cell_assr, marker_assr)
    
  }

  ggg_nn <- make_knn(gene_dists,
                     k = k_g,
                     decr = FALSE,
                     loops = loops)

  if(isFALSE(calc_gene_cell_kNN)){
    gcg_nn <- Matrix::t(cgg_nn)

  } else if(isTRUE(calc_gene_cell_kNN)){
    gcg_nn <- make_knn(gene_cell_assr,
                       k = k_gc,
                       decr = TRUE,
                       loops = loops)
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
#' @param SNN_prune numeric. Value between 0-1. Sets cutoff of acceptable jaccard 
#' similarity scores for neighborhood overlap of vertices in SNN. Edges with values 
#' less than this will be set as 0. The default value is 1/15.
#' @param select_genes TRUE/FALSE. Should genes be selected by wether they have
#' an edge in the cell-gene kNN graph?
#' @param prune_overlap TRUE/FALSE. If TRUE edges to genes that share less
#' than `overlap` of genes with the nearest neighbours of the cell are removed.
#' @param overlap Numeric between 0 and 1. Overlap cutoff if
#' prune_overlap = TRUE.
#' @param calc_cell_gene_kNN TRUE/FALSE. If TRUE a cell-gene graph is calculated
#' by choosing the `k_gc` nearest cells for each gene. If FALSE the cell-gene
#' graph is transposed.
#' @param mode The type of neighboring vertices to use for calculating similarity
#'  scores(Jaccard Index). Three options: "out", "in" and "all":
#' * "out": Selecting neighbouring vertices by out-going edges;
#' * "in": Selecting neighbouring vertices by in-coming edges;
#' * "all": Selecting neigbouring vertices by both in-coming and out-going edges.
#' @inheritParams create_bigraph
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
                       SNN_prune = 1/15,
                       loops = FALSE,
                       mode = "out",
                       select_genes = TRUE,
                       prune_overlap = TRUE,
                       overlap = 0.2,
                       calc_gene_cell_kNN = FALSE,
                       marker_genes = NULL) {
  
  
  stopifnot(all(c("cc", "gg", "cg", "gc") %in% names(distances)))
  stopifnot(mode %in% c("out", "in", "all"))
  # S = caobj@U %*% diag(caobj@D) %*% t(caobj@V)
  adj <- create_bigraph(cell_dists = distances[["cc"]],
                        gene_dists = distances[["gg"]],
                        cell_gene_assr = distances[["cg"]],
                        # cell_gene_assr = t(S),
                        gene_cell_assr = distances[["gc"]],
                        # gene_cell_assr = S,
                        k_c = k_c,
                        k_g = k_g,
                        k_cg = k_cg,
                        k_gc = k_gc,
                        loops = loops,
                        overlap = overlap,
                        prune_overlap = prune_overlap,
                        select_genes = select_genes,
                        calc_gene_cell_kNN = calc_gene_cell_kNN,
                        marker_genes = marker_genes)
  
  if(!is(adj, "dgCMatrix")){
    adj <- as(adj, "dgCMatrix")  
  }
  
  snn.matrix <- ComputeSNNasym(adj, SNN_prune, mode = mode)
  
  ## to coincide with output of "igraph"
  diag(snn.matrix) = 1
  
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
#' @param rand_seed Random seed.
#' 
#' @return 
#' vector of type `factor.` Assigned clusters of cells and genes. 
#' The names of cells and genes are saved in the names of the 
#' vector. 
#' 
run_leiden <- function(SNN, 
                       resolution = 1,
                       n.int = 10, 
                       rand_seed = 2358,
                       dense = TRUE) {
  if(isTRUE(dense)){
    SNN = as.matrix(SNN)
  }
  
  clusters <- leiden::leiden(object = SNN,
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

#' Calculate Normalized Graph Laplacian
#' @description 
#' Calculate Normalized graph laplacian of the input adjacency matrix, the graph laplacian
#' L = I-D^(-1/2)AD^(-1/2), where D is a diagonal matrix with row sums of A as entries.
#' @param adj The adjacency matix of type 'matrix/array'
#' @return 
#' A sparsed normalized graph laplacian matrix of type "dgCmatrix".
#' 
NormLaplacian = function(adj){
  
  D = Matrix::rowSums(adj)
  mat_D = Matrix::sparseMatrix(i = seq_len(length(D)),
                       j = seq_len(length(D)),
                       x = 1/sqrt(D),
                       dims =c(nrow(adj), ncol(adj)) )
  mat_D[is.na(mat_D)] = 0
  mat_D[is.infinite(mat_D)] = 0
  
  I = Matrix::sparseMatrix(i = seq_len(length(D)),
                           j = seq_len(length(D)),
                           x = 1,
                           dims =c(nrow(adj), ncol(adj)) )
  
  L = I - mat_D %*% adj %*% mat_D
  
  return(L)
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
  # cat('eigengap idx',idx,'\n')
  if ( idx <2){
    selected_eigen = v[, 1:2]
  }else{
    selected_eigen = v[, 1:idx]
  }
  
  return(selected_eigen)
}

#' An integration of several runs of skmens with different random seeds
#' @description 
#' This function will select the optimal clustering result from several skmeans::skmeans runs
#' with different random seeds, The clustering result with smallest within-cluster-sum
#' of squared distances will be selected.
#' @param method see skmeans::skmeans function
#' @param m A number not less than 1 controlling the softness of the partition.(see skmeans::skmeans 'm')
#' @param weights A numeric vector of non-negative case weights.(see skmeans::skmeans 'm')
#' @param control A list of control parameters for different methods.
#' @param num.seeds number of trials with random seeds
#' @return 
#' The optimal skmeans clustering result
#' @export
SKMeans <- function (x, k, 
                     method = NULL, 
                     m = 1, 
                     weights = 1, 
                     control = list(), 
                     num.seeds = 10) 
{
  if (mode(x) == "numeric") 
    # x <- data.frame(new.x = x)
  KM <- skmeans::skmeans(x = x, 
                         k = k, 
                         method = method,
                         m = m,
                         weights = weights,
                         control = control)
  wss = do.call(sum, lapply(list(1:length(KM$cluster)), function(i){
    sum((x[i,]-KM$prototypes[KM$cluster[i],])^2)
  }))
  for (i in 2:num.seeds) {
    newKM <- skmeans::skmeans(x = x, 
                              k = k, 
                              method = method,
                              m = m,
                              weights = weights,
                              control = control)
    
    newwss = do.call(sum, lapply(list(1:length(newKM$cluster)), function(i){
      sum((x[i,]-newKM$prototypes[newKM$cluster[i],])^2)
    }))
    if (newwss <wss) {
      KM <- newKM
      wss = newwss
    }
  }
  # xmean <- apply(x, 2, mean)
  # centers <- rbind(KM$centers, xmean)
  # bss1 <- as.matrix(dist(centers)^2)
  # KM$betweenss <- sum(as.vector(bss1[nrow(bss1), ]) * c(KM$size, 
  #                                                       0))
  return(KM)
}

#' run spectral clustering
#' @description 
#' This function is designed for detecting clusters from input graph adjacency 
#' matrix by using spectral clustering with normalized graph laplacian.
#' @param SNN The adjacency matrix of graph
#' @param use_gap TRUE/FALSE. If TRUE, 'eigengap' method will be used to find the
#' most important eigenvector automatically, and the number of output clusters 
#' equals number of selected eigenvectors. If FALSE, 'nclust'(integer) should be given. 
#' The eigenvectors corresponding with the smallest 'nclust' eigenvalues will be
#' selcted and 'nclust' clusters will be detected by skmeans.
#' @param python TRUE/FALSE. If TRUE, pytorch function will be used to do eigenvalue
#' decompositon, else R-base function 'svd' will be used for calculation. 
#' @return 
#' The clustering results
#' 
run_spectral <- function(SNN, 
                         use_gap = TRUE, 
                         nclust = NULL,
                         clust.method = 'kmeans',
                         iter.max=10, 
                         num.seeds=10,
                         return.eig = TRUE,
                         dims = 30) {
  diag(SNN) = 0
  L = NormLaplacian(SNN)
  
  if (dims > ncol(SNN)){
    
    dims = ncol(SNN)
    
  }
  
    
  SVD <- irlba::irlba(L, nv =dims, smallest = TRUE) # eigenvalues in a decreasing order
  names(SVD)[1:3] <- c("D", "U", "V")
  
  idx = order(SVD$D, decreasing = FALSE)
  eigenvalues = SVD$D[idx]
  eigenvectors = SVD$U[,idx]
  
  # cat('SVD for graph laplacian is done....\n')
  
  if (use_gap == FALSE){
    
    if (is.null(nclust)){
      
        stop('Number of selected eigenvectors of lapacian is required, change value of nclust as an integer!')
    
      # }else if (nclust > dims){
      #   
      #   stop('Number of dims should be larger than number of clusters (nclust)')
      
      }else{
      
        eig = eigenvectors[,1:nclust] # in an increasing order
        # eig = eigenvectors
        
    }
    
    } else if (use_gap == TRUE){
      
      
    eig = eigengap(eigenvalues, eigenvectors)# in an increasing order
    nclust = ncol(eig)
    
  }
  
  if (clust.method == 'skmeans'){
    
    clusters = SKMeans(eig, k = nclust, num.seeds = num.seeds)$cluster
  
  }else if (clust.method == 'kmeans'){
    
  clusters = RcmdrMisc::KMeans(eig, centers = ncol(eig), iter.max=iter.max, num.seeds= num.seeds)$cluster

  }else if (clust.method == 'GMM'){
    
    gmm = ClusterR::GMM(eig,
                        gaussian_comps = ncol(eig),
                        dist_mode = "maha_dist",
                        seed_mode = "random_subset",
                        km_iter = 30,
                        em_iter = 30,
                        verbose = F)   
    
    cluster_probs = predict_GMM(data = eig,
                                 CENTROIDS = gmm$centroids, 
                                 COVARIANCE = gmm$covariance_matrices,
                                 WEIGHTS = gmm$weights)
    
    cluster_probs <- cluster_probs[-which(names(cluster_probs)=="log_likelihood")]
    rownames(cluster_probs$cluster_proba) <- rownames(SNN)
    colnames(cluster_probs$cluster_proba) <- paste("BC", seq_len(ncol(cluster_probs$cluster_proba)))
    cluster_probs$cluster_labels <- as.factor(cluster_probs$cluster_labels)
    names(cluster_probs$cluster_labels) <- rownames(SNN)
    
    return(cluster_probs)
    
  }else{
  stop('clustering method should be chosen from kmeans and skmeans!')
  }
  
  clusters <- as.factor(clusters)
  names(clusters) <- rownames(SNN)
  
  if (return.eig){
    
    if (is.null(dims)){
      dims = min(30, ncol(SNN))
    }
    
    return(list(clusters = clusters,
                eigen = eigenvectors[,1:dims]))
  }else{
    return(clusters)
  }
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
#' @inheritParams create_bigraph
#' @inheritParams create_SNN
#' @inheritParams run_leiden
#' @inheritParams run_spectral
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
                        SNN_prune = 1/15,
                        loops = FALSE,
                        mode = "out",
                        select_genes = TRUE,
                        prune_overlap = TRUE,
                        overlap = 0.2,
                        calc_gene_cell_kNN = FALSE,
                        resolution = 1,
                        marker_genes = NULL,
                        n.int = 10,
                        rand_seed = 2358,
                        use_gap = TRUE,
                        nclust = NULL,
                        python = TRUE,
                        clust.method = 'kmeans',
                        iter.max=10, 
                        num.seeds=10,
                        return.eig = TRUE,
                        sc.dims = NULL,
                        leiden.dense = TRUE,
                        appx = TRUE) {
  
  call_params <- as.list(match.call())
  names(call_params)[1] <- "Call"
  
  distances <- calc_distances(caobj = caobj, appx = appx)
  
  SNN <- create_SNN(caobj = caobj, 
                    distances = distances,
                    k_c = k_sym,
                    k_g = k_sym,
                    k_cg = k_asym,
                    k_gc = k_asym,
                    loops = loops,
                    mode = mode,
                    SNN_prune = SNN_prune,
                    select_genes = select_genes,
                    prune_overlap = prune_overlap,
                    overlap = overlap,
                    calc_gene_cell_kNN = calc_gene_cell_kNN,
                    marker_genes = marker_genes)
  
  if (algorithm == "leiden"){
    
    clusters <- run_leiden(SNN = SNN, 
                           resolution = resolution,
                           n.int = n.int, 
                           rand_seed = rand_seed,
                           dense = leiden.dense)
    
    eigen <- matrix()
    
  } else if (algorithm == "spectral"){
    
    if (is.null(sc.dims)){
      sc.dims = length(caobj@D)
    }
      clusters <- run_spectral(SNN = SNN,
                               use_gap = use_gap,
                               nclust = nclust,
                               python = python,
                               clust.method = clust.method,
                               iter.max = iter.max, 
                               num.seeds = num.seeds,
                               return.eig = return.eig,
                               dims = sc.dims)


    if (isTRUE(return.eig)){
      
      eigen <- clusters$eigen
      rownames(eigen) <- rownames(SNN)
      
    } else {
      
      eigen <- matrix()
    }
    
    if (clust.method == "GMM"){
      
      cell_idx <- which(names(clusters$cluster_labels) %in% rownames(caobj@prin_coords_cols))
      gene_idx <- which(names(clusters$cluster_labels) %in% rownames(caobj@prin_coords_rows))
      
      cell_clusters <- clusters$cluster_labels[cell_idx]
      gene_clusters <- clusters$cluster_labels[gene_idx]
      
      cell_idx <- which(rownames(clusters$cluster_proba) %in% rownames(caobj@prin_coords_cols))
      gene_idx <- which(rownames(clusters$cluster_proba) %in% rownames(caobj@prin_coords_rows))
      cell_prob <- clusters$cluster_proba[cell_idx,]
      gene_prob <- clusters$cluster_proba[gene_idx,]
      
      caclust_res <- do.call(new_caclust, list("cell_clusters" = cell_clusters,
                                               "gene_clusters" = gene_clusters,
                                               "SNN" = SNN,
                                               "eigen" = eigen,
                                               "parameters" = call_params,
                                               "cell_prob" = cell_prob,
                                               "gene_prob" = gene_prob))
      return(caclust_res)
    } else {
      clusters <- as.factor(clusters$clusters)
    }

      

    
  } else{
    stop("algorithm should choose from 'leiden' and 'spectral'!")
  }

  cell_idx <- which(names(clusters) %in% rownames(caobj@prin_coords_cols))
  gene_idx <- which(names(clusters) %in% rownames(caobj@prin_coords_rows))
  
  cell_clusters <- clusters[cell_idx]
  gene_clusters <- clusters[gene_idx]
  
  
  caclust_res <- do.call(new_caclust, list("cell_clusters" = cell_clusters,
                                           "gene_clusters" = gene_clusters,
                                           "SNN" = SNN,
                                           "eigen" = eigen,
                                           "parameters" = call_params))
  return(caclust_res)
}


#' Plots UMAP depicting both cells and genes.
#' 
#' @description 
#' TODO
#' 
#' @param caobj A cacomp object with principal and standard coordinates 
#' calculated.
#' @param caclust_obj results from biclustering of class "caclust"
#' @param k_umap Number of nearest neighbours to use from the SNN graph for
#' UMAP
#' @param rand_seed Random seed for UMAP.
#' 
#' @return 
#' data frame containing the UMAP coordinates of cells and genes.
#' 
#' @export 
run_biUMAP <- function(caobj,
                              caclust_obj,
                              k_umap,
                              rand_seed = 2358,
                              algorithm = 'SNNdist',
                              coords = 1,...){

  
  stopifnot(is(caobj, "cacomp"))
  stopifnot(is(caclust_obj, "caclust"))
  stopifnot(algorithm %in% c("SNNgraph", "SNNdist", "spectral", "ca", "ca_assR"))
  
  custom.config = umap::umap.defaults
  custom.config$random_state = rand_seed
  
  if (algorithm == 'SNNgraph'){
  
    SNN <- get_snn(caclust_obj)
    
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

    snn_umap_graph = umap::umap.knn(indexes = SNN_idx,
                                    distances = SNN_jacc)

    assym <- rbind(caobj@std_coords_cols, caobj@prin_coords_rows)
    assym <- assym[rownames(assym) %in% rownames(SNN),]

    caclust_umap = umap::umap(assym,
                              config = custom.config,
                              n_neighbors = k_umap,
                              knn = snn_umap_graph)
    
    umap_coords <- as.data.frame(caclust_umap$layout)
    
    
  }else if (algorithm == "SNNdist"){
    
    SNNdist <- as.matrix(1-get_snn(caclust_obj))
    
    reticulate::source_python(system.file("inst/python/umap.py", package = "CAclust"), envir = globalenv())
    
    umap_coords = python_umap(dm = SNNdist,
                              metric = "precomputed",
                              n_neighbors = as.integer(k_umap))
    
    umap_coords <- as.data.frame(umap_coords)
    rownames(umap_coords) <- colnames(SNNdist)
  
  }else if (algorithm == 'spectral'){
    
    eigen = get_eigen(caclust_obj)
    
    if(is.na(eigen)) stop("Spectral clustering not run.")

    
    caclust_umap = umap::umap(eigen, 
                              config = custom.config,
                              n_neighbors = k_umap)
    
    umap_coords <- as.data.frame(caclust_umap$layout)
    
  }else if (algorithm == 'ca'){
    # eigen = rbind(caobj@prin_coords_cols, caobj@prin_coords_rows)
    if (coords == 1){
      
      if(sum(!is.null(caobj@prin_coords_rows), !is.null(caobj@std_coords_cols)) != 2){
        stop("Principal and/or standard coordinates not found, ",
             "please run ca_coords() first!")
      }
      eigen = rbind(caobj@prin_coords_cols, caobj@prin_coords_rows)
      
    } else if (coords == 2){
      if(sum(!is.null(caobj@prin_coords_cols), !is.null(caobj@std_coords_rows)) != 2){
        stop("Principal and/or standard coordinates not found, ",
             "please run ca_coords() first!")
      }
      eigen = rbind(caobj@std_coords_cols, caobj@std_coords_rows)
    }else if (coords == 3){
      if(sum(!is.null(caobj@U), !is.null(caobj@V)) != 2){
        stop("Singular eigenvectors not found, ",
             "please run ca_coords() first!")
      }
      eigen = rbind(caobj@V, caobj@U)
    } else {
      stop("princ_coords must be either 1 for rows or 2 for columns.")
    }
    

    idx1 = which(rownames(eigen)  %in% names(gene_clusters(caclust_obj)))
    idx2 = which(rownames(eigen) %in% names(cell_clusters(caclust_obj)))
    
    eigen = eigen[c(idx1,idx2),]
    SNN <- get_snn(caclust_obj)
    eigen <- eigen[rownames(eigen) %in% rownames(SNN),]
    
    
    caclust_umap = umap::umap(eigen, 
                              config = custom.config,
                              metric = 'cosine',
                              n_neighbors = k_umap)
    umap_coords <- as.data.frame(caclust_umap$layout)
    
  }else if (algorithm == 'ca_assR'){
    
    SNN <- get_snn(caclust_obj)
    
    cc <- calc_assR(caobj = caobj, direction = "cells-cells")
    cc <- cc[rownames(cc) %in% rownames(SNN), colnames(cc) %in% colnames(SNN)]
    
    cg <- calc_assR(caobj = caobj, direction = "cells-genes")
    cg <- cg[rownames(cg) %in% rownames(SNN), colnames(cg) %in% colnames(SNN)]
    
    gg <- calc_assR(caobj = caobj, direction = "genes-genes")
    gg <- gg[rownames(gg) %in% rownames(SNN), colnames(gg) %in% colnames(SNN)]
    
    gc <- calc_assR(caobj = caobj, direction = "genes-cells")
    gc <- gc[rownames(gc) %in% rownames(SNN), colnames(gc) %in% colnames(SNN)]
    
    assR <- rbind(cbind(cc,cg), cbind(gc,gg))

    
    
    reticulate::source_python(system.file("python/umap.py", package = "CAclust"), envir = globalenv())
    
    umap_coords = python_umap(dm = assR,
                              metric = "precomputed",
                              n_neighbors = as.integer(k_umap))
    
    umap_coords <- as.data.frame(umap_coords)
    rownames(umap_coords) <- colnames(assR)
    
  } else {
    stop()
  }
  
  
  cellc <- cell_clusters(caclust_obj)
  genec <- gene_clusters(caclust_obj)
  
  
  colnames(umap_coords) <- c("x", "y")
  umap_coords$name <- rownames(umap_coords)
  
  umap_coords$type <- NA
  umap_coords$type[umap_coords$name %in% names(cellc)] <- "cell" 
  umap_coords$type[umap_coords$name %in% names(genec)] <- "gene" 

  umap_coords$cluster <- NA
  cell_idx <- na.omit(match(names(cellc), umap_coords$name))
  gene_idx <- na.omit(match(names(genec), umap_coords$name))
  
  umap_coords$cluster[cell_idx] <- cell_clusters(caclust_obj)
  umap_coords$cluster[gene_idx] <- gene_clusters(caclust_obj)
  umap_coords$cluster <- as.factor(umap_coords$cluster)

  umap_coords <- umap_coords %>% dplyr::arrange(desc(type))

  
  return(umap_coords)
}





aR_metric <- function(matrix, origin, target){
  
  assR <- matrix[,origin] %*% matrix[,target]
  assR <- drop(assR)
  return(assR)
  # assR <- caobj@std_coords_rows %*% t(caobj@prin_coords_cols)
}

assign_clusters_GMM <- function(caclust_obj, type = "genes", cutoff=0.5){
  if(type == "genes"){
    prob_slot <- "gene_prob"
  } else if (type == "cells"){
    prob_slot <- "cells_prob"
  }
  
  probs <- slot(caclust_obj, prob_slot)
  probs <- probs > cutoff
  
  return(probs)
}

