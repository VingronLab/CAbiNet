#' @include classes.R
NULL



#' Calculate Normalized Graph Laplacian
#' @description
#' Calculate Normalized graph laplacian of the input adjacency matrix, the graph laplacian
#' L = I-D^(-1/2)AD^(-1/2), where D is a diagonal matrix with row sums of A as entries.
#' @param adj The adjacency matix of type 'matrix/array'.
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



#' An integration of several runs of skmens with different random seeds
#' @description
#' This function will select the optimal clustering result from several
#' skmeans::skmeans runs with different random seeds, The clustering result with the smallest
#' within-cluster-sum of squared distances will be selected.
#' @param k Integer. Number of cluters to detect for skmeans.
#' @param x Matrix. This function will cluster the rows of the input matrix.
#' @param num_seeds Integer. Number of trials with random seeds
#' @inheritParams skmeans::skmeans
#'
#' @return
#' Returns an object inheriting from classes skmeans and pclust (see ?skmeans) which
#' gives the local optimal skmeans clustering result within several trials.
#'
optimal_skm <- function (x,
                        k,
                        num_seeds = 10,
                        method = NULL,
                        m = 1,
                        weights = 1,
                        control = list()){

  if(!is(x, 'matrix')){
    x = as.matrix(x)
  }

  if (num_seeds < 1){
    stop('num_seeds should be larger than 0.')
  }

  wcs <- Inf

  for (i in 1:num_seeds) {
    newres <- skmeans::skmeans(x = x,
                               k = k,
                               method = method,
                               m = m,
                               weights = weights,
                               control = control)

    newwcs <- do.call(sum,
                      lapply(list(1:length(newres$cluster)),
                                  function(i){
                                    sum((x[i,]-newres$prototypes[newres$cluster[i],])^2)
                                    }))

    if (newwcs <= wcs) {
      res <- newres
      wcs <- newwcs
    }
  }

  return(res)
}

#' An integration of several runs of kmens with different random seeds
#' @description
#' This function will select the optimal clustering result from several kmeans runs
#' with different random seeds, The clustering result with the smallest
#' within-cluster-sum of squared distances will be selected.
#' @param k Integer. Number of cluters to detect for kmeans.
#' @param x Matrix. This function will cluster the rows of the input matrix.
#' @param iter_max Integer. Number of iterations.
#' @param num_seeds Integer. Number of trials with random seeds
#' @param ... Further arguments handed to stats::kmeans
#' @return
#' Returns an object of class "kmeans" with is a list with several components
#' (see ?kmeans) which guves the local optimal kmeans clustering result within
#' #num_seeds trials.

optimal_km <- function (x,
                        k,
                        num_seeds = 10,
                        iter_max = 10,
                        ...){

  if(!is(x, 'matrix')){
    x = as.matrix(x)
  }

  if (num_seeds < 1){
    stop('num_seeds should be larger than 0.')
  }

  wcs <- Inf
  res <- NULL

  for (i in 1:num_seeds) {
    newres <- stats::kmeans(x = x,
                            centers = k,
                            iter.max = iter_max,
                            ...)

    newwcs <- sum(newres$withinss)

    if (newwcs <= wcs) {
      res <- newres
      wcs <- newwcs
    }
  }

  return(res)
}

#' Run spectral clustering
#'
#' @family biclustering
#' @description
#' Spectral clustering algorithm with normalized graph laplacian.
#'
#' @param caclust Caclust-class object.
#' @param dims Integer. Number of dimensions to choose from SVD of graph laplacian.
#' @param use_gap Logical, TRUE/FALSE. If TRUE, 'eigengap' method will be used to find the
#' most important eigenvector automatically, and the number of output clusters
#' equals number of selected eigenvectors. If FALSE, 'nclust'(integer) should be specified.
#' The eigenvectors corresponding with the smallest 'nclust' eigenvalues will be
#' selcted and 'nclust' clusters will be detected by skmeans/kmeans/GMM.
#' @param nclust Integer. Number of clusters.
#' @param spectral_method character. Name of the method to cluster the eigenvectors.
#' Can be on of the following 3:
#' * "kmeans": k-means clustering
#' * "skmeans": spherical k-means clustering
#' * "GMM": Gaussian-Mixture-Model fuzzy clustering.
#' @param iter_max Number of iterations for k-means clustering and GMM.
#' @param num_seeds Number of times k-means clustering is repeated.
#' @param return_eig Logical. Whether or not to return eigenvectors and store them in caclust-object.
#'
#' @return
#' The clustering results of type 'caclust'.
#' @md
#' @export
run_spectral <- function(caclust,
                         dims = 30,
                         use_gap = TRUE,
                         nclust = NULL,
                         spectral_method = 'kmeans',
                         iter_max = 10,
                         num_seeds = 10,
                         return_eig = TRUE) {

  call_params <- as.list(match.call())
  names(call_params)[1] <- "run_spectral"


  stopifnot(is(caclust, "caclust"))
  if (is.empty(caclust@SNN)){
    stop("No SNN graph found. Please run make_SNN() first!")
  }

  SNN <- caclust@SNN

  diag(SNN) = 0
  L = NormLaplacian(SNN)

  if (dims > ncol(SNN)){

    dims = ncol(SNN)

  }

  SVD <- irlba::irlba(L, nv = dims, smallest = TRUE) # eigenvalues in a decreasing order
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

  if (spectral_method == 'skmeans'){

    clusters = optimal_skm(eig, k = nclust, num_seeds = num_seeds)$cluster

  }else if (spectral_method == 'kmeans'){

    clusters = optimal_km(eig,
                         k = nclust,
                         iter_max = iter_max,
                         num_seeds = num_seeds)$cluster

  }else if (spectral_method == 'GMM'){

    gmm = ClusterR::GMM(eig,
                        gaussian_comps = ncol(eig),
                        dist_mode = "maha_dist",
                        seed_mode = "random_subset",
                        km_iter = iter_max,
                        em_iter = iter_max,
                        verbose = F)

    gmm_res = ClusterR::predict_GMM(data = eig,
                                 CENTROIDS = gmm$centroids,
                                 COVARIANCE = gmm$covariance_matrices,
                                 WEIGHTS = gmm$weights)

    gmm_res <- gmm_res[-which(names(gmm_res)=="log_likelihood")]

    cluster_proba <- gmm_res$cluster_proba
    clusters <- gmm_res$cluster_labels

    rownames(cluster_proba) <- rownames(SNN)
    colnames(cluster_proba) <- paste("BC", seq_len(ncol(cluster_proba)))

    # cell_idx <- which(rownames(cluster_proba) %in% rownames(caobj@prin_coords_cols))
    # gene_idx <- which(rownames(cluster_proba) %in% rownames(caobj@prin_coords_rows))

    cell_prob <- cluster_proba[caclust@cell_idxs,]
    gene_prob <- cluster_proba[caclust@gene_idxs,]

    caclust@cell_prob <- cell_prob
    caclust@gene_prob <- gene_prob



  }else{
  stop('clustering method should be chosen from kmeans and skmeans!')
  }

  # cell_idx <- which(names(clusters) %in% rownames(caobj@prin_coords_cols))
  # gene_idx <- which(names(clusters) %in% rownames(caobj@prin_coords_rows))

  clusters <- as.factor(clusters)
  names(clusters) <- rownames(SNN)

  cell_clusters <- clusters[caclust@cell_idxs]
  gene_clusters <- clusters[caclust@gene_idxs]

  caclust@cell_clusters <- cell_clusters
  caclust@gene_clusters <- gene_clusters

  if (return_eig){

    if (is.null(dims)){
      dims = min(30, ncol(SNN))
    }

    eigenv <- eigenvectors[,1:dims]

    rownames(eigenv) <- rownames(SNN)


  }else{
    eigenv <- matrix()
  }

  caclust@eigen <- eigenv

  stopifnot(validObject(caclust))

  return(caclust)
}




#' Leiden clustering on bigraph
#'
#' @family biclustering
#' @description
#' This function takes a caclust object with precomputed SNN-graph and
#' clusters cells and genes simultaneously.
#'
#' @param caclust caclust object with SNN calculated.
#' @param resolution float number. Resolution for leiden algorithm.
#' @param n.int Integer. Number of iterations for leiden algorithm.
#' @param rand_seed integer. Random seed.
#' @param cast_to_dense logical. Should the SNN-graph be converted to a dense
#' @param leiden_pack character. Optional values are 'igraph'(default) and 'leiden', the package used for leiden clustering.
#' matrix before running leiden clustering?
#' Casting to dense speeds up the leiden algorithm.
#'
#' @return
#' Object of type "caclust" with cell and gene clusters saved.
#' @export
run_leiden <- function(caclust,
                       resolution = 1,
                       n.int = 10,
                       rand_seed = 2358,
                       cast_to_dense = TRUE,
                       leiden_pack = 'igraph') {

  call_params <- as.list(match.call())
  names(call_params)[1] <- "run_leiden"

  stopifnot(is(caclust, "caclust"))

  if (is.empty(caclust@SNN)){
    stop("No SNN graph found. Please run make_SNN() first!")
  }

  if (leiden_pack == 'leiden'){
    if (is(caclust@SNN, "dgCMatrix") & isTRUE(cast_to_dense)){
      SNN <- as.matrix(caclust@SNN)
    }else{
      SNN <- caclust@SNN
    }


    suppressWarnings({
      ## due to the issue with packge leiden_0.4.3 and reticulate_1.26,
      ## we have to use suppressWarnings to ensure that the function runs smoothly.
    clusters <- leiden::leiden(object = SNN,
                               resolution_parameter = resolution,
                               partition_type = "RBConfigurationVertexPartition",
                               initial_membership = NULL,
                               weights = NULL,
                               node_sizes = NULL,
                               n_iterations = n.int,
                               seed = rand_seed)
    })
  }else if(leiden_pack == 'igraph'){

    SNN <- caclust@SNN

    g = igraph::graph_from_adjacency_matrix(
                                            SNN,
                                            mode = "undirected",
                                            weighted = TRUE,
                                            diag = TRUE,
                                            add.colnames = NULL,
                                            add.rownames = NA)
    set.seed(rand_seed)
    clusters = igraph::cluster_leiden(g,
                                      n_iterations = n.int,
                                      objective_function = 'modularity',
                                      weights = igraph::E(g)$weight,
                                      resolution_parameter =resolution)

    clusters = unlist(igraph::membership(clusters))
  }else{
    stop('leiden_pack should be either leiden or igraph!')
  }

  clusters <- as.factor(clusters)
  names(clusters) <- rownames(SNN)

  cell_clusters <- clusters[caclust@cell_idxs]
  gene_clusters <- clusters[caclust@gene_idxs]

  caclust@cell_clusters <- cell_clusters
  caclust@gene_clusters <- gene_clusters
  caclust@parameters <- append(caclust@parameters, call_params)


  stopifnot(validObject(caclust))
  return(caclust)

}



#' Run biclustering
#'
#' @description
#' Convenient wrapper around `make_SNN` and `run_leiden`/`run_spectral`.
#' `run_caclust` takes a cacomp object and biclusters cells and genes.
#'
#' @param caobj A cacomp object with principal and standard coordinates
#' calculated.
#' @param k Either an integer (same k for all subgraphs) or a vector of
#' exactly four integers specifying in this order: the k_c for the cell-cell
#' kNN-graph, k_g for the gene-gene kNN-graph, k_cg for the cell-gene
#' kNN-graph, k_gc for the gene-cell kNN-graph.
#' @param algorithm Character. Algorithm for clustering. Options are "leiden" or "spectral". Defalut: 'leiden'.
#' @inheritParams create_bigraph
#' @inheritParams make_SNN
#' @inheritParams run_leiden
#' @inheritParams run_spectral
#'
#' @return
#' Returns caclust objec with clustering results stored. The cell and gene
#' clusters can be accessed via `cell_clusters(obj)`/`gene_clusters(obj)`.
#'
#' @md
#' @export
run_caclust <- function(caobj,
                        k,
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
                        spectral_method = 'kmeans',
                        iter_max = 10,
                        num_seeds = 10,
                        return_eig = TRUE,
                        dims = NULL,
                        cast_to_dense = TRUE,
                        method = BiocNeighbors::KmknnParam(),
                        BPPARAM = BiocParallel::SerialParam(),
                        leiden_pack = 'leiden') {

  call_params <- as.list(match.call())
  names(call_params)[1] <- "Call"


  stopifnot("Invalid k! Should be either a single integer or of lenght 4!" =
              (length(k) == 1 | length(k) == 4))

  caclust <- make_SNN(caobj = caobj,
                      k = k,
                      loops = loops,
                      mode = mode,
                      SNN_prune = SNN_prune,
                      select_genes = select_genes,
                      prune_overlap = prune_overlap,
                      overlap = overlap,
                      calc_gene_cell_kNN = calc_gene_cell_kNN,
                      marker_genes = marker_genes)

  if (algorithm == "leiden"){

    caclust <- run_leiden(caclust = caclust,
                           resolution = resolution,
                           n.int = n.int,
                           rand_seed = rand_seed,
                           cast_to_dense = cast_to_dense,
                         leiden_pack = leiden_pack)

  } else if (algorithm == "spectral"){

    if (is.null(dims)){
      dims = length(caobj@D)
    }
    caclust <- run_spectral(caclust = caclust,
                               use_gap = use_gap,
                               nclust = nclust,
                               spectral_method = spectral_method,
                               iter_max = iter_max,
                               num_seeds = num_seeds,
                               return_eig = return_eig,
                               dims = dims)

  } else{
    stop("algorithm should choose from 'leiden' and 'spectral'!")
  }

  return(caclust)
}


#' Add caclust object to SingleCellExperiment object
#' @description
#' Add caclust clustering results to SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' @param caclust CAclust::caclust object
#' @param caclust_meta_name column name not listed in colData(sce),
#' rowData(sce), or metadata(sce)
#' @returns
#' SingleCellExperiment with caclust object stored.
#' @export
#'
add_caclust_sce <- function(sce, caclust, caclust_meta_name = 'caclust'){
  cell.clust <- cell_clusters(caclust)
  gene.clust <- gene_clusters(caclust)
  idx <- rownames(sce) %in% gene.clust
  matched_genes <- match(rownames(sce)[idx], names(gene.clust))

  SummarizedExperiment::colData(sce)[[caclust_meta_name]] <- cell.clust
  SummarizedExperiment::rowData(sce)[[caclust_meta_name]] <- 'not_in_caclust'
  SummarizedExperiment::rowData(sce)[[caclust_meta_name]][idx] <- gene.clust[matched_genes]

  S4Vectors::metadata(sce)[[caclust_meta_name]] <- caclust

  return(sce)
}


#' check_caobj_sce
#' @description
#' Check if cacomp object is already added to SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' @param cacomp_meta_name Character. Name of cacomp slpt in sce object.
#' @returns
#' TRUE if cacomp object is stored, FALSE otherwise.
#'
check_caobj_sce <- function(sce, cacomp_meta_name = 'CA'){

  ix <- cacomp_meta_name %in% SingleCellExperiment::reducedDimNames(sce)

  return(ix)

}



#' caclust
#' @family biclustering
#' @description
#' `caclust()` performs biclustering on either a "cacomp" or
#' "SingleCellExperiment" object.
#' @param obj A cacomp object or SingleCellExperiment object
#' @inheritParams run_caclust
#' @param ... further arguments
#' @details
#' Convenient wrapper around `make_SNN` and `run_leiden`/`run_spectral`.
#' `run_caclust` takes a cacomp object and biclusters cells and genes.
#' @return
#' A caclust object or SingleCellExperiment object
#' @export
setGeneric("caclust", function(obj,
                               k,
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
                               spectral_method = 'kmeans',
                               iter_max = 10,
                               num_seeds = 10,
                               return_eig = TRUE,
                               dims = NULL,
                               cast_to_dense = TRUE,
                               ...){
  standardGeneric("caclust")
})


#'
#' @rdname caclust
#' @export
setMethod(f = "caclust",
          signature(obj = "cacomp"),
          function(obj,
                   k,
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
                   spectral_method = 'kmeans',
                   iter_max = 10,
                   num_seeds = 10,
                   return_eig = TRUE,
                   dims = NULL,
                   cast_to_dense = TRUE,
                   method = BiocNeighbors::KmknnParam(),
                   BPPARAM = BiocParallel::SerialParam(),
                   leiden_pack = 'igraph',
                   ...){

          caclust_res <- run_caclust(caobj = obj,
                                     k = k,
                                     algorithm = algorithm,
                                     SNN_prune = SNN_prune,
                                     loops = loops,
                                     mode = mode,
                                     select_genes = select_genes,
                                     prune_overlap = prune_overlap,
                                     overlap =overlap,
                                     calc_gene_cell_kNN = calc_gene_cell_kNN,
                                     resolution = resolution ,
                                     marker_genes =marker_genes,
                                     n.int = n.int,
                                     rand_seed = rand_seed,
                                     use_gap = use_gap,
                                     nclust = nclust,
                                     spectral_method = spectral_method ,
                                     iter_max = iter_max,
                                     num_seeds = num_seeds,
                                     return_eig =return_eig,
                                     dims = dims,
                                     cast_to_dense = cast_to_dense,
                                     method = method,
                                     BPPARAM = BPPARAM,
                                     leiden_pack = leiden_pack,
                                     ...)
          return(caclust_res)

})


#'
#' @rdname caclust
#' @param cacomp_meta_name Character. The name of cacomp object stored in
#' metadata(SingleCellExperiment object). Default: 'caobj'.
#' @param caclust_meta_name the name of caclust object stored in
#' metadata(SingleCellExperiment object). Default: 'caclust.'
#' @export
setMethod(f = "caclust",
          signature(obj = "SingleCellExperiment"),
          function(obj,
                   k,
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
                   spectral_method = 'kmeans',
                   iter_max = 10,
                   num_seeds = 10,
                   return_eig = TRUE,
                   dims = NULL,
                   cast_to_dense = TRUE,
                   method = BiocNeighbors::KmknnParam(),
                   BPPARAM = BiocParallel::SerialParam(),
                   leiden_pack = 'igraph',
                   ...,
                   caclust_meta_name = 'caclust',
                   cacomp_meta_name = 'CA'){

            correct <- check_caobj_sce(obj, cacomp_meta_name = cacomp_meta_name)

            if(isFALSE(correct)){
              stop("No 'CA' dimension reduction object found. ",
                   "Please run cacomp(sce_obj, top, coords = FALSE, ",
                   "return_input=TRUE) first.")
            }


            if (isTRUE(caclust_meta_name %in% names(S4Vectors::metadata(obj)))){
              warning('The given meta_name or "caclust" is already in colData(sce)/rowData(sce)/metadata(sce), the slot will be overwritten!')
            }

            caobj <- APL::as.cacomp(obj)

            caclust_res <- run_caclust(caobj = caobj,
                                       k = k,
                                       algorithm = algorithm,
                                       SNN_prune = SNN_prune,
                                       loops = loops,
                                       mode = mode,
                                       select_genes = select_genes,
                                       prune_overlap = prune_overlap,
                                       overlap =overlap,
                                       calc_gene_cell_kNN = calc_gene_cell_kNN,
                                       resolution = resolution ,
                                       marker_genes =marker_genes,
                                       n.int = n.int,
                                       rand_seed = rand_seed,
                                       use_gap = use_gap,
                                       nclust = nclust,
                                       spectral_method = spectral_method ,
                                       iter_max = iter_max,
                                       num_seeds = num_seeds,
                                       return_eig =return_eig,
                                       dims = dims,
                                       cast_to_dense = cast_to_dense,
                                       method = method,
                                       BPPARAM = BPPARAM,
                                       leiden_pack = leiden_pack,
                                         ...)
            obj <- add_caclust_sce(sce = obj,
                                   caclust = caclust_res,
                                   caclust_meta_name = caclust_meta_name)

          return(obj)

})


#' Assign cluster to cells/genes
#' @description
#' Based on a probability cutoff genes/cells are assigned to all clusters for
#' which they have a probability higher than 'cutoff'.
#' @param caclust_obj A caclust object
#' @param type Either "cell" or "gene".
#' @param cutoff Probability cutoff.
#'
#' @returns
#' logical matrix indicating which gene belongs to which cluster.
#' @export
assign_clusters_GMM <- function(caclust_obj, type = "genes", cutoff=0.5){

  if(type == "genes"){
    prob_slot <- "gene_prob"
  } else if (type == "cells"){
    prob_slot <- "cells_prob"
  }

  probs <- slot(caclust_obj, prob_slot)
  stopifnot("No probabilities in caclust object." = !is.empty(probs))

  probs <- probs > cutoff

  return(probs)
}
