#' @include classes.R
NULL



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



#' An integration of several runs of skmens with different random seeds
#' @description 
#' This function will select the optimal clustering result from several 
#' skmeans::skmeans runs
#' with different random seeds, The clustering result with the smallest 
#' within-cluster-sum of squared distances will be selected.
#' @inheritParams skmeans::skmeans
#' @param num_seeds number of trials with random seeds
#' 
#' @return 
#' The optimal skmeans clustering result
#' TODO return
SKMeans <- function (x, 
                     k, 
                     method = NULL, 
                     m = 1, 
                     weights = 1, 
                     control = list(), 
                     num_seeds = 10){
  
  # if (mode(x) == "numeric")  x <- data.frame(new.x = x)
  
  ###############
  ## IMPORTANT ##
  ###############
  # This code seems to be heavily inspired by RcmdrMisc::KMeans
  # https://github.com/cran/RcmdrMisc/blob/master/R/cluster.R
  # Do we need a citation? 
  
  KM <- skmeans::skmeans(x = x, 
                         k = k, 
                         method = method,
                         m = m,
                         weights = weights,
                         control = control)
  wss <- do.call(sum, lapply(list(1:length(KM$cluster)), function(i){
    sum((x[i,]-KM$prototypes[KM$cluster[i],])^2)
  }))
  
  for (i in 2:num_seeds) {
    newKM <- skmeans::skmeans(x = x, 
                              k = k, 
                              method = method,
                              m = m,
                              weights = weights,
                              control = control)
    
    newwss <- do.call(sum, lapply(list(1:length(newKM$cluster)), function(i){
      sum((x[i,]-newKM$prototypes[newKM$cluster[i],])^2)
    }))
    if (newwss <wss) {
      KM <- newKM
      wss <- newwss
    }
  }
  # xmean <- apply(x, 2, mean)
  # centers <- rbind(KM$centers, xmean)
  # bss1 <- as.matrix(dist(centers)^2)
  # KM$betweenss <- sum(as.vector(bss1[nrow(bss1), ]) * c(KM$size, 
  #                                                       0))
  return(KM)
}

#' Run spectral clustering
#'
#' @description 
#' This function is designed for detecting clusters from input graph adjacency 
#' matrix by using spectral clustering with normalized graph laplacian.
#' 
#' @param caclust caclust object.
#' @param dims integer. Number of dimensions to compute during SVD for spectral clustering.
#' @param use_gap TRUE/FALSE. If TRUE, 'eigengap' method will be used to find the
#' most important eigenvector automatically, and the number of output clusters 
#' equals number of selected eigenvectors. If FALSE, 'nclust'(integer) should be given. 
#' The eigenvectors corresponding with the smallest 'nclust' eigenvalues will be
#' selcted and 'nclust' clusters will be detected by skmeans.
#' @param spectral_method character. Name of the method to cluster the eigenvectors.
#' Can be on of the following 3:
#' * "kmeans": k-means clustering
#' * "skmeans": spherical k-means clustering
#' * "GMM": Gaussian-Mixture-Model fuzzy clustering.
#' @param iter_max Number of iterations for k-means clustering and GMM.
#' @param num_seeds Number of times k-means clustering is repeated.
#' @param return_eig Whether or not to return eigenvectors.
#' 
#' @return 
#' The clustering results
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
    
    clusters = SKMeans(eig, k = nclust, num_seeds = num_seeds)$cluster
  
  }else if (spectral_method == 'kmeans'){
    
    clusters = RcmdrMisc::KMeans(eig,
                                 centers = ncol(eig),
                                 iter.max = iter_max,
                                 num.seeds = num_seeds)$cluster

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
#' @description 
#' This function takes a caclust object with precomputed SNN-graph and 
#' clusters cells and genes simultaneously.
#' 
#' @param caclust caclust object with SNN calculated.
#' @param resolution integer. Resolution for leiden algorithm.
#' @param n.int integer. Number of iterations for leiden algorithm.
#' @param rand_seed Random seed.
#' @param cast_to_dense logical. Should the SNN-graph be converted to a dense 
#' matrix before running leiden clustering?
#' Casting to dense speeds up the leiden algorithm.
#' 
#' @return 
#' Object of type "caclust" with cell and gene clusters saved.
#' 
run_leiden <- function(caclust,
                       resolution = 1,
                       n.int = 10, 
                       rand_seed = 2358,
                       cast_to_dense = TRUE) {
  
  call_params <- as.list(match.call())
  names(call_params)[1] <- "run_leiden"
  
  stopifnot(is(caclust, "caclust"))
  
  if (is.empty(caclust@SNN)){
    stop("No SNN graph found. Please run make_SNN() first!")
  }
  
  if (is(caclust@SNN, "dgCMatrix") & isTRUE(cast_to_dense)){
    SNN <- as.matrix(caclust@SNN)
  } else {
    SNN <- caclust@SNN
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
  
  # 
  # cell_idx <- which(names(clusters) %in% rownames(caobj@prin_coords_cols))
  # gene_idx <- which(names(clusters) %in% rownames(caobj@prin_coords_rows))
  
  cell_clusters <- clusters[caclust@cell_idxs]
  gene_clusters <- clusters[caclust@gene_idxs]
  
  caclust@cell_clusters <- cell_clusters
  caclust@gene_clusters <- gene_clusters
  caclust@parameters <- append(caclust@parameters, call_params)
  
  # caclust_res <- do.call(new_caclust, list("cell_clusters" = cell_clusters,
  #                                          "gene_clusters" = gene_clusters,
  #                                          "SNN" = SNN,
  #                                          "eigen" = matrix(),
  #                                          "parameters" = call_params))
  # 
  
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
#' @param algorithm Algorithm for clustering. Options are "leiden" or "spectral".
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
                        cast_to_dense = TRUE) {
  
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
                           cast_to_dense = cast_to_dense)
    
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


#' Add caclust clustering results to SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' @param caclust caclust::caclust object
#' @param caclust_meta_name column name not listed in colData(sce), rowData(sce), or metadata(sce)
#' @export
#' 
add_caclust_sce <- function(sce, caclust, caclust_meta_name = 'caclust'){
  cell.clust <- cell_clusters(caclust)
  gene.clust <- gene_clusters(caclust)
  idx <- rownames(sce) %in% gene.clust
  matched_genes <- match(rownames(sce)[idx], names(gene.clust))
  
  colData(sce)[[caclust_meta_name]] <- cell.clust
  rowData(sce)[[caclust_meta_name]] <- 'not_in_caclust'
  rowData(sce)[[caclust_meta_name]][idx] <- gene.clust[matched_genes]
  
  metadata(sce)[[caclust_meta_name]] <- caclust
  
  return(sce)
}

#' Add cacomp obj results to SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' @param caobj caclust::caclust object
#' @param cacomp_meta_name column name not listed in colData(sce), rowData(sce), or metadata(sce)
#' @export
#' 
add_caobj_sce <- function(sce, caobj, cacomp_meta_name = 'caobj'){
  
  if (isTRUE(cacomp_meta_name %in% colnames(colData(sce))) | 
      isTRUE(cacomp_meta_name %in% colnames(rowData(sce))) |
      isTRUE(cacomp_meta_name %in% names(metadata(sce)))){
    stop('The given cacomp_meta_name is already in colData(sce)/rowData(sce)/metadata(sce), change meta_name')
  }
  
  metadata(sce)[[cacomp_meta_name]] <- caobj
  
  return(sce)
}

#' check if cacomp object is already added to SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' 
check_caobj_sce <- function(sce, cacomp_meta_name = 'caobj'){
  
 ix <- 'cacomp' %in% unlist(lapply(metadata(sce), class))
 if(isFALSE(ix)){
   stop('cacomp object is missing from SingleCellExperiment object sce, plese make sure you have
        run cacomp and add_caobj_sce in ahead!')
 }
 
 ix <- cacomp_meta_name %in% names(metadata(sce))
 if(isFALSE(ix)){
   stop('cacomp object with the given name is missing from SingleCellExperiment object sce, plese make 
   sure you have run add_caobj_sce with the same cacomp_meta_name!')
 }
  
}



#'
#' @description
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#' @param obj A cacomp object or SingleCellExperiment object  
#' @param cacomp_meta_name the name of cacomp object stored in metadata(SingleCellExperiment object)
#' @param caclust_meta_name the name of caclust object stored in metadata(SingleCellExperiment object)
#' @inheritParams run_caclust
#' @details
#' TODO
#' @return
#' an caclust object or SingleCellExperiment objects
#' @export
setGeneric("caclust", function(obj,
                               k_sym,
                               k_asym = k_sym,
                               cacomp_meta_name = 'caobj',
                               caclust_meta_name = 'caclust',
                               ...){
  standardGeneric("caclust")
})


#
#' @rdname caclust
#' @inheritParams run_caclust
#' @export
setMethod(f = "caclust",
          signature(obj = "caclust"),
          function(obj, 
                   k_sym,
                   ...){
          
          caclust_res <- run_caclust(caobj = obj,
                                 k_sym,
                                 ...)
          return(caclust_res)
            
})

#
#' @rdname caclust
#' @param obj SingleCellExperiment object
#' @param k_sym neighour of nearest neighours, see run_caclust  function
#' @param cacomp_meta_name the name of cacomp object stored in metadata(SingleCellExperiment object)
#' @param caclust_meta_name the name of caclust object stored in metadata(SingleCellExperiment object)
#' @export
setMethod(f = "caclust",
          signature(obj = "SingleCellExperiment"),
          function(obj, 
                   k_sym,
                   cacomp_meta_name = 'caobj',
                   caclust_meta_name = 'caclust',
                   ...){
          
            check_caobj_sce(obj, cacomp_meta_name = cacomp_meta_name)
            if (isTRUE(caclust_meta_name %in% names(metadata(sce)))){
              stop('The given meta_name or "caclust" is already in colData(sce)/rowData(sce)/metadata(sce), change meta_name')
            }
            
            caobj <- metadata(obj)[[cacomp_meta_name]]
            
            caclust_res <- run_caclust(caobj = caobj,
                                         k_sym,
                                         ...)
            obj <- add_caclust_sce(sce = obj, 
                                   caclust = caclust_res,
                                   caclust_meta_name = caclust_meta_name)
          
          return(obj)
          
})


          