


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

#' An integration of several runs of skmens with different random seeds
#' @description 
#' This function will select the optimal clustering result from several 
#' skmeans::skmeans runs
#' with different random seeds, The clustering result with smallest 
#' within-cluster-sum
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
                        spectral_method = 'kmeans',
                        iter.max=10, 
                        num.seeds=10,
                        return.eig = TRUE,
                        sc.dims = NULL,
                        leiden.dense = TRUE) {
  
  call_params <- as.list(match.call())
  names(call_params)[1] <- "Call"
  
  distances <- calc_distances(caobj = caobj)
  
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



#' caclust
#' @description
#' run biclustering algorithm, compatible for matrix and sce input
#' @name caclust
#' @rdname caclust
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


          