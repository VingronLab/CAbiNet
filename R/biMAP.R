#' @include classes.R
NULL

#' Run UMAP embedding for cell-gene graph built up by caclust.
#' 
#' @description 
#' This function takes cacomp and caclust object as input to calculate UMAP embedding 
#' of cell-gene graph in several different ways:
#' * 'SNNdist'(Default): run UMAP on the distance matrix of cell-gene SNN graph
#' built up by caclust, which is '1-adj(SNN)'.
#' * 'spectral': run UMAP on the selected eigenvectors of cell-gene graph 
#' laplacian (only eligiable when algorithm is set as 'spectral' in 'caclust' 
#' function)
#' * 'ca': run UMAP on the singular vectors from Correspondence Analysis.
#' @rdname run_biMAP
#' @param obj results from biclustering of class "caclust"
#' @param caobj A cacomp object with principal and standard coordinates 
#' calculated. Only needs to be supplied when using method "ca". 
#' @param k integer. Number of nearest neighbours to use to compute UMAP.
#' @param rand_seed integer. Random seed for UMAP.
#' @param method Can be either "SNNdist", "spectral" or "ca". When using "ca",
#' a "cacomp" object has to be provided for `caobj`.
#'
#' @return 
#' caclust object with biMAP coordinates stored in the `bimap` slot.
#' 
#' @md
run_biMAP <- function(obj,
                      caobj = NULL,
                      k = 30,
                      rand_seed = 2358,
                      method = 'SNNdist'){
  
  stopifnot(is(obj, "caclust"))
  stopifnot(method %in% c("SNNdist", "spectral", "ca"))
  
  if (method == "SNNdist"){
    
    SNNdist <- as.matrix(1-get_snn(obj))
    
    reticulate::source_python(system.file("python/umap.py", package = "CAclust"), envir = globalenv())
    
    umap_coords <- python_umap(dm = SNNdist,
                                metric = "precomputed",
                                n_neighbors = as.integer(k))
      
    umap_coords <- as.data.frame(umap_coords)
    rownames(umap_coords) <- colnames(SNNdist)
    
  }else if (method == 'spectral'){
    eigen = get_eigen(obj)
    
    if(is.na(eigen)) stop("Spectral clustering not run.")
    custom.config = umap::umap.defaults
    custom.config$random_state = rand_seed
    
    caclust_umap = umap::umap(eigen, 
                              config = custom.config,
                              n_neighbors = k,
                              metric = 'cosine')
    
    umap_coords <- as.data.frame(caclust_umap$layout)
    
  }else if (method == 'ca'){
    
    stopifnot(!is.null(caobj))
    stopifnot(is(caobj, "cacomp"))
    
    SNN <- get_snn(obj)
    
    eigen = rbind(caobj@V, caobj@U)
    # eigen <- rbind(caobj@std_coords_cols, caobj@prin_coords_rows)
    custom.config = umap::umap.defaults
    custom.config$random_state = rand_seed
    
    eigen <- eigen[rownames(eigen) %in% rownames(SNN),]
    
    
    caclust_umap = umap::umap(eigen, 
                              config = custom.config,
                              metric = 'cosine',
                              n_neighbors = k)
    umap_coords <- as.data.frame(caclust_umap$layout)
    
  } else {
    stop()
  }
  
  
  cellc <- cell_clusters(obj)
  genec <- gene_clusters(obj)
  
  
  colnames(umap_coords) <- c("x", "y")
  umap_coords$name <- rownames(umap_coords)
  
  umap_coords$type <- NA
  umap_coords$type[umap_coords$name %in% names(cellc)] <- "cell" 
  umap_coords$type[umap_coords$name %in% names(genec)] <- "gene" 
  
  umap_coords$cluster <- NA
  cell_idx <- na.omit(match(names(cellc), umap_coords$name))
  gene_idx <- na.omit(match(names(genec), umap_coords$name))
  
  umap_coords$cluster[cell_idx] <- cell_clusters(obj)
  umap_coords$cluster[gene_idx] <- gene_clusters(obj)
  umap_coords$cluster <- as.factor(umap_coords$cluster)
  
  umap_coords <- umap_coords %>% dplyr::arrange(desc(type))
  
  obj@bimap <- umap_coords
  return(obj)
}

# #' Add cacomp obj results to SingleCellExperiment object
# #' @param sce SingleCellExperiment object
# #' @param umap_coords data.frame with coordinates of genes and cells
# #' @param  biMAP_meta_name name not listed in colData(sce), rowData(sce), or metadata(sce)
# #' @export
# #' 
# add_biMAP_sce <- function(sce, umap_coords, biMAP_meta_name = 'biMAP'){
#   
#   S4Vectors::metadata(sce)[[biMAP_meta_name]] <- umap_coords
#   
#   return(sce)
# }




#' Compute biMAP
#' 
#' @description
#' The function takes either a `caclust` or `SingleCellExperiment`as input and 
#' stores the biMAP in the "bimap" slot in the caclust object. If a
#' SingleCellExperiment was provided the caclust object is stored in its 
#' metadata.
#' @name biMAP
#' @rdname biMAP
#' @param obj A caclust object or SingleCellExperiment object 
#' @param method Can be either "SNNdist" or "spectral".
#' @inheritParams run_biMAP
#' @param ... Further arguments
#' @details
#' The biMAP cell and gene embeddings can be calculated via different methods 
#' as controlled by the parameter `method`:
#' * 'SNNdist'(Default): run UMAP on the distance matrix of cell-gene SNN graph
#' built up by caclust, which is '1-adj(SNN)'.
#' * 'spectral': run UMAP on the selected eigenvectors of cell-gene graph 
#' laplacian (only eligiable when algorithm is set as 'spectral' in 'caclust' 
#' function)
#' @return
#' A caclust object or SingleCellExperiment object.
#' 
#' @md
#' @export
setGeneric("biMAP", function(obj,
                             k = 30,
                             rand_seed = 2358,
                             method = 'SNNdist',
                             ...){
  standardGeneric("biMAP")
})



#' @rdname biMAP
#' @export
setMethod(f = "biMAP",
          signature(obj = "caclust"),
          function(obj,
                   k = 30,
                   rand_seed = 2358,
                   method = 'SNNdist',
                   ...){
            
            stopifnot(method %in% c("SNNdist", "spectral"))
            
            obj <- run_biMAP(obj = obj,
                             caobj = NULL,
                             k = k,
                             rand_seed = rand_seed,
                             method = method)
            return(obj)
            
          })

#' @rdname biMAP
#' @param caclust_meta_name The name of caclust object stored in metadata
#' (SingleCellExperiment object).
#' @export
setMethod(f = "biMAP",
          signature(obj = 'SingleCellExperiment'),
          function(obj,
                   k = 30,
                   rand_seed = 2358,
                   method = 'SNNdist',
                   caclust_meta_name = 'caclust',
                   ...){
            
            stopifnot(method %in% c("SNNdist", "spectral"))
            
            caclust_obj <- S4Vectors::metadata(obj)[[caclust_meta_name]]
            
            caclust_obj <- run_biMAP(obj = caclust_obj,
                                     caobj = NULL,
                                     k = k,
                                     rand_seed = rand_seed,
                                     method = method)
            
            S4Vectors::metadata(obj)[[caclust_meta_name]] <- caclust_obj
            # TODO
            # allow adding multi-bimap coordinate slots to caclust with slot 
            # names 'biMAP_'+algorithm, eg. 'biMAP_SNNdist'
            
            
            return(obj)
            
          })



#' Compute biMAP basedon CA results
#'
#' @description
#' The function takes either a `caclust` and a `cacomp` object as input and 
#' computes the biMAP embedings for cells and genes on the basis of the singular
#' vectors of CA.
#' @name ca_biMAP
#' @rdname ca_biMAP
#'
#' @param obj A `caclust` object or `SingleCellExperiment` object with `caclust`
#' and `cacomp` objects stored in the metadata.
#' @param caobj A `cacomp` object.
#' @inheritParams run_biMAP
#' @param ... Further arguments
#' @details
#' The biMAP embeddings are computed on the basis of the singular vectors 
#' from Correspondence Analysis.
#' If a `SingleCellExperiment` object with `caclust` and `cacomp`
#' objects stored is provided the argument `caobj` is not required.
#'
#' @return
#' an caclust object or SingleCellExperiment objects
#' @export
setGeneric("ca_biMAP", function(obj,
                                caobj,
                                k = 30,
                                rand_seed = 2358,
                                ...){
  standardGeneric("ca_biMAP")
})



#' @rdname ca_biMAP
#' @export
setMethod(f = "ca_biMAP",
          signature(obj = "caclust", caobj = "cacomp"),
          function(obj,
                   caobj,
                   k = 30,
                   rand_seed = 2358,
                   ...){
            
            obj <- run_biMAP(obj = obj,
                             caobj = caobj,
                             k = 30,
                             rand_seed = 2358,
                             method = 'ca')
            return(obj)
            
          })



#' @rdname ca_biMAP
#' @param cacomp_meta_name the name of cacomp object stored in 
#' metadata(obj)
#' @param caclust_meta_name the name of caclust object stored in metadata(obj)
#' @export
setMethod(f = "ca_biMAP",
          signature(obj = 'SingleCellExperiment'),
          function(obj,
                   caobj = NULL,
                   k = 30,
                   rand_seed = 2358,
                   caclust_meta_name = "caclust",
                   cacomp_meta_name = "cacomp",
                   ...){
            
            
            caobj <- S4Vectors::metadata(obj)[[cacomp_meta_name]]
            
            caclust_obj <- S4Vectors::metadata(obj)[[caclust_meta_name]]
            
            caclust_obj <- run_biMAP(obj = caclust_obj,
                                     caobj = caobj,
                                     k = k,
                                     rand_seed = rand_seed,
                                     method = "ca")
            
            S4Vectors::metadata(obj)[[caclust_meta_name]] <- caclust_obj
            # TODO
            # allow adding multi-bimap coordinate slots to caclust with slot 
            # names 'biMAP_'+algorithm, eg. 'biMAP_SNNdist'
            
            
            return(obj)
            
          })

