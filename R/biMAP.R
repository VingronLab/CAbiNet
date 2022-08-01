#' @include classes.R
NULL

#' Run UMAP embedding for cell-gene graph built up by caclust.
#' 
#' @description 
#' This function takes cacomp and caclust object as input to calculate UMAP embedding 
#' of cell-gene graph in several different ways:
#' * 'SNNdist'(Default): run UMAP on the distance matrix of cell-gene SNN graph built up by caclust, which is '1-adj(SNN)'.
#' * 'ca': run UMAP on the singular vectors from Correspondence Analysis.
#' * 'spectral': run UMAP on the selected eigenvectors of cell-gene graph laplacian (only eligiable when algorithm is set as 'spectral' in 'caclust' function)
#' @rdname run_biMAP
#' @param caobj A cacomp object with principal and standard coordinates 
#' calculated.
#' @param caclust_obj results from biclustering of class "caclust"
#' @param k_umap integer. Number of nearest neighbours to use from the SNN graph for
#' UMAP
#' @param rand_seed integer. Random seed for UMAP.
#' @param method Can be either "SNNgraph", "SNNdist", "spectral", "ca" or "ca_assR".
#'
#' @return 
#' caclust object with biMAP coordinates stored in the `bimap` slot.
#' 
#' @md
#' @export 
run_biMAP <- function(caclust_obj,
                      caobj = NULL,
                      k_umap = 30,
                      rand_seed = 2358,
                      method = 'SNNdist'){
  
  stopifnot(is(caclust_obj, "caclust"))
  stopifnot(method %in% c("SNNgraph", "SNNdist", "spectral", "ca", "ca_assR"))
  
   if (method == "SNNdist"){
    
    SNNdist <- as.matrix(1-get_snn(caclust_obj))
    
    reticulate::source_python(system.file("python/umap.py", package = "CAclust"), envir = globalenv())
    
    umap_coords <- python_umap(dm = SNNdist,
                                metric = "precomputed",
                                n_neighbors = as.integer(k_umap))
      
    umap_coords <- as.data.frame(umap_coords)
    rownames(umap_coords) <- colnames(SNNdist)
    
  }else if (method == 'spectral'){
    eigen = get_eigen(caclust_obj)
    
    if(is.na(eigen)) stop("Spectral clustering not run.")
    custom.config = umap::umap.defaults
    custom.config$random_state = rand_seed
    
    caclust_umap = umap::umap(eigen, 
                              config = custom.config,
                              n_neighbors = k_umap,
                              metric = 'cosine')
    
    umap_coords <- as.data.frame(caclust_umap$layout)
    
  }else if (method == 'ca'){
    stopifnot(!is.null(caobj))
    stopifnot(is(caobj, "cacomp"))
    
    SNN <- get_snn(caclust_obj)
    
    eigen = rbind(caobj@V, caobj@U)
    # eigen <- rbind(caobj@std_coords_cols, caobj@prin_coords_rows)
    custom.config = umap::umap.defaults
    custom.config$random_state = rand_seed
    
    eigen <- eigen[rownames(eigen) %in% rownames(SNN),]
    
    
    caclust_umap = umap::umap(eigen, 
                              config = custom.config,
                              metric = 'cosine',
                              n_neighbors = k_umap)
    umap_coords <- as.data.frame(caclust_umap$layout)
    
  }else if (method == 'ca_assR'){
    stopifnot(!is.null(caobj))
    stopifnot(is(caobj, "cacomp"))
    
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
    
    umap_coords <- python_umap(dm = assR,
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
  
  caclust_obj@bimap <- umap_coords
  return(caclust_obj)
}

#' Add cacomp obj results to SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' @param umap_coords data.frame with coordinates of genes and cells
#' @param  biMAP_meta_name name not listed in colData(sce), rowData(sce), or metadata(sce)
#' @export
#' 
add_biMAP_sce <- function(sce, umap_coords, biMAP_meta_name = 'biMAP'){
  
  S4Vectors::metadata(sce)[[biMAP_meta_name]] <- umap_coords
  
  return(sce)
}




#' biMAP
#' @description
#' This function takes cacomp and caclust objects or SingleCellExperiment object as input to calculate UMAP embedding 
#' of cell-gene graph in several different ways:
#' * 'SNNdist'(Default): run UMAP on the distance matrix of cell-gene SNN graph built up by caclust, which is '1-adj(SNN)'.
#' * 'ca': run UMAP on the singular vectors from Correspondence Analysis.
#' * 'spectral': run UMAP on the selected eigenvectors of cell-gene graph laplacian (only eligiable when algorithm is set as 'spectral' in 'caclust' function)
#' @name biMAP
#' @rdname biMAP
#' @param obj A cacomp object or SingleCellExperiment object  
#' @param caclust_obj the name of caclust object stored in metadata(SingleCellExperiment object)
#' @param cacomp_meta_name the name of cacomp object stored in metadata(SingleCellExperiment object)
#' @param caclust_meta_name the name of caclust object stored in metadata(SingleCellExperiment object)
#' @param biMAP_meta_name character. Slot name of biMAP coordinates in metadata of sce object. Default: biMAP_$method.
#' @inheritParams run_biMAP
#' @details
#' TODO
#' @return
#' an caclust object or SingleCellExperiment objects
#' @export
setGeneric("biMAP", function(obj,
                             caclust_obj,
                             k_umap,
                             cacomp_meta_name = 'CA',
                             caclust_meta_name = 'caclust',
                             biMAP_meta_name = NULL,
                             method = 'SNNdist',
                             message = TRUE,
                               ...){
  standardGeneric("biMAP")
})


#' @description
#' TODO
#' @rdname biMAP
#' @return
#' @export
setMethod(f = "biMAP",
          signature(obj = "cacomp", caclust_obj = "caclust"),
          function(obj, 
                   caclust_obj,
                   k_umap,
                   ...){
            
            umap_coords <- run_biMAP(caobj = obj,
                                     caclust_obj = caclust_obj,
                                     k_umap,
                                     ...)
            return(umap_coords)
            
          })

#' @description
#' TODO
#' @rdname biMAP
#' @inheritParams run_biMAP
#' @return
#' @export
setMethod(f = "biMAP",
          signature(obj = 'SingleCellExperiment'),
          function(obj, 
                   caclust_obj = NULL,
                   k_umap,
                   cacomp_meta_name = 'CA',
                   caclust_meta_name = 'caclust',
                   biMAP_meta_name = NULL,
                   method = 'SNNdist',
                   message = TRUE,
                   ...){
            
            if(is.null(biMAP_meta_name )){
              biMAP_meta_name = paste0('biMAP_', method)
            }
            
            check_caobj_sce(obj, cacomp_meta_name = cacomp_meta_name)
            if (isFALSE(caclust_meta_name %in% names(S4Vectors::metadata(obj)))){
              stop('The caclust_meta_name in not found in metadata(sce obj), change meta_name')
            }
            if (isTRUE(biMAP_meta_name %in% names(S4Vectors::metadata(obj)))){
              warning('The biMAP coordinates with name biMAP_meta_name already exist in metadata(sce obj), the old data will be overwritten!')
            }
            
            caobj <- S4Vectors::metadata(obj)[[cacomp_meta_name]]
            
            caclust_obj <- S4Vectors::metadata(obj)[[caclust_meta_name]]
            
            caclust_obj <- run_biMAP(caobj = caobj,
                                     caclust_obj = caclust_obj,
                                     k_umap,
                                     ...)
            
            S4Vectors::metadata(obj)[[caclust_meta_name]] <- caclust_obj
            #TODO
            #allow adding multi-bimap coordinate slots to caclust with slot names 'biMAP_'+algorithm, eg. 'biMAP_SNNdist'
            
            
            if(isTRUE(message)){
              message('bimap coordinates data.frame is added to metadata(sce obj) with name ', biMAP_meta_name, '.\n')
            }
            
            
            return(obj)
            
          })




