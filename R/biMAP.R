
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
run_biMAP <- function(caobj,
                      caclust_obj,
                      k_umap,
                      rand_seed = 2358,
                      algorithm = 'SNNdist'){
  
  stopifnot(is(caobj, "cacomp"))
  stopifnot(is(caclust_obj, "caclust"))
  stopifnot(algorithm %in% c("SNNgraph", "SNNdist", "spectral", "ca", "ca_assR"))
  
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
    
    custom.config = umap::umap.defaults
    custom.config$random_state = rand_seed
    
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
    
    reticulate::source_python(system.file("python/umap.py", package = "CAclust"), envir = globalenv())
    
    umap_coords = python_umap(dm = SNNdist,
                              metric = "precomputed",
                              n_neighbors = as.integer(k_umap))
    
    umap_coords <- as.data.frame(umap_coords)
    rownames(umap_coords) <- colnames(SNNdist)
    
  }else if (algorithm == 'spectral'){
    
    eigen = get_eigen(caclust_obj)
    
    if(is.na(eigen)) stop("Spectral clustering not run.")
    custom.config = umap::umap.defaults
    custom.config$random_state = rand_seed
    
    caclust_umap = umap::umap(eigen, 
                              config = custom.config,
                              n_neighbors = k_umap,
                              metric = 'cosine')
    
    umap_coords <- as.data.frame(caclust_umap$layout)
    
  }else if (algorithm == 'ca'){
    
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

