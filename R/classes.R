
#' Checks if caclust-class was constructed correctly
#' 
#' @param object caclust object
#' 
#' @return 
#' If object is a valid caclust object returns TRUE, otherwise the errors.
#' 
check_caclust <- function(object){
  
  errors <- character()
  
  if (!is(object, "caclust")){
    msg <- "'object' is not of class caclust!"
    errors <- c(errors, msg)
  }
  
#  if(!identical(levels(object@cell_clusters), levels(object@gene_clusters))){
#    msg <- "factor levels for cells and genes are not the same!"
#    errors <- c(errors, msg)
#  }
  
  if(isTRUE(!is.empty(object@SNN)) & nrow(object@SNN) != ncol(object@SNN)){
    msg <- "SNN number of rows not equal to number of columns!"
    errors <- c(errors, msg)
  }
  
  if(isTRUE(!is.empty(object@cell_clusters)) & 
            is.null(names(object@cell_clusters))){
    msg <- "cell clusters are not named!"
    errors <- c(errors, msg)
  }
  
  if(isTRUE(!is.empty(object@gene_clusters)) & 
     is.null(names(object@gene_clusters))){
    msg <- "gene clusters are not named!"
    errors <- c(errors, msg)
  }
  
  if(isTRUE(!is.empty(object@SNN)) & 
     !identical(rownames(object@SNN), colnames(object@SNN))){
    msg <- "SNN rownames not identical to colnames!"
    errors <- c(errors, msg)
  }
  
  if (isTRUE(!is.empty(object@SNN)) &
      isTRUE(!is.empty(object@cell_clusters)) &
      isTRUE(!is.empty(object@gene_clusters))){
    
    if(!all(c(names(object@cell_clusters),
              names(object@gene_clusters)) %in% rownames(object@SNN))){
      msg <- "not all cell & gene names in SNN!"
      errors <- c(errors, msg)
    }
    
    if(!all(rownames(object@SNN) %in% c(names(object@cell_clusters),
                                        names(object@gene_clusters)))){
      msg <- "not all SNN names in cell & gene names!"
      errors <- c(errors, msg)
    }
  }
  
  
  if (length(errors) == 0) TRUE else errors
}

#' An S4 class that contains all elements needed for CA.
#' @name caclust-class
#' @rdname caclust-class
#' @description
#' Class to store biclustering clustering results.
#'
#' @slot cell_clusters factors. The assigned cell clusters with cell names in 
#' the names attribute.
#' @slot gene_clusters factors. The assigned gene clusters with gene names in 
#' the names attribute.
#' @slot SNN sparse shared nearest neighbours matrix. Values indicate the
#' jaccard similarity.
#' @slot parameters List of used parameters and function name with which results
#' were generated.
#' @slot eigen matrix, Slot for storing eigenvectors from spectral clustering
#' @slot cell_prob matrix. Matrix that stores the probabilities that a cell belongs
#' to a cluster. Only filled when running spectral clustering with GMM.
#' @slot gene_prob matrix. Matrix that stores the probabilities that a gene belongs
#' to a cluster. Only filled when running spectral clustering with GMM.
#' @slot cell_idxs integer. Indices of the cells in the SNN adjacency matrix.
#' @slot gene_idxs integer. Indices of the genes in the SNN adjacency matrix.
#' @slot bimap data.frame. Data frame storing the biMAP coordinates (x, y) and 
#' the type (cell or gene) as well as the assigned clusters.
#' 
#' @export
setClass("caclust",
         representation(
           cell_clusters = "factor",
           gene_clusters = "factor",
           SNN = "dgCMatrix",
           eigen = 'matrix',
           parameters = "list",
           cell_prob = "matrix",
           gene_prob = "matrix",
           cell_idxs = "integer",
           gene_idxs = "integer",
           bimap = "data.frame"
         ),
         prototype(
           cell_clusters = factor(),
           gene_clusters = factor(),
           SNN = as(matrix(0,0,0), "dgCMatrix"),
           eigen = matrix(0,0,0),
           parameters = list(),
           cell_prob = matrix(0,0,0),
           gene_prob = matrix(0,0,0),
           cell_idxs = NA_integer_,
           gene_idxs = NA_integer_,
           bimap = data.frame()),
         # contains = "cacomp",
         validity = check_caclust
)



#' Create new caclust object
#' 
#' @param ... arguments forwarded to new. Should be named slots of a caclust
#' object
#' @description 
#' TODO
new_caclust <- function(...) new("caclust",...)


#' Try converting to caclust object from list or biclustlib.
#'
#' @description
#' Converts object to caclust object if possible.
#'
#' @details
#' TODO
#'
#' @return
#' A caclust object.
#'
#' @param obj object
#' @export
setGeneric("as.caclust", function(obj) {
  standardGeneric("as.caclust")
})


#' Convert from list to caclust object.
#' @param obj A list with elements named same as caclust slots.
#' @rdname as.caclust
#' 
#' @export
setMethod(f = "as.caclust",
          signature=(obj="list"),
          function(obj){
  
  obj <- do.call(new_caclust, obj)         
  return(obj)
})

#' Get cell clusters
#' @param object caclust object
#' @export
cell_clusters <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "cell_clusters")
}

#' Get gene clusters
#' @param object caclust object
#' @export
gene_clusters <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "gene_clusters")
}

#' Get SNN graph
#' @param object caclust object
#' @export
get_snn <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "SNN")
}

#' Get eigenvectors from spectral clustering
#' @param object caclust object
#' @export
get_eigen <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "eigen")
}

#' Get parameters for generating clustering
#' @param object caclust object
#' @export
get_params <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "parameters")
}

#' Get probablities of cell clusters
#' @param object caclust object
#' @export
get_cell_prob <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "cell_prob")
}

#' Get probablities of gene clusters
#' @param object caclust object
#' @export
get_gene_prob <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "gene_prob")
}


#' Get biMAP coordinates
#' @param object caclust object
#' @export
get_bimap <- function(object){
  stopifnot(is(object, "caclust"))
  slot(object, "bimap")
}


#' Print caclust object in console
#' @param object a caclust object
show.caclust <- function(object){
  
  stopifnot(is(object, "caclust"))
  
  ncells <- length(object@cell_idxs) 
  ngenes <- length(object@gene_idxs) 
  cat("caclust object with", ncells, "cells and", ngenes, "genes.")
  
  if (!is.empty(cell_clusters(object)) & !is.empty(gene_clusters(object))){
    stopifnot(identical(levels(cell_clusters(object)),
                        levels(gene_clusters(object))))
    cat("\n")
    cat(length(levels(gene_clusters(object))), "clusters found.")
    cat("\nClustering results:\n\n")
    df <- data.frame("cluster" = levels(cell_clusters(object)),
                     "ncells" = summary(cell_clusters(object), maxsum = Inf),
                     "ngenes" = summary(gene_clusters(object), maxsum = Inf))
    
    print(df, row.names = FALSE, right = FALSE)
    
  } else {
    cat("\nNo biclustering run yet.\n\n")
  }

  # print(df[seq_len(min(n_rows, nrows(df))),], row.names = FALSE, right = FALSE)
}

#' @rdname show.caclust
#' @export
setMethod(f = "show",
          signature(object = "caclust"),
          function(object) {
  show.caclust(object)
})


#' Converts CAclust results to biclustlib results
#'
#'@param caclust A caclust object
#'
#'@return 
#' An object of type "Biclust".
#'@export
convert_to_biclust <- function(caclust){
  
  stopifnot(is(caclust, "caclust"))
  
  cc <- cell_clusters(caclust)
  gc <- gene_clusters(caclust)
  params <- get_params(caclust)
  
  ctypes = sort(unique(cc))
  gtypes = sort(unique(gc))
  bitypes = union(ctypes, gtypes)
  
  Number = length(bitypes)

  if (Number == 0){
    NumberxCol = matrix(0)
    RowxNumber = matrix(0)
  }else{
    NumberxCol = do.call(rbind, lapply(bitypes, function(x){cc == x}))
    RowxNumber = do.call(cbind, lapply(bitypes, function(x){gc == x}))
  }
  
  rownames(RowxNumber) <- names(gc)
  colnames(RowxNumber) <- paste0("BC", bitypes)

  rownames(NumberxCol) <- paste0("BC", bitypes)
  colnames(NumberxCol) <- names(cc)

  bic <- new("Biclust", "Parameters" = params,
                       "RowxNumber" = RowxNumber,
                       "NumberxCol" = NumberxCol,
                       "Number" = Number,
                       "info" = list("Output of CAclust package"))
  
  return(bic)
}


#' Helper function to check if object is empty.
#' @param x object
#' @return TRUE if x has length 0 and is not NULL. FALSE otherwise
is.empty <- function(x) return(isTRUE(length(x) == 0 & !is.null(x)))

#' Remove clusters only consisting of cells/genes
#' 
#' @description 
#' Takes an object of class biclust and removes all clusters that only consist
#' of cells or genes.
#' 
#' @param obj biclustering results from caclust or biclust or a SingleCellExperiment object.
#' @param ... further arguments
#' @return 
#' caclust/biclust object with monoclusters removed.
#' 
#' @export
setGeneric("rm_monoclusters", function(obj,
                                       ...) {
  standardGeneric("rm_monoclusters")
})


#' @rdname rm_monoclusters
#' @export
setMethod(f = "rm_monoclusters",
          signature=(obj="caclust"),
          function(obj){
        
        cc <- unique(cell_clusters(obj))
        gc <- unique(gene_clusters(obj))
        
        keep <- sort(as.character(intersect(cc,gc)))
        
        if(length(keep > 0)){
          
          obj@cell_clusters <- obj@cell_clusters[as.character(obj@cell_clusters) %in% keep]
          obj@cell_clusters <- droplevels(obj@cell_clusters)

          obj@gene_clusters <- obj@gene_clusters[as.character(obj@gene_clusters) %in% keep]
          obj@gene_clusters <- droplevels(obj@gene_clusters)

        }
        
        stopifnot(!is.null(names(obj@cell_clusters)))
        stopifnot(!is.null(names(obj@gene_clusters)))
        
        if(!is.empty(obj@SNN)){
          
          selr <- which(rownames(obj@SNN) %in% c(names(obj@cell_clusters),
                                                 names(obj@gene_clusters)))
          
          selc <- which(colnames(obj@SNN) %in% c(names(obj@cell_clusters),
                                                 names(obj@gene_clusters)))
          obj@SNN <- obj@SNN[selr, selc]
          
        }
        
        if(!is.empty(obj@bimap)){
         warning("The biMAP embedding is not valid after removing mono-clusters",
                 "and should be recomputed!") 
        }

        # Current implementation does not subset biMAP coords.
        # if(!is.empty(obj@bimap)){
        # 
        #   sel <- which(obj@bimap$cluster %in% keep)
        # 
        #   obj@bimap <- obj@bimap[sel,]
        #   
        # }

        return(obj)      
                       
})


#' @rdname rm_monoclusters
#' @export
setMethod(f = "rm_monoclusters",
          signature=(obj="Biclust"),
          function(obj){
            
            keep <- colSums(obj@RowxNumber) > 0 & rowSums(obj@NumberxCol) > 0
            
            if(any(!keep)){
              
              obj@RowxNumber <- obj@RowxNumber[,keep, drop=FALSE]
              obj@NumberxCol <- obj@NumberxCol[keep, ,drop=FALSE]
              obj@Number <- sum(keep)
              
            }
            
            return(obj)            
            
})


#' @rdname rm_monoclusters
#' @param subset_sce logical. Subset the SCE object to the cells which have co-clusters genes if is TRUE.
#' @export
setMethod(f = "rm_monoclusters",
          signature=(obj="SingleCellExperiment"),
          function(obj,
                   subset_sce = FALSE){
            
            if (isFALSE('caclust' %in% names(metadata(obj)))){
              stop('caclust not found from input SingleCellExperiment object metadata, run caclust function first!')
            }
            
            clust = metadata(obj)[['caclust']]
            clust = rm_monoclusters(clust)
            metadata(obj)[['caclust']] = clust
            
            if (isTRUE(subset_sce)){
              idxg = rownames(obj) %in% names(gene_clusters(clust))
              idxc = colnames(obj) %in% names(cell_clusters(clust))
              
              obj = obj[idxg, idxc]
            }
            
            return(obj)      
            
          })
