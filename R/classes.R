
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
  
  if(nrow(object@SNN) != ncol(object@SNN)){
    msg <- "SNN number of rows not equal to number of columns!"
    errors <- c(errors, msg)
  }
  
  if(is.null(names(object@cell_clusters))){
    msg <- "cell clusters are not named!"
    errors <- c(errors, msg)
  }
  
  if(is.null(names(object@gene_clusters))){
    msg <- "gene clusters are not named!"
    errors <- c(errors, msg)
  }
  
  if(!identical(rownames(object@SNN), colnames(object@SNN))){
    msg <- "SNN rownames not identical to colnames!"
    errors <- c(errors, msg)
  }
  
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
#' @slot cell_prob matrix.
#' @slot gene_prob matrix.
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
           gene_prob = "matrix"
         ),
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

#' Print caclust object in console
#' @param object a caclust object
show.caclust <- function(object, n_rows = 10){
  
  stopifnot(is(object, "caclust"))
  
  ncells <- length(cell_clusters(object)) 
  ngenes <- length(gene_clusters(object)) 
  cat("caclust object with", ncells, "cells and", ngenes, "genes.")

  stopifnot(identical(levels(cell_clusters(object)),
                      levels(gene_clusters(object))))
  cat("\n")
  cat(length(levels(gene_clusters(object))), "clusters found.")
  cat("\nClustering results:\n\n")
  df <- data.frame("cluster" = levels(cell_clusters(object)),
                   "ncells" = summary(cell_clusters(object), maxsum = Inf),
                   "ngenes" = summary(gene_clusters(object), maxsum = Inf))
  
  print(df, row.names = FALSE, right = FALSE)
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
#' @param bic biclustering results from caclust or biclust
#' 
#' @return 
#' caclust/biclust object with monoclusters removed.
#' 
#' @export
setGeneric("rm_monoclusters", function(bic) {
  standardGeneric("rm_monoclusters")
})


#' @rdname rm_monoclusters
#' @export
setMethod(f = "rm_monoclusters",
          signature=(bic="caclust"),
          function(bic){
        
        cc <- unique(cell_clusters(bic))
        gc <- unique(gene_clusters(bic))
        
        keep <- sort(as.character(intersect(cc,gc)))
        
        if(length(keep > 0)){
          
          bic@cell_clusters <- bic@cell_clusters[as.character(bic@cell_clusters) %in% keep]
          bic@cell_clusters <- droplevels(bic@cell_clusters)
          # bic@cell_clusters <- factor(as.character(bic@cell_clusters), levels = keep)
          # bic@cell_clusters <- na.omit(bic@cell_clusters)
          
          bic@gene_clusters <- bic@gene_clusters[as.character(bic@gene_clusters) %in% keep]
          bic@gene_clusters <- droplevels(bic@gene_clusters)
          # bic@gene_clusters <- factor(as.character(bic@gene_clusters), levels = keep)
          # bic@gene_clusters <- na.omit(bic@gene_clusters)
          
        }
        
        stopifnot(!is.null(names(bic@cell_clusters)))
        stopifnot(!is.null(names(bic@gene_clusters)))
        
        if(!is.empty(bic@SNN)){
          
          selr <- which(rownames(bic@SNN) %in% c(names(bic@cell_clusters),
                                                 names(bic@gene_clusters)))
          
          selc <- which(colnames(bic@SNN) %in% c(names(bic@cell_clusters),
                                                 names(bic@gene_clusters)))
          bic@SNN <- bic@SNN[selr, selc]
          
        }
        
        return(bic)      
                       
})


#' @rdname rm_monoclusters
#' @export
setMethod(f = "rm_monoclusters",
          signature=(bic="Biclust"),
          function(bic){
            
            keep <- colSums(bic@RowxNumber) > 0 & rowSums(bic@NumberxCol) > 0
            
            if(any(!keep)){
              
              bic@RowxNumber <- bic@RowxNumber[,keep, drop=FALSE]
              bic@NumberxCol <- bic@NumberxCol[keep, ,drop=FALSE]
              bic@Number <- sum(keep)
              
            }
            
            return(bic)            
            
})



