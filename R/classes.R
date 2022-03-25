
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
  
  if(!identical(levels(object@cell_clusters), levels(object@gene_clusters))){
    msg <- "factor levels for cells and genes are not the same!"
    errors <- c(errors, msg)
  }
  
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
#' @export
setClass("caclust",
         representation(
           cell_clusters = "factor",
           gene_clusters = "factor",
           SNN = "dgCMatrix",
           eigen = 'matrix',
           parameters = "list"
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

#' Print caclust object in console
#' @param object a caclust object
show.caclust <- function(object){
  
  stopifnot(is(object, "caclust"))
  
  ncells <- length(cell_clusters(object)) 
  ngenes <- length(gene_clusters(object)) 
  cat("caclust object with", ncells, "cells and", ngenes, "genes.")

  stopifnot(identical(levels(cell_clusters(object)),
                      levels(gene_clusters(object))))
  
  cat("\nClustering results:\n\n")
  df <- data.frame("cluster" = levels(cell_clusters(object)),
                   "ncells" = summary(cell_clusters(object)),
                   "ngenes" = summary(gene_clusters(object)))
  print(df, row.names = FALSE, right = FALSE)
  
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
  
  ctypes = unique(cc)
  gtypes = unique(gc)
  bitypes = intersect(ctypes, gtypes)
  
  Number = length(bitypes)
  if (Number == 0){
    NumberxCol = matrix(0)
    RowxNumber = matrix(0)
  }else{
    NumberxCol = do.call(rbind, lapply(bitypes, function(x){cc == x}))
    RowxNumber = do.call(cbind, lapply(bitypes, function(x){gc == x}))
  }
  bic <- new("Biclust", "Parameters" = params,
                       "RowxNumber" = RowxNumber,
                       "NumberxCol" = NumberxCol,
                       "Number" = Number,
                       "info" = list("Output of CAclust package"))
  
  return(bic)
}

