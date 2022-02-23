
#' Checks if caclust-class was constructed correctly
check_caclust <- function(object){
  
  errors <- character()
  
  if (!is(object, "caclust")){
    
    msg <- "'object' is not of class caclust!"
    errors <- c(errors, msg)
    
  }
  
  #TODO
  
  if (length(errors) == 0) TRUE else errors
}

#' An S4 class that contains all elements needed for CA.
#' @name caclust-class
#' @rdname caclust-class
#' @description
#' TODO
#'
#' @slot cell_clusters
#' @slot gene_clusters
#' @slot SNN
#' @export
setClass("caclust",
         representation(
           cell_clusters = "factor",
           gene_clusters = "factor",
           SNN = "dgCMatrix"
         ),
         validity = check_caclust
)



#' Create new caclust object
new_caclust <- function(...) new("caclust",...)

cell_clusters <- function(object){
  slot(object, "cell_clusters")
}

gene_clusters <- function(object){
  slot(object, "gene_clusters")
}

get_snn <- function(object){
  slot(object, "SNN")
}

show.caclust <- function(object){
  
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
setMethod(f = "show", signature(object = "caclust"), function(object) {
  show.caclust(object)
})