
#' Helper function to check if object is empty.
#' @param x object
#' @return TRUE if x has length 0 and is not NULL. FALSE otherwise
is.empty <- function(x) return(isTRUE(length(x) == 0 & !is.null(x)))


