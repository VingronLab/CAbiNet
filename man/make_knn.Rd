% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.R
\name{make_knn}
\alias{make_knn}
\title{Make a kNN graph}
\usage{
make_knn(dists, k, decr = TRUE, loops = FALSE)
}
\arguments{
\item{dists}{distance matrix}

\item{k}{integer. Number of nearest neighbours.}

\item{decr}{boolean. Whether the the values in `dists` should be sorted into
decreasing (TRUE) or increasing (FALSE) order.}

\item{loops}{TRUE/FALSE. If TRUE self-loops are allowed, otherwise not.}
}
\value{
adjacency matrix of kNN graph: sparse matrix of type "dgCMatrix".
}
\description{
Given a distance matrix, `make_knn()` builds up a  k-nearest-neighbours graph
and outputs the adjacency matrix.
}
