% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{determine_overlap}
\alias{determine_overlap}
\title{Determines the overlap of chosen genes between nearest neighbour cells.
It works faster on sparse matrix than on dense matrix.}
\usage{
determine_overlap(cg_adj, cc_adj)
}
\arguments{
\item{cg_adj}{cell-gene adjacency matrix}

\item{cc_adj}{cell-cell adjacency matrix}
}
\value{
Matrix (cells x genes) showing the overlap of each gene between the nearest
neighbour cells.
}
\description{
Determines how many genes are shared (have an edge) between
neighbouring cells for all cells in the dataset.
}
