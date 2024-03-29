% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_overlap}
\alias{calc_overlap}
\title{WIP replacement for `determine_overlap` function}
\usage{
calc_overlap(cc_adj, cg_adj, threshold)
}
\arguments{
\item{cc_adj}{sparse matrix (dgCMatrix), cell-cell adjacency matrix}

\item{cg_adj}{sparse matrix (dgCMatrix), cell-gene adjacency matrix}

\item{threshold}{numeric value between 0 and 1. The cutoff of cell-nqighour-overlapping of each gene.}
}
\description{
c++ implementation for calculating the cell-neighour-overlapping among k nearest cell neighbours of each gene.
}
