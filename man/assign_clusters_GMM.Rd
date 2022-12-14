% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{assign_clusters_GMM}
\alias{assign_clusters_GMM}
\title{Assign cluster to cells/genes}
\usage{
assign_clusters_GMM(caclust_obj, type = "genes", cutoff = 0.5)
}
\arguments{
\item{caclust_obj}{A caclust object}

\item{type}{Either "cell" or "gene".}

\item{cutoff}{Probability cutoff.}
}
\value{
logical matrix indicating which gene belongs to which cluster.
}
\description{
Based on a probability cutoff genes/cells are assigned to all clusters for
which they have a probability higher than 'cutoff'.
}
