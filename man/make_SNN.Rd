% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.R
\name{make_SNN}
\alias{make_SNN}
\title{Create SNN-graph from caobj}
\usage{
make_SNN(
  caobj,
  k,
  SNN_prune = 1/15,
  loops = FALSE,
  mode = "out",
  select_genes = TRUE,
  prune_overlap = TRUE,
  overlap = 0.2,
  calc_gene_cell_kNN = FALSE,
  marker_genes = NULL,
  method = BiocNeighbors::KmknnParam(),
  BPPARAM = BiocParallel::SerialParam()
)
}
\arguments{
\item{caobj}{A cacomp object with standard and principal coordinates
calculated.}

\item{k}{Either an integer (same k for all subgraphs) or a vector of
exactly four integers specifying in this order:
\itemize{
\item k_c for the cell-cell kNN-graph
\item k_g for the gene-gene kNN-graph
\item k_cg for the cell-gene kNN-graph
\item k_gc for the gene-cell kNN-graph.
}}

\item{SNN_prune}{numeric. Value between 0-1. Sets cutoff of acceptable jaccard
similarity scores for neighborhood overlap of vertices in SNN. Edges with values
less than this will be set as 0. The default value is 1/15.}

\item{loops}{TRUE/FALSE. If TRUE self-loops are allowed, otherwise not.}

\item{mode}{The type of neighboring vertices to use for calculating similarity
scores(Jaccard Index). Three options: "out", "in" and "all":
\itemize{
\item "out": Selecting neighbouring vertices by out-going edges;
\item "in": Selecting neighbouring vertices by in-coming edges;
\item "all": Selecting neigbouring vertices by both in-coming and out-going edges.
}}

\item{select_genes}{TRUE/FALSE. Should genes be selected by whether they have
an edge in the cell-gene kNN graph?}

\item{prune_overlap}{TRUE/FALSE. If TRUE edges to genes that share less
than \code{overlap} of genes with the nearest neighbours of the cell are removed.
Pruning is only performed if select_genes = TRUE.}

\item{overlap}{Numeric between 0 and 1. Overlap cutoff applied if
prune_overlap = TRUE.}

\item{calc_gene_cell_kNN}{TRUE/FALSE. If TRUE a cell-gene graph is calculated
by choosing the \code{k_gc} nearest cells for each gene. If FALSE the cell-gene
graph is transposed.}

\item{marker_genes}{character. Optional. Names of known marker genes that
should be excempt from any pruning on the graph and be kept.}

\item{method}{\link[BiocNeighbors:BiocNeighborParam]{BiocNeighbors::BiocNeighborParam} object specifying the
algorithm to use. see Details.}

\item{BPPARAM}{\link{BiocParallel} settings parameter. By default single core
\code{\link[BiocParallel:SerialParam-class]{BiocParallel::SerialParam()}} but other parameters can be passed.}
}
\value{
A sparse adjacency Matrix of type "dgCMatrix". The values in the matrix
are the Jaccard similarity between nodes in the graph. The range between 0
and 1, with 0 meaning that no edges are shared between nodes, wheras 1 means
all edges are shared between nodes.
}
\description{
Builds a shared nearest neighbour graph (SNN) from a "cacomp" object.
}
\seealso{
Other biclustering: 
\code{\link{caclust}},
\code{\link{run_leiden}()},
\code{\link{run_spectral}()}
}
\concept{biclustering}
