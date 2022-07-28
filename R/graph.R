



#' Make a kNN graph
#' 
#' @description Given a distance matrix `make_knn()` outputs the adjacency
#' matrix of the k-nearest-neighbours graph.
#'
#' @param dists distance matrix
#' @param k numeric. Number of nearest neighbours.
#' @param decr boolean. Whether the the values in `dists` should be sorted into
#' decreasing (TRUE) or increasing (FALSE) order.
#' @param loops TRUE/FALSE. If TRUE self-loops are allowed, otherwise not.
#'
#' @return
#' adjacency matrix of kNN graph: sparse matrix of type "dgCMatrix".
#'
#' @export
make_knn <- function(dists,
                     k,
                     decr = TRUE,
                     loops = FALSE) {
  
  n.row <- nrow(dists)
  n.col <- ncol(dists)
  
  if (n.col < k) {
    warning(
      "k set larger than number of genes Setting k to number of genes - 1.",
      call. = FALSE
    )
    k <- n.col - 1
  }
  knn.mat <- matrix(data = 0, ncol = k, nrow = n.row)
  # knd.mat <- knn.mat
  
  if (isTRUE(loops)){
    nns <- seq_len(k)
    
  } else if (isFALSE(loops)) {
    nns <- seq_len(k)+1
  } else {
    stop()
  }
  
  for (i in 1:n.row) {
    knn.mat[i, ] <- order(dists[i, ], decreasing = decr)[nns]
    # knd.mat[i, ] <- dists[i, knn.mat[i, ]]
  }
  # nn.ranked <- knn.mat[,seq_len(k)]
  
  # convert nn.ranked into a Graph
  j <- as.numeric(t(knn.mat))
  i <- ((1:length(j)) - 1) %/% k + 1
  
  nn.matrix <- Matrix::sparseMatrix(i = i,
                                    j = j,
                                    x = 1,
                                    dims = c(nrow(x = dists), ncol(x = dists)))
  rownames(nn.matrix) <- rownames(dists)
  colnames(nn.matrix) <- colnames(dists)
  return(nn.matrix)
}






#' Combine kNN graphs to large cell-gene adjecency matrix
#'
#' @description
#' Builds a single adjacency matrix consisting of cells and genes from 4 
#' seperate sub kNN-graphs.
#' 
#' @param cell_dists cell-cell euclidean distances
#' @param gene_dists gene-gene euclidean distances
#' @param cell_gene_assr cell-gene association ratios
#' @param gene_cell_assr gene-cell association ratios
#' @param k_c k for cell-cell kNN
#' @param k_g k for gene-gene kNN
#' @param k_cg k for cell-gene kNN
#' @param k_gc k for gene-cell kNN
#' @param loops TRUE/FALSE. If TRUE self-loops are allowed, otherwise not.
#' @param select_genes TRUE/FALSE. Should genes be selected by whether they have
#' an edge in the cell-gene kNN graph?
#' @param prune_overlap TRUE/FALSE. If TRUE edges to genes that share less
#' than `overlap` of genes with the nearest neighbours of the cell are removed.
#' Pruning is only performed if select_genes = TRUE.
#' @param overlap Numeric between 0 and 1. Overlap cutoff applied if
#' prune_overlap = TRUE.
#' @param calc_gene_cell_kNN TRUE/FALSE. If TRUE a cell-gene graph is calculated
#' by choosing the `k_gc` nearest cells for each gene. If FALSE the cell-gene
#' graph is transposed.
#' @param marker_genes character. Optional. Names of known marker genes that 
#' should be excempt from any pruning on the graph and be kept.
#'
#' @return
#' Adjacency matrix of type `dgCMatrix`. 
#' The combined adjacency matrix consists of the cell-cell graph, gene-gene 
#' graph and cell-gene/gene-cell graph.
#'
#' @export
create_bigraph <- function(cell_dists,
                           gene_dists,
                           cell_gene_assr,
                           gene_cell_assr,
                           k_c,
                           k_g,
                           k_cg,
                           k_gc,
                           loops = FALSE,
                           select_genes = TRUE,
                           prune_overlap = TRUE,
                           overlap = 0.2,
                           calc_gene_cell_kNN = FALSE,
                           marker_genes = NULL) {
  
  cgg_nn <- make_knn(cell_gene_assr,
                     k = k_cg,
                     decr = TRUE,
                     loops = loops)
  
  if (!is.null(marker_genes)){
    stopifnot(is(marker_genes, "character"))
    
    idx <- which(colnames(cgg_nn) %in% marker_genes)
    marker_knn <- cgg_nn[,idx]
    cgg_nn <- cgg_nn[,-idx]
    
    marker_dists <- gene_dists[idx,]
    marker_assr <- gene_cell_assr[idx,]
    
    gene_dists <- gene_dists[-idx, -idx]
    gene_cell_assr <- gene_cell_assr[-idx,]
  }
  
  if(isTRUE(select_genes)){
    
    idx <- Matrix::colSums(cgg_nn) > 0
    cgg_nn <- cgg_nn[,idx]
    gene_dists <- gene_dists[idx,idx]
    gene_cell_assr <- gene_cell_assr[idx,]
  }
  
  
  ccg_nn <- make_knn(cell_dists,
                     k = k_c,
                     decr = FALSE,
                     loops = loops)
  
  
  if(isTRUE(select_genes) & isTRUE(prune_overlap)){
    
    overlap_mat <- determine_overlap(cg_adj = cgg_nn,
                                     cc_adj = ccg_nn)
    
    cgg_nn[overlap_mat < overlap] <- 0
    
    idx <- Matrix::colSums(cgg_nn) > 0
    cgg_nn <- cgg_nn[,idx]
    gene_dists <- gene_dists[idx,idx]
    gene_cell_assr <- gene_cell_assr[idx,]
    
  }
  
  
  if(!is.null(marker_genes)){
    
    cgg_nn <- cbind(cgg_nn, marker_knn)
    
    marker_dists <- marker_dists[,c(colnames(gene_dists), rownames(marker_dists))]
    gene_dists <- cbind(rbind(gene_dists, marker_dists[,colnames(gene_dists)]), t(marker_dists))
    gene_cell_assr <- rbind(gene_cell_assr, marker_assr)
    
  }
  
  ggg_nn <- make_knn(gene_dists,
                     k = k_g,
                     decr = FALSE,
                     loops = loops)
  
  if(isFALSE(calc_gene_cell_kNN)){
    gcg_nn <- Matrix::t(cgg_nn)
    
  } else if(isTRUE(calc_gene_cell_kNN)){
    gcg_nn <- make_knn(gene_cell_assr,
                       k = k_gc,
                       decr = TRUE,
                       loops = loops)
  } else {
    stop("calc_cell_gene_kNN has to be either TRUE or FALSE!")
  }
  
  
  
  GSG_1 <- cbind(ccg_nn, cgg_nn)
  GSG_2 <- cbind(gcg_nn, ggg_nn)
  
  GSG <- rbind(GSG_1, GSG_2)
  return(GSG)
}





#' Create SNN-graph from caobj
#'
#' @description 
#' Builds a shared nearest neighbour graph (SNN) from a "cacomp" object.
#' 
#' @param caobj A cacomp object with standard and principal coordinates 
#' calculated.
#' @param k Either an integer (same k for all subgraphs) or a vector of 
#' exactly four integers specifying in this order: 
#' * k_c for the cell-cell kNN-graph
#' * k_g for the gene-gene kNN-graph
#' * k_cg for the cell-gene kNN-graph
#' * k_gc for the gene-cell kNN-graph.
#' @param SNN_prune numeric. Value between 0-1. Sets cutoff of acceptable jaccard 
#' similarity scores for neighborhood overlap of vertices in SNN. Edges with values 
#' less than this will be set as 0. The default value is 1/15.
#' @param mode The type of neighboring vertices to use for calculating similarity
#'  scores(Jaccard Index). Three options: "out", "in" and "all":
#' * "out": Selecting neighbouring vertices by out-going edges;
#' * "in": Selecting neighbouring vertices by in-coming edges;
#' * "all": Selecting neigbouring vertices by both in-coming and out-going edges.
#' @inheritParams create_bigraph
#' 
#' @returns 
#' A sparse adjacency Matrix of type "dgCMatrix". The values in the matrix
#' are the Jaccard similarity between nodes in the graph. The range between 0
#' and 1, with 0 meaning that no edges are shared between nodes, wheras 1 means 
#' all edges are shared between nodes.
#' 
#' @md 
#' @export
make_SNN <- function(caobj,
                       k,
                       SNN_prune = 1/15,
                       loops = FALSE,
                       mode = "out",
                       select_genes = TRUE,
                       prune_overlap = TRUE,
                       overlap = 0.2,
                       calc_gene_cell_kNN = FALSE,
                       marker_genes = NULL) {
  
  if (length(k) == 1){
    k_c <- k_g <- k_cg <- k_gc <- k
  } else if (length(k) == 4){
    k_c <- k[1]
    k_g <- k[2]
    k_cg <- k[3]
    k_gc <- k[4]
  } else {
    stop("invalid k.")
  }
  
  distances <- calc_distances(caobj = caobj)
  
  stopifnot(all(c("cc", "gg", "cg", "gc") %in% names(distances)))
  stopifnot(mode %in% c("out", "in", "all"))
  
  adj <- create_bigraph(cell_dists = distances[["cc"]],
                        gene_dists = distances[["gg"]],
                        cell_gene_assr = distances[["cg"]],
                        gene_cell_assr = distances[["gc"]],
                        k_c = k_c,
                        k_g = k_g,
                        k_cg = k_cg,
                        k_gc = k_gc,
                        loops = loops,
                        overlap = overlap,
                        prune_overlap = prune_overlap,
                        select_genes = select_genes,
                        calc_gene_cell_kNN = calc_gene_cell_kNN,
                        marker_genes = marker_genes)
  
  if(!is(adj, "dgCMatrix")){
    adj <- as(adj, "dgCMatrix")  
  }
  
  
  snn.matrix <- ComputeSNNasym(adj, SNN_prune, mode = mode)
  
  ## to coincide with output of "igraph"
  diag(snn.matrix) = 1
  
  rownames(snn.matrix) <- rownames(adj)
  colnames(snn.matrix) <- rownames(adj)
  
  cidxs <- which(rownames(snn.matrix) %in% rownames(caobj@std_coords_cols))
  gidxs <- which(rownames(snn.matrix) %in% rownames(caobj@std_coords_rows))
  
  caclust <- new("caclust",
                 SNN=snn.matrix,
                 cell_idxs = cidxs,
                 gene_idxs = gidxs)
  return(caclust)
}