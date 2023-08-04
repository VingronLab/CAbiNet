

#' Combine kNN graphs to large cell-gene adjecency matrix
#' @description
#' Builds a single adjacency matrix consisting of cells and genes from 4
#' seperate sub kNN-graphs.
#' @md
#' @param caobj A cacomp object with standard and principal coordinates
#' calculated.
#' @param k_c k for cell-cell kNN, integer.
#' @param k_g k for gene-gene kNN, integer.
#' @param k_cg k for cell-gene kNN, integer.
#' @param k_gc k for gene-cell kNN, interger.
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
#' @param method [BiocNeighbors::BiocNeighborParam] object specifying the
#'  algorithm to use. see Details.
#' @param BPPARAM [BiocParallel] settings parameter. By default single core 
#' [BiocParallel::SerialParam()] but other parameters can be passed.
#' 
#' @details 
#' `method` should be a kNN algorithm defined by 
#' [BiocNeighbors::BiocNeighborParam]. For exact kNN search use 
#' `BiocNeighbors::KmknnParam()` or `BiocNeighbors::VptreeParam()`.
#'
#' @return
#' Adjacency matrix of type `dgCMatrix`.
#' The combined adjacency matrix consists of the cell-cell graph, gene-gene
#' graph and cell-gene/gene-cell graph.
#'
#' @export
create_bigraph <- function(caobj,
                         k_c,
                         k_g,
                         k_cg,
                         k_gc,
                         loops = FALSE,
                         select_genes = TRUE,
                         prune_overlap = TRUE,
                         overlap = 0.2,
                         calc_gene_cell_kNN = FALSE,
                         marker_genes = NULL,
                         method = BiocNeighbors::KmknnParam(),
                         BPPARAM = BiocParallel::SerialParam()){
    
    # apply vector augmentation for MIP search via euclidean distance.
    Xt <- add_zero_dim(caobj@std_coords_cols)
    Qt <- augment_vector(caobj@prin_coords_rows)
    
    cgg_nn <- BiocNeighbors::queryKNN(X = Qt,
                       query = Xt,
                       k = k_cg,
                       get.distance = FALSE,
                       BNPARAM = method,
                       BPPARAM = BPPARAM)$index
    
    cgg_nn <- indx_to_spmat(indx_mat = cgg_nn,
                            row_names = rownames(caobj@std_coords_cols),
                            col_names = rownames(caobj@prin_coords_rows))
    

    # gene_idx <- seq_len(nrow(caobj@prin_coords_rows))
    
    if (!is.null(marker_genes)){
        stopifnot(is(marker_genes, "character"))
        
        idx <- which(colnames(cgg_nn) %in% marker_genes)
        
        if (length(idx) == 0){
            warning("Marker genes not found in the data.")
            marker_genes <- NULL
            
        } else {
            
            if(length(idx) < length(marker_genes)){
                warning("Not all marker genes are in the provided data.")
                marker_genes <- marker_genes[marker_genes %in% colnames(cgg_nn)]
            }
            
            marker_knn <- cgg_nn[,idx, drop=FALSE]
            cgg_nn <- cgg_nn[,-idx, drop = FALSE]
            

        }
        
    }
    

    ccg_nn = BiocNeighbors::findKNN(caobj@prin_coords_cols,
                     k=k_c,
                     get.distance = FALSE,
                     BNPARAM=method,
                     BPPARAM = BPPARAM)$index
    
    if (isTRUE(loops)){
        ccg_nn <- cbind(seq_len(nrow(ccg_nn)),
                        ccg_nn[, -ncol(ccg_nn), drop = FALSE])
    }
    
    ccg_nn <- indx_to_spmat(indx_mat = ccg_nn, 
                            row_names = rownames(caobj@prin_coords_cols), 
                            col_names = rownames(caobj@prin_coords_cols))

    
    if (isTRUE(select_genes)){
      
      idx <- Matrix::colSums(cgg_nn) > 0
      cgg_nn <- cgg_nn[,idx, drop = FALSE]

      if (isTRUE(prune_overlap)){
        
        overlap_mat <- calc_overlap( cc_adj = ccg_nn,
                                     cg_adj = cgg_nn)
        
        # For the case overlap = 1, all the genes are supposed to removed such that
        # the algorithm allows for clustering for cells without genes.
        cgg_nn[overlap_mat <= overlap] <- 0
        idx <- Matrix::colSums(cgg_nn) > 0
        cgg_nn <- cgg_nn[,idx, drop = FALSE]

      }
    } 
    
    

    if(!is.null(marker_genes)){
        
        cgg_nn <- cbind(cgg_nn, marker_knn)
        
        # marker_dists <- marker_dists[,c(colnames(gene_dists), rownames(marker_dists))]
        # gene_dists <- cbind(rbind(gene_dists, marker_dists[,colnames(gene_dists)]), t(marker_dists))
        # gene_cell_assr <- rbind(gene_cell_assr, marker_assr)
        
    }
    
    # gene_idx <- which(rownames(caobj@prin_coords_rows) %in% colnames(cgg_nn))
    
    gene_idx <- match(colnames(cgg_nn),
                      rownames(caobj@prin_coords_rows),
                      nomatch = NA_integer_)
    stopifnot(!any(is.na(gene_idx)))
    
    ggg_nn = BiocNeighbors::findKNN(caobj@prin_coords_rows[gene_idx,],
                     k=k_g,
                     get.distance = FALSE,
                     BNPARAM=method,
                     BPPARAM = BPPARAM)$index  
    
    if (isTRUE(loops)){
        ggg_nn <- cbind(seq_len(nrow(ggg_nn)),
                        ggg_nn[, -ncol(ggg_nn), drop = FALSE])
    }
    
    ggg_nn <- indx_to_spmat(indx_mat = ggg_nn, 
                            row_names = rownames(caobj@prin_coords_rows)[gene_idx], 
                            col_names = rownames(caobj@prin_coords_rows)[gene_idx])
    
    
    # ggg_nn <- make_knn(gene_dists,
    #                    k = k_g,
    #                    decr = FALSE,
    #                    loops = loops)
    
    if(isFALSE(calc_gene_cell_kNN)){
        gcg_nn <- Matrix::t(cgg_nn)
        
    } else if(isTRUE(calc_gene_cell_kNN)){
        
        Xt <- augment_vector(caobj@prin_coords_cols)
        Qt <- add_zero_dim(caobj@std_coords_rows[gene_idx,])
        
        gcg_nn <- BiocNeighbors::queryKNN(X = Xt,
                           query = Qt,
                           k = k_gc,
                           get.distance = FALSE,
                           BNPARAM=method,
                           BPPARAM = BPPARAM)$index
        
        gcg_nn <- indx_to_spmat(indx_mat = gcg_nn, 
                                row_names = rownames(caobj@std_coords_rows)[gene_idx], 
                                col_names = rownames(caobj@prin_coords_cols))
        
        # gcg_nn <- make_knn(gene_cell_assr,
        #                    k = k_gc,
        #                    decr = TRUE,
        #                    loops = TRUE)
    } else {
        stop("calc_cell_gene_kNN has to be either TRUE or FALSE!")
    }
    
    
    
    GSG_1 <- cbind(ccg_nn, cgg_nn)
    GSG_2 <- cbind(gcg_nn, ggg_nn)
    
    GSG <- rbind(GSG_1, GSG_2)
    return(GSG)
    
}


# 
# # TODO: ONLY MOVE ON IF YOU ARE SURE ABOUT THE INDICES.
# create_bigraph_biocneighbors_indxmat <- function(caobj,
#                                                    k_c,
#                                                    k_g,
#                                                    k_cg,
#                                                    k_gc,
#                                                    loops = FALSE,
#                                                    select_genes = TRUE,
#                                                    prune_overlap = TRUE,
#                                                    overlap = 0.2,
#                                                    calc_gene_cell_kNN = FALSE,
#                                                    marker_genes = NULL,
#                                                    method = BiocNeighbors::KmknnParam(),
#                                                    BPPARAM = BiocParallel::SerialParam()){
#     
#     # apply vector augmentation for MIP search via euclidean distance.
#     Xt <- add_zero_dim(caobj@std_coords_cols)
#     Qt <- augment_vector(caobj@prin_coords_rows)
#     
#     cgg_nn <- BiocNeighbors::queryKNN(X = Qt,
#                                       query = Xt,
#                                       k = k_cg,
#                                       get.distance = FALSE,
#                                       BNPARAM=method,
#                                       BPPARAM = BPPARAM)$index
#     
#     # (we cannot guarantee during pruning etc that the number of neighbors is the same)
#     # cgg_nn <- as.list(data.frame(t(cgg_nn)))     # convert to list of indices 
# 
#     rownames(cgg_nn) <- rownames(caobj@std_coords_cols)
# 
#     
#     if (!is.null(marker_genes)){
#         stopifnot(is(marker_genes, "character"))
#         
#         idx <- which(rownames(caobj@prin_coords_rows) %in% marker_genes)
#         
#         if (length(idx) == 0){
#             warning("Marker genes not found in the data.")
#             marker_genes <- NULL
#             
#         } else {
#             
#             if(length(idx) < length(marker_genes)){
#                 warning("Not all marker genes are in the provided data.")
#                 marker_genes <- marker_genes[marker_genes %in% rownames(caobj@prin_coords_rows)]
#             }
#             
#             marker_knn <- base::apply(cgg_nn, 1, function(x) idx[idx %in% x])
#             
#             # We keep it in the adjacency matrix bc calc_overlap expects it.
#             # setdiff removes the marker genes from the graph, we add it later.
#             # cgg_nn <- lapply(seq_len(length(cgg_nn)), function(x) setdiff(cgg_nn[[x]], marker_knn[[x]]))      
#             # names(cgg_nn) <- rownames(caobj@std_coords_cols)
#             
#         }
#         
#     }
#     
#     ccg_nn = BiocNeighbors::findKNN(caobj@prin_coords_cols,
#                                     k=k_c,
#                                     get.distance = FALSE,
#                                     BNPARAM=method,
#                                     BPPARAM = BPPARAM)$index
#     
#     rownames(ccg_nn) <- rownames(caobj@prin_coords_cols)
#     
#     if (isTRUE(loops)){
#       ccg_nn <- cbind(seq_len(nrow(ccg_nn)),
#                       ccg_nn[, -ncol(ccg_nn), drop = FALSE])
#     }
#     
#     # ccg_nn <- as.list(data.frame(t(ccg_nn)))
# 
#     if (isTRUE(select_genes) & isTRUE(prune_overlap)){
# 
#       # FIXME: Change calc_overlap for index matrices.
#       # TODO: Ensure that calc_overlap removes genes.
#       overlap_mat <- calc_overlap( cc_adj = ccg_nn,
#                                    cg_adj = cgg_nn)
# 
#     }
#     
#     # add marker genes back in BEFORE we get gene_idx!
#     if(!is.null(marker_genes)){
# 
#         cgg_nn <- do.call("rbind", 
#                    lapply(seq_len(nrow(cgg_nn)),
#                           function(x) add2knn(cgg_nn[x,], marker_knn[[x]]))
#                    )
#     
#         rownames(cgg_nn) <- rownames(caobj@std_coords_cols)
#       
#     }
#     
#     if (isTRUE(select_genes)){
#       # indices of genes with an edge to a cell.
#       # If we subset to the genes that have an edge to a cell.
#       gene_idx <- sort(unique(as.numeric(cgg_nn)))
#       
#     } else {
#       # If we do not subset.
#       gene_idx <- seq_len(nrow(caobj@prin_coords_rows))
#     }
#     
#     
#     ggg_nn = BiocNeighbors::findKNN(caobj@prin_coords_rows[gene_idx,],
#                                     k=k_g,
#                                     get.distance = FALSE,
#                                     BNPARAM=method,
#                                     BPPARAM = BPPARAM)$index  
#     
#     rownames(ggg_nn) <- rownames(caobj@prin_coords_rows[gene_idx,])
#     
#     if (isTRUE(loops)){
#       ggg_nn <- cbind(seq_len(nrow(ggg_nn)),
#                       ggg_nn[, -ncol(ggg_nn), drop = FALSE])
#     }
#     
#     # ggg_nn <- as.list(data.frame(t(ggg_nn)))
#     
#     
#     if(isFALSE(calc_gene_cell_kNN)){
#       
#       # for each gene that has an edge to a cell (not all genes!!):
#       # check to which cells it has an edge and put their indices in list.
#       # TODO: Check if correct.
#       gcg_nn <- lapply(gene_idx,
#                        function(y) which(vapply(cgg_nn,
#                                                 function(x) y %in% x,
#                                                 TRUE)))
#       max_cells <- max(lengths(gcg_nn))
#       
#       gcg_nn <- lapply(gcg_nn,
#                        function(x) c(x, rep(NA_integer_, max_cells-length(x))))
#       
#       gcg_nn <- do.call("rbind", gcg_nn)
#       
#       rownames(gcg_nn) <- rownames(caobj@prin_coords_rows)[gene_idx]
#       
#       
#     } else if(isTRUE(calc_gene_cell_kNN)){
#       
#       Xt <- augment_vector(caobj@prin_coords_cols)
#       Qt <- add_zero_dim(caobj@std_coords_rows[gene_idx,])
#       
#       # calculated cell (!) indices for the subsetted genes.
#       gcg_nn <- BiocNeighbors::queryKNN(X = Xt,
#                                         query = Qt,
#                                         k = k_gc,
#                                         get.distance = FALSE,
#                                         BNPARAM=method,
#                                         BPPARAM = BPPARAM)$index
#       
#       rownames(gcg_nn) <- rownames(caobj@std_coords_rows[gene_idx,])
#       
#       # gcg_nn <- as.list(data.frame(t(gcg_nn)))
#       
#         
#     } else {
#       stop("calc_cell_gene_kNN has to be either TRUE or FALSE!")
#     }
#     
#     # Reindex cell-gene graph so that the gene indixes refer only to the gene
#     # that are actually left after the pruning. This ensures continuous indxs 
#     # and that the adj matrix GSG is self contained (only references itself).
#     
#     # cgg_nn <- lapply(cgg_nn, function(x) match(x, gene_idx, nomatch = NA))
#     
#     cgg_nn <- do.call("rbind", 
#                         apply(cgg_nn,
#                               1,
#                               function(x) match(x,
#                                                 gene_idx,
#                                                 nomatch = NA),
#                               simplify = FALSE)
#                 )
#     
#     
#     # stopifnot(!anyNA(cgg_nn))
#     
#     ncells = nrow(ccg_nn)
# 
#     stopifnot(nrow(ccg_nn) == nrow(cgg_nn))
#     
#     GSG_1 <- cbind(ccg_nn, cgg_nn + ncells)
#     
# 
#     rownames(GSG_1) <- rownames(ccg_nn)
#     
#     stopifnot(nrow(gcg_nn) == nrow(ggg_nn))
#     # stopifnot(nrow(gcg_nn) == nrow(ggg_nn))
#     
#     GSG_2 <- cbind(gcg_nn, ggg_nn + ncells)
# 
#     rownames(GSG_2) <- rownames(ggg_nn)
#     
#     stopifnot(ncol(GSG_1) == ncol(GSG_2))
#     
#     GSG <- rbind(GSG_1, GSG_2)
# 
# 
#     return(GSG)
#     
# }




#' Create SNN-graph from caobj
#'
#' @family biclustering
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

  stopifnot(mode %in% c("out", "in", "all"))

  adj <- create_bigraph(caobj = caobj,
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
