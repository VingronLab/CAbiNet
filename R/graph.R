
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

    org_cellnames = rownames(caobj@std_coords_cols)
    org_genenames = rownames(caobj@prin_coords_rows)

    cgg_nn <- indx_to_spmat(indx_mat = cgg_nn,
                            row_names = org_cellnames,
                            col_names = org_genenames)


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
                            row_names = org_cellnames,
                            col_names = org_cellnames)


    if (isTRUE(select_genes)){

      idx <- Matrix::colSums(cgg_nn) > 0
      cgg_nn <- cgg_nn[,idx, drop = FALSE]
      new_genenames = colnames(cgg_nn)


      if (isTRUE(prune_overlap)){
        # This cpp function maps to the input spare matrices directly, and modify their values on site.
        # The weight of edges samller than overlap in cgg_nn will be set as 0 directly
 
        cgg_nn <- calc_overlap( cc_adj = ccg_nn,
                                cg_adj = cgg_nn,
                                threshold = overlap)

        # For the case overlap = 1, all the genes are supposed to removed such that
        # the algorithm allows for clustering for cells without genes.
        # cgg_nn[overlap_mat <= overlap] <- 0 # this step is done by cpp function to reduce the copy between cpp and R
        idx <- Matrix::colSums(cgg_nn) > 0
        cgg_nn <- cgg_nn[,idx, drop = FALSE]
        colnames(cgg_nn) <- new_genenames[idx]
        rownames(cgg_nn) <- org_cellnames

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
                      org_genenames,
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
                            row_names = org_genenames[gene_idx],
                            col_names = org_genenames[gene_idx])


    if(isFALSE(calc_gene_cell_kNN)){
        gcg_nn <- Matrix::t(cgg_nn)

        if (k_gc != k_cg){
          warning('The given values of k_gc and k_cg are different, But the calc_cell_gene_kNN is FALSE,
          then the gene-cell graph adjacency matrix will be calculated as the transpose of cell-gene graph adjacency matrix.
          This will ignore the given k_gc value. If you want to give k_gc and k_cg different values, set calc_cell_gene_kNN as TRUE.')
        }

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
                                row_names = org_genenames[gene_idx],
                                col_names = org_cellnames)

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
                       marker_genes = NULL,
                       method = BiocNeighbors::KmknnParam(),
                       BPPARAM = BiocParallel::SerialParam()) {

  if (length(k) == 1){
    k_c <- k_g <- k_cg <- k_gc <- k
  } else if (length(k) == 4){
    k_c <- k[1]
    k_g <- k[2]
    k_cg <- k[3]
    k_gc <- k[4]
  } else {
    stop("Invalid k. k should be either an interger or a vector with four integers. See ?make_SNN.")
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
                        marker_genes = marker_genes,
                        method = BiocNeighbors::KmknnParam(),
                        BPPARAM = BiocParallel::SerialParam())


  if(!is(adj, "dgCMatrix")){
    adj <- as(adj, "dgCMatrix")
  }


  snn.matrix <- ComputeSNNasym(adj, prune = SNN_prune, mode = mode)
  ## use memory mapping instead of copying

  # ## to coincide with output of "igraph"
  Matrix::diag(snn.matrix) = 1

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
