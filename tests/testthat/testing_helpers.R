
#Rcpp::sourceCpp("../../src/testing_funs.cpp")


#' Make a kNN graph
#'
#' @description Given a distance matrix, `make_knn()` builds up a  k-nearest-neighbours graph
#' and outputs the adjacency matrix.
#'
#' @param dists distance matrix
#' @param k integer. Number of nearest neighbours.
#' @param decr boolean. Whether the the values in `dists` should be sorted into
#' decreasing (TRUE) or increasing (FALSE) order.
#' @param loops TRUE/FALSE. If TRUE self-loops are allowed, otherwise not.
#'
#' @return
#' adjacency matrix of kNN graph: sparse matrix of type "dgCMatrix".
#'
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
    
    if (n.col == k & isFALSE(loops)) {
        warning(
            "k equal to ncol and loops = FALSE. Setting k to ncol - 1.",
            call. = FALSE
        )
        k <- n.col - 1
    }
    

    knn.mat <- matrix(data = 0, ncol = k, nrow = n.row)
    # knd.mat <- knn.mat
    
    if (isTRUE(loops)){
        nns <- seq_len(k)
        
    } else if (isFALSE(loops)) {
        if (n.col ==  k){
            nns <- seq_len(k)[-1]
        } else{
            nns <- seq_len(k)+1
        }
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
#' Depracted! Please use create_bigraph instead.
#' 
#' @description
#' Builds a single adjacency matrix consisting of cells and genes from 4
#' seperate sub kNN-graphs.
#'
#' @param cell_dists cell-cell euclidean distances matrix
#' @param gene_dists gene-gene euclidean distances matrix
#' @param cell_gene_assr cell-gene association ratios matrix
#' @param gene_cell_assr gene-cell association ratios matrix
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
#'
#' @return
#' Adjacency matrix of type `dgCMatrix`.
#' The combined adjacency matrix consists of the cell-cell graph, gene-gene
#' graph and cell-gene/gene-cell graph.
#'
create_bigraph_manual <- function(cell_dists,
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
                       loops = TRUE)
    
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
            
            marker_knn <- cgg_nn[,idx, drop = FALSE]
            cgg_nn <- cgg_nn[,-idx, drop = FALSE]
            
            marker_dists <- gene_dists[idx,,drop = FALSE]
            marker_assr <- gene_cell_assr[idx,,drop = FALSE]
            
            gene_dists <- gene_dists[-idx, -idx, drop = FALSE]
            gene_cell_assr <- gene_cell_assr[-idx,, drop = FALSE]
            
        }
        
    }
    
    if(isTRUE(select_genes)){
        
        idx <- Matrix::colSums(cgg_nn) > 0
        cgg_nn <- cgg_nn[,idx, drop = FALSE]
        gene_dists <- gene_dists[idx,idx, drop = FALSE]
        gene_cell_assr <- gene_cell_assr[idx,,drop = FALSE]
    }
    
    
    ccg_nn <- make_knn(cell_dists,
                       k = k_c,
                       decr = FALSE,
                       loops = loops)
    
    
    if(isTRUE(select_genes) & isTRUE(prune_overlap)){

        overlap_mat <- calc_overlap_deprecated( cc_adj = ccg_nn,
                                               cg_adj = cgg_nn)
        print('original overlap_mat:')
        print(overlap_mat)
        # For the case overlap = 1, all the genes are supposed to removed such that
        # the algorithm allows for clustering for cells without genes.
        cgg_nn[overlap_mat <= overlap] <- 0
        idx <- Matrix::colSums(cgg_nn) > 0
        cgg_nn <- cgg_nn[,idx, drop = FALSE]
        gene_dists <- gene_dists[idx,idx, drop = FALSE]
        gene_cell_assr <- gene_cell_assr[idx,,drop = FALSE]
        
    }
    
    
    if(!is.null(marker_genes)){
        
        cgg_nn <- cbind(cgg_nn, marker_knn)
        
        marker_dists <- marker_dists[,c(colnames(gene_dists), rownames(marker_dists)), drop = FALSE]
        gene_dists <- cbind(rbind(gene_dists, marker_dists[,colnames(gene_dists), drop = FALSE]), t(marker_dists))
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
                           loops = TRUE)
    } else {
        stop("calc_cell_gene_kNN has to be either TRUE or FALSE!")
    }
    
    
    
    GSG_1 <- cbind(ccg_nn, cgg_nn)
    GSG_2 <- cbind(gcg_nn, ggg_nn)
    
    GSG <- rbind(GSG_1, GSG_2)
    return(GSG)
}




