

#' Plot biUMAP
#' 
#' @param umap_coords data frame as outputted by `run_biUMAP_*`
#' @param color_by Either "type" or "cluster". "type" colors by the type 
#' (cell or gene) while "cluster" colors by the assigned cluster.
#' @param metadata optional. data.frame that should have either gene or cell names
#' (or both) in as rownames or a column named `name` and a column with the same name as 
#' `color_by`.
#' @param type Either "scatter", "contour" or "hex".
#' @return 
#' ggplot of UMAP
#' 
#' @export
plot_biUMAP <- function(umap_coords,
                        color_by = "type",
                        metadata=NULL,
                        type = "scatter",
                        point_size = 1,
                        hex_n = 40,
                        contour_n = 0.02,
                        color_genes = FALSE){
  
  if(!is.null(metadata)){
    
    if(!is(metadata,"data.frame")){
      metadata <- as.data.frame(metadata)
    }
    stopifnot(color_by %in% colnames(metadata))
    
    if(!"name" %in% colnames(metadata)){
      metadata$name <- rownames(metadata)
    }
    
    sel <- metadata$name %in% umap_coords$name
    metadata <- metadata[sel,]
    
    sel <- umap_coords$name %in% metadata$name
    matched_names <- match(umap_coords$name[sel], metadata$name)
    umap_coords[,color_by] <- "not_in_metadata"
    umap_coords[sel, color_by] <- as.character(metadata[matched_names, color_by])
    
  }
  
  cats <- length(unique(umap_coords[,color_by]))
  
  if (cats <= 9){
    
    colors <- RColorBrewer::brewer.pal(cats, "Set1")
    names(colors) <- sort(unique(umap_coords[,color_by]))
    
  } else if (cats <= 12) {
    
    colors <- RColorBrewer::brewer.pal(cats, "Set3")
    names(colors) <- sort(unique(umap_coords[,color_by]))
    
    
  } else {
    colors <- Polychrome::createPalette(N = cats,
                                        seedcolors = c("#00ffff", "#ff00ff", "#ffff00"), 
                                        range = c(10, 60))
    names(colors) <- sort(unique(umap_coords[,color_by]))
  }
  
  
  if (type == "contour" ){
    
    
    umap_cells <- dplyr::filter(umap_coords, type == "cell")
    umap_genes <- dplyr::filter(umap_coords, type == "gene")    
    
    xrange <- max(umap_cells$x)-min(umap_cells$x)
    yrange <- max(umap_cells$y)-min(umap_cells$y)
    
    if (isTRUE(color_genes)){
      color_by_genes <- color_by
      gene_colors <- colors
    } else {
      color_by_genes <- "type"
      
      if (is(color_genes, "character")){
        gene_colors <- c("gene" = color_genes)
      } else {
        gene_colors <- c("gene" = "#7393B3")
        
      }
    }
    
    interact <- paste0(
      "Type: ", umap_genes$type, "\n",
      "Name: ", umap_genes$name, "\n",
      "Cluster: ", umap_genes$cluster
    )
    p <- ggplot() +
      geom_density_2d(data = umap_cells,
                      mapping = aes_(x = ~x,
                                     y = ~y,
                                     colour = as.name(color_by)),
                      contour_var = "ndensity",
                      breaks = seq(0, 1.0, length.out = contour_n),
      ) +
      geom_point(data = umap_genes,
                 mapping = aes_(x = ~x,
                                y = ~y,
                                fill = as.name(color_by_genes),
                                text = quote(interact)
                 ),
                 color = "black",
                 size = point_size,
                 shape = 21,
                 # stroke = 0.25,
                 alpha = 1 ) +
      scale_fill_manual(values = gene_colors) +
      scale_color_manual(values = colors) +
      theme_bw()
    
  } else if (type == "hex"){
    
    umap_cells <- dplyr::filter(umap_coords, type == "cell")
    umap_genes <- dplyr::filter(umap_coords, type == "gene")
    
    xrange <- max(umap_cells$x)-min(umap_cells$x)
    yrange <- max(umap_cells$y)-min(umap_cells$y)
    
    bin_size <- c(xrange/hex_n, yrange/hex_n)
    
    if (isTRUE(color_genes)){
      color_by_genes <- color_by
    } else {
      color_by_genes <- "type"
      
      if (is(color_genes, "character")){
        colors<- c(colors, "gene" = color_genes)
      } else {
        colors<- c(colors, "gene" = "#7393B3")
        
      }
    }
    
    interact <- paste0(
      "Type: ", umap_genes$type, "\n",
      "Name: ", umap_genes$name, "\n",
      "Cluster: ", umap_genes$cluster
    )
    p <- ggplot() +
      geom_hex(data = umap_cells,
               mapping = aes_(x = ~x,
                              y = ~y,
                              fill = as.name(color_by),
                              alpha = quote(..count..)),
               # bins = bin_n,
               binwidth = bin_size,
               alpha = 0.7,
               color = "black") +
      geom_point(data = umap_genes,
                 mapping = aes_(x = ~x,
                                y = ~y,
                                fill = as.name(color_by_genes),
                                text = quote(interact)
                 ),
                 color = "black",
                 shape = 21,
                 size = point_size,
                 # stroke = 0.25,
                 alpha = 1 ) +
      scale_fill_manual(values = colors) +
      theme_bw()
    
    
  } else if (type == "scatter"){
    
    interact <- paste0(
      "Type: ", umap_coords$type, "\n",
      "Name: ", umap_coords$name, "\n",
      "Cluster: ", umap_coords$cluster)
    
    p <- ggplot(umap_coords, aes_(x=~x, y=~y, color = as.name(color_by),
                                  text = quote(interact))) +
      geom_point(alpha = 0.7, size = point_size) +
      scale_color_manual(values = colors) +
      theme_bw()
  }
  
  
  
  return(p)
  
}

#' plot biUMAP with gene expression
#' 
#' @param umap_coords
#' @param sce
#' @param feature
#' @param color_cells_by
#' @param assay
#' 
feature_biUMAP <- function(umap_coords, sce, feature = NULL, color_cells_by="expression", assay = "logcounts"){
  stopifnot(length(feature)<=1)
  
  if(color_cells_by == "expression") {
    isExpr <- FALSE
    if(!is.null(feature)) lgnd <- feature
  }else{
    isExpr <- TRUE
    lgnd <- color_cells_by
  }
  
  cell_idx <- which(umap_coords$type == "cell")
  
  if(!is.null(feature)){
    
    stopifnot(isTRUE(feature %in% umap_coords$name))
    cnts <- SummarizedExperiment::assay(sce, assay)
    umap_coords$expression <- NA
    umap_coords[cell_idx,]$expression <- cnts[feature, umap_coords$name[cell_idx]]
    
  }
  
  
  ggplot()+
    geom_point(umap_coords[umap_coords$type == "gene",],
               mapping=aes_(~x, ~y, text = paste0(
                                       "Type: ", quote(type), "\n",
                                       "Name: ", quote(name), "\n",
                                       "Cluster: ", quote(cluster))),color ="grey", alpha = 0.5) +
    geom_point(umap_coords[umap_coords$type == "cell",],
               mapping=aes_(~x, ~y, color = as.name(color_cells_by), text = paste0(
                                                    "Type: ", quote(type), "\n",
                                                    "Name: ", quote(name), "\n",
                                                    "Cluster: ", quote(cluster)))) +
    geom_point(data = na.omit(umap_coords[feature,c("name", "x","y")]),
               aes_(~x, ~y),
               color = "red") +
    geom_text_repel(data = na.omit(umap_coords[feature,c("name", "x","y")]),
                    aes_(~x, ~y, label= ~name),
                    color = "red") +
    viridis::scale_color_viridis(name=lgnd, discrete = isExpr) +
    theme_bw()
    
}

#' Shuffle rows of a data frame for better plotting.
#' @param df data.frame
#' @export
mix <- function(df){
  df <- df[sample(seq_len(nrow(df)), size = nrow(df)),]
}
  
#' Plot of 2D CA projection of the data.
#'
#' @description
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#'
#' @details
#' Choosing type "plotly" will generate an interactive html plot with the 
#' package plotly.
#' Type "ggplot" generates a static plot.
#' Depending on whether `princ_coords` is set to 1 or 2 either
#' the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted 
#' (assymetric biplot).
#' Labels for rows and columns should be stored in the row and column names 
#' respectively.
#' @return
#' Plot of class "plotly" or "ggplot".
#'
#' @param caobj An object of class "cacomp" with the relevant standardized and 
#' principal coordinates calculated,
#'  or alternatively an object of class "Seurat" or "SingleCellExperiment" 
#'  with a dim. reduction named "CA" saved.
#' @param caclust caclust object containing clustering results.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param princ_coords Integer. If 1 then principal coordinates are used for 
#' the rows,
#' if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label 
#' should be added
#' (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label 
#' should be added
#' (label should be stored in colnames).
#' Default NULL (no columns).
#' @param type String. Type of plot to draw. Either "ggplot" or "plotly". 
#' Default "plotly".
#' @param ... Further arguments.
#' @export
setGeneric("bicplot", function(caobj,
                               caclust,
                                 xdim = 1,
                                 ydim = 2,
                                 coords = 1,
                                 row_labels = NULL,
                                 col_labels = NULL,
                                 type = "plotly",
                                 ...) {
  standardGeneric("bicplot")
})


#TODO lots of new dependencies. Better way?
#' @rdname bicplot
#' @export
setMethod(f = "bicplot",
          signature(caobj = "cacomp", caclust = "caclust"),
          function(caobj, 
                   caclust,
                   xdim = 1,
                   ydim = 2,
                   coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   type = "plotly",
                   ...){
            
            if (!is(caobj,"cacomp")){
              stop("Not a CA object. Please run cacomp() first!")
            }
            
            ngenes <- nrow(caobj@std_coords_rows)
            ncells <- nrow(caobj@std_coords_cols)
            gene.idx <- which(rownames(caobj@prin_coords_rows) %in% names(gene_clusters(caclust)))
            cell.idx <- which(rownames(caobj@prin_coords_cols) %in% names(cell_clusters(caclust)))
            
            cells <- data.frame(clusters = rep(NA, ncells)) 
            genes <- data.frame(clusters = rep(NA, ngenes))
            
            cells[cell.idx,] <- as.vector(cell_clusters(caclust))
            genes[gene.idx,] <- as.vector(gene_clusters(caclust))
            
            if (coords == 1){
              
              if(sum(!is.null(caobj@prin_coords_rows), !is.null(caobj@std_coords_cols)) != 2){
                stop("Principal and/or standard coordinates not found, ",
                     "please run ca_coords() first!")
              }
              rows <- cbind(caobj@prin_coords_rows, genes)
              cols <- cbind(caobj@std_coords_cols, cells)
            } else if (coords == 2){
              if(sum(!is.null(caobj@prin_coords_cols), !is.null(caobj@std_coords_rows)) != 2){
                stop("Principal and/or standard coordinates not found, ",
                     "please run ca_coords() first!")
              }
              rows <- cbind(caobj@std_coords_rows, genes)
              cols <- cbind(caobj@prin_coords_cols, cells)
            }else if (coords == 3){
              if(sum(!is.null(caobj@U), !is.null(caobj@V)) != 2){
                stop("Singular eigenvectors not found, ",
                     "please run ca_coords() first!")
              }
              rows <- cbind(caobj@U, genes)
              cols <- cbind(caobj@V, cells)
            } else {
              stop("princ_coords must be either 1 for rows or 2 for columns.")
            }
            
            rows <- as.data.frame(rows)
            cols <- as.data.frame(cols)
            if (type == "ggplot"){
              
              # rows <- as.data.frame(rows)
              # cols <- as.data.frame(cols)
              #
              rnmx <- colnames(rows)[xdim]
              rnmy <- colnames(rows)[ydim]
              cnmx <- colnames(cols)[xdim]
              cnmy <- colnames(cols)[ydim]
              
              p <- ggplot2::ggplot()+
                ggplot2::geom_point(data=rows,
                                    ggplot2::aes_(x = as.name(rnmx), y = as.name(rnmy),
                                                  colour= rows$clusters),
                                    # colour = "#0066FF",
                                    alpha = 0.7, 
                                    shape = 1) +
                ggplot2::geom_point(data=cols,
                                    ggplot2::aes_(x = as.name(cnmx), y = as.name(cnmy),
                                                  colour= cols$clusters),
                                    # colour = "#990000",
                                    shape = 4) +
                ggplot2::theme_bw()
              
              if (!is.null(row_labels)){
                p <- p +
                  ggplot2::geom_point(data=rows[row_labels,],
                                      ggplot2::aes_(x = as.name(rnmx),
                                                    y = as.name(rnmy),
                                                    colour= rows$clusters),
                                      # colour = "#FF0000",
                                      shape = 16) +
                  ggrepel::geom_text_repel(data=rows[row_labels,],
                                           ggplot2::aes_(x = as.name(rnmx),
                                                         y = as.name(rnmy),
                                                         colour= rows[row_labels,]$clusters,
                                                         label=rownames(rows[row_labels,])),
                                           # colour = "#FF0000",
                                           max.overlaps = Inf)
              }
              if (!is.null(col_labels)){
                p <- p +
                  ggplot2::geom_point(data=cols[col_labels,],
                                      ggplot2::aes_(x = as.name(cnmx),
                                                    y = as.name(cnmy),
                                                    colour= cols[col_labels,]$clusters),
                                      # colour = "#990000",
                                      shape = 1) +
                  ggrepel::geom_text_repel(data=cols[col_labels,],
                                           ggplot2::aes_(x = as.name(cnmx),
                                                         y = as.name(cnmy),
                                                         colour= cols[col_labels,]$clusters,
                                                         label=rownames(cols[col_labels,])),
                                           # colour = "#990000",
                                           max.overlaps = Inf)
              }
            } else if (type == "plotly"){
              p <- plotly::plot_ly(type='scatter',
                                   source='plot2D',
                                   mode='markers') %>%
                plotly::add_trace(x = cols[,xdim],
                                  y = cols[,ydim],
                                  mode = 'markers',
                                  text = rownames(cols),
                                  textposition = "left",
                                  opacity = 1,
                                  marker = list(color = cols$clusters, # '#990000',
                                                symbol = 'x',
                                                size = 5),
                                  name = 'Columns',
                                  hoverinfo = 'text',
                                  type = 'scatter') %>%
                plotly::add_trace(x = rows[,xdim],
                                  y = rows[,ydim],
                                  mode = 'markers',
                                  text = rownames(rows),
                                  opacity = 0.7,
                                  marker = list(color = rows$clusters, #'#0066FF',
                                                symbol = 'circle-open',
                                                size = 2.5),
                                  name = 'genes',
                                  hoverinfo = 'text',
                                  type = 'scatter')
              
              if (!is.null(row_labels)){
                p <- p %>%
                  plotly::add_trace(x = rows[row_labels, xdim],
                                    y = rows[row_labels, ydim],
                                    mode = 'markers+text',
                                    text = rownames(rows)[row_labels],
                                    textposition = "left",
                                    textfont = list(color= rows$clusters), #'#FF0000'),
                                    marker = list(symbol = 'circle',
                                                  color ='#FF0000',
                                                  size = 5),
                                    name = 'marked row(s)',
                                    hoverinfo = 'text',
                                    type = 'scatter')
              }
              
              if (!is.null(col_labels)){
                p <- p %>%
                  plotly::add_trace(x = cols[col_labels, xdim],
                                    y = cols[col_labels, ydim],
                                    mode = 'markers+text',
                                    text = rownames(cols)[col_labels],
                                    textposition = "left",
                                    textfont = list(color= cols$clusters), #'#990000'),
                                    marker = list(symbol = 'circle-open',
                                                  color = cols$clusters, #'#990000',
                                                  size = 6.5),
                                    name = 'marked column(s)',
                                    hoverinfo = 'text',
                                    type = 'scatter')
              }
              
              p <- p %>%
                plotly::layout(autosize = TRUE,
                               title = '2D CA plot',
                               showlegend = FALSE,
                               xaxis = list(title = paste0('Dim', xdim)),
                               yaxis = list(title = paste0('Dim', ydim)))
              
            }
            
            return(p)
            
          })



