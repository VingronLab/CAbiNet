

#' Plot biUMAP
#' 
#' @param umap_coords data frame as outputted by `run_biUMAP_*`
#' @param color_by Either "type" or "cluster". "type" colors by the type 
#' (cell or gene) while "cluster" colors by the assigned cluster.
#' 
#' @return 
#' ggplot of UMAP
#' 
#' @export
plot_biUMAP <- function(umap_coords, color_by = "type"){
  
  if(color_by == "type"){
    p <- ggplot(umap_coords, aes(x=x, y=y, color = type,
                                 text = paste0(
                                   "Type: ", type, "\n",
                                   "Name: ", name, "\n",
                                   "Cluster: ", cluster))) +
      geom_point(alpha = 0.4) +
      theme_bw()
  } else if (color_by == "cluster"){
    
    p <- ggplot(umap_coords, aes(x=x, y=y, color = cluster,
                                 text = paste0(
                                   "Type: ", type, "\n",
                                   "Name: ", name, "\n",
                                   "Cluster: ", cluster)))+
      geom_point(alpha = 0.4) +
      theme_bw()
  } else {
    stop("color_by has to be either 'type' or 'cluster'.")
  }
  
  
  return(p)
  
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
#' @param obj An object of class "cacomp" with the relevant standardized and 
#' principal coordinates calculated,
#'  or alternatively an object of class "Seurat" or "SingleCellExperiment" 
#'  with a dim. reduction named "CA" saved.
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
#' @examples
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Run correspondence analysis
#' ca <- cacomp(obj = cnts, princ_coords = 3)
#'
#' ca_biplot(ca)
setGeneric("bicplot", function(obj,
                                 xdim = 1,
                                 ydim = 2,
                                 princ_coords = 1,
                                 row_labels = NULL,
                                 col_labels = NULL,
                                 type = "plotly",
                                 ...) {
  standardGeneric("bicplot")
})



#' @rdname bicplot
#' @export
setMethod(f = "bicplot",
          signature=(obj="cabicomp"),
          function(obj,
                   xdim = 1,
                   ydim = 2,
                   princ_coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   type = "plotly",
                   ...){
            
            if (!is(obj,"cacomp")){
              stop("Not a CA object. Please run cacomp() first!")
            }
            
            if (princ_coords == 1){
              
              if(sum(!is.null(obj@prin_coords_rows), !is.null(obj@std_coords_cols)) != 2){
                stop("Principal and/or standard coordinates not found, ",
                     "please run ca_coords() first!")
              }
              rows <- obj@prin_coords_rows
              cols <- obj@std_coords_cols
            } else if (princ_coords == 2){
              if(sum(!is.null(obj@prin_coords_cols), !is.null(obj@std_coords_rows)) != 2){
                stop("Principal and/or standard coordinates not found, ",
                     "please run ca_coords() first!")
              }
              rows <- obj@std_coords_rows
              cols <- obj@prin_coords_cols
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
                                    ggplot2::aes_(x = as.name(rnmx), y = as.name(rnmy)),
                                    # colour = "#0066FF",
                                    alpha = 0.7, 
                                    shape = 1) +
                ggplot2::geom_point(data=cols,
                                    ggplot2::aes_(x = as.name(cnmx), y = as.name(cnmy)),
                                    colour = "#990000",
                                    shape = 4) +
                ggplot2::theme_bw()
              
              if (!is.null(row_labels)){
                p <- p +
                  ggplot2::geom_point(data=rows[row_labels,],
                                      ggplot2::aes_(x = as.name(rnmx),
                                                    y = as.name(rnmy)),
                                      colour = "#FF0000",
                                      shape = 16) +
                  ggrepel::geom_text_repel(data=rows[row_labels,],
                                           ggplot2::aes_(x = as.name(rnmx),
                                                         y = as.name(rnmy),
                                                         label=rownames(rows[row_labels,])),
                                           colour = "#FF0000",
                                           max.overlaps = Inf)
              }
              if (!is.null(col_labels)){
                p <- p +
                  ggplot2::geom_point(data=cols[col_labels,],
                                      ggplot2::aes_(x = as.name(cnmx),
                                                    y = as.name(cnmy)),
                                      colour = "#990000",
                                      shape = 1) +
                  ggrepel::geom_text_repel(data=cols[col_labels,],
                                           ggplot2::aes_(x = as.name(cnmx),
                                                         y = as.name(cnmy),
                                                         label=rownames(cols[col_labels,])),
                                           colour = "#990000",
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
                                  marker = list(color = '#990000',
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
                                  marker = list(color ='#0066FF',
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
                                    textfont = list(color='#FF0000'),
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
                                    textfont = list(color='#990000'),
                                    marker = list(symbol = 'circle-open',
                                                  color ='#990000',
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



