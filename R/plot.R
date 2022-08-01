
#' Mixes the colors of two clusters proportionally.
#' 
#' @param df data.frame of cells with clusters in `color_by` and assigned 
#' hex bin in `hexbin`.
#' @param colors colors to be mixed.
#' @param cell Which hexbin to mix colors in.
#' @param color_by Column name where the clusters/groups are stored in `df`.
#' 
#' @returns 
#' Mixed color as hex code.
#' 
mix_rgb <- function(df, colors, cell, color_by){
  rgbcols <- col2rgb(colors)
  
  sel <- which(df$hexbin == cell)
  n_clust <- df[sel,color_by]
  n_clust <- table(as.character(n_clust))
  prop <- as.numeric(n_clust)
  names(prop) <- names(n_clust)
  prop <- prop/sum(prop)
  
  rgb_new <- sweep(rgbcols[,names(prop), drop=FALSE], MARGIN =2, FUN = "*", prop)  
  rgb_new <- rowSums(rgb_new)
  rgb_new <- rgb(red = rgb_new["red"], 
                 green = rgb_new["green"],
                 blue = rgb_new["blue"], 
                 maxColorValue = 255)
  return(rgb_new)
}

#' biMAP plotter
#' 
#' @description 
#' This function plots the different variants of biMAP
#'
#' @rdname biMAP_plotter
#' @param caclust caclust object with biMAP stored.
#' @param color_by Either "type" or "cluster". "type" colors by the type 
#' (cell or gene) while "cluster" colors by the assigned cluster. If meta_df is
#' provided a column from the data.frame can be specified.
#' @param meta_df optional. data.frame that should have either gene or cell names
#' (or both) as rownames or a column named `name` and a column with the same name as 
#' `color_by`.
#' @param type Determines the kind of plot outputted. Either "scatter" 
#' (scatter plot), "contour" (density contour plot) or "hex" (binned hexagonals).
#' @param cell_size integer. Size of the cells in the plot.
#' @param gene_size integer. Size of the genes in the plot.
#' @param cell_alpha numeric. Alpha value of cells (between 0 and 1).
#' @param gene_alpha numeric. Alpha value of genes (between 0 and 1).
#' @param hex_n integer. Number of hexagonal bins that partition the x-axis.
#' @param min_bin integer. Only bins with more than `min_bin` cells are plotted.
#' @param show_density logical. Should the alpha value for plot type "hex" be
#' adjusted by the number of cells in the bin?
#' @param contour_n integer. Number of contour lines to plot.
#' @param color_genes logical. If TRUE colors genes by color_by.
#' @param label_groups logical. If TRUE puts the group label on the median 
#' coordinates of the point.
#' @param group_label_size integer. Size of the group label.
#' @param label_marker_genes logical. If TRUE, names of genes are displayed.
#' 
#' @returns
#' biMAP plot as ggplot object.
#' 
biMAP_plotter <- function(caclust,
                          color_by = "cluster",
                          meta_df = NULL,
                          type = "scatter",
                          cell_size = 1,
                          gene_size = 3,
                          cell_alpha = 0.5,
                          gene_alpha = 1,
                          hex_n = 40,
                          min_bin = 2,
                          show_density = FALSE,
                          contour_n = 5,
                          color_genes = FALSE,
                          label_groups = TRUE,
                          group_label_size = 4,
                          label_marker_genes = FALSE){
  
  stopifnot(is(caclust, "caclust"))
  
  if (is.empty(caclust@bimap)){
    stop("Please run biMAP() first!")
  }
  
  umap_coords <- caclust@bimap
  
  if(!is.null(meta_df)){
    
    if(!is(meta_df,"data.frame")){
      meta_df <- as.data.frame(meta_df)
    }
    if(color_by %in% c(colnames(meta_df), colnames(umap_coords))){
      stop('color_by is found both from meta_df and bimap coordinates dataframe at the same time, 
           plot_biMAP is confused about which to use, please rename the column name in meta_df.')
    }
    stopifnot((color_by %in% colnames(meta_df)) | (color_by %in% colnames(umap_coords)))
    
    if(!"name" %in% colnames(meta_df)){
      meta_df$name <- rownames(meta_df)
    }
    
    sel <- meta_df$name %in% umap_coords$name
    meta_df <- meta_df[sel,]
    
    sel <- umap_coords$name %in% meta_df$name
    matched_names <- match(umap_coords$name[sel], meta_df$name)
    umap_coords[,color_by] <- "not_in_meta_df"
    umap_coords[sel, color_by] <- as.character(meta_df[matched_names, color_by])
    
  }
  
  cats <- length(unique(umap_coords[,color_by]))
  
  if (cats <= 9){
    
    colors <- suppressWarnings(RColorBrewer::brewer.pal(cats, "Set1"))
    colors <- colors[seq_len(cats)]
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
  
  umap_cells <- dplyr::filter(umap_coords, type == "cell")
  umap_genes <- dplyr::filter(umap_coords, type == "gene")   
  
  interact_cells <- paste0(
    "Type: ", umap_cells$type, "\n",
    "Name: ", umap_cells$name, "\n",
    "Cluster: ", umap_cells$cluster)
  
  
  interact_genes <- paste0(
    "Type: ", umap_genes$type, "\n",
    "Name: ", umap_genes$name, "\n",
    "Cluster: ", umap_genes$cluster)
  
  if (type == "contour" ){
    
    suppressWarnings({
      p <- contour_plot(umap_cells = umap_cells,
                        umap_genes = umap_genes,
                        color_by = color_by,
                        color_genes = color_genes,
                        colors = colors,
                        interact_genes = interact_genes,
                        gene_size = gene_size,
                        gene_alpha = gene_alpha,
                        contour_n = contour_n)
    })
    
  } else if (type == "hex"){
    
    suppressWarnings({
      p <- hex_plot(umap_cells = umap_cells,
                    umap_genes = umap_genes,
                    color_by = color_by,
                    colors = colors,
                    color_genes = color_genes,
                    hex_n = hex_n,
                    min_bin = min_bin,
                    show_density = show_density,
                    cell_alpha = cell_alpha,
                    gene_size = gene_size,
                    gene_alpha = gene_alpha,
                    interact_genes = interact_genes)
    })
    
  } else if (type == "scatter"){
    
    suppressWarnings({
      p <- scatter_plot(umap_cells = umap_cells, 
                        umap_genes = umap_genes,
                        color_by = color_by,
                        colors = colors,
                        color_genes = color_genes,
                        cell_size = cell_size,
                        cell_alpha = cell_alpha,
                        gene_size = gene_size,
                        gene_alpha = gene_alpha,
                        interact_genes = interact_genes,
                        interact_cells = interact_cells)
    })
    
  } else {
    stop("type must be either 'scatter', 'hex' or 'contour'")
  }
  

  if(label_groups){
    
    label_coords <- umap_cells %>%
      dplyr::group_by_at(color_by) %>%
      dplyr::summarize(fraction_of_group = dplyr::n(),
                       median_x = stats::median(x = x),
                       median_y = stats::median(x = y)) %>%
      dplyr::mutate(label = as.character(.data[[color_by]]))
    
    p <- p + ggrepel::geom_text_repel(data = label_coords,
                                      mapping = ggplot2::aes(x = median_x,
                                                             y = median_y,
                                                             label = label),
                                      size=I(group_label_size),
                                      fontface = "bold")
  }

  if(label_marker_genes){
    
    p <- p + ggrepel::geom_text_repel(data = umap_genes,
                                      mapping = ggplot2::aes(x = x,
                                                             y = y,
                                                             label = name),
                                      max.overlaps = 12,
                                      size=I(group_label_size-1))
  }
  p <- p + theme_bw()
  
  
  return(p)
  
}


#' internal helper to plot scatter biMAP
#' 
#' @param umap_cells biMAP coordinates data.frame subsetted to cells.
#' @param umap_genes biMAP coordinates data.frame subsetted to genes.
#' @param colors named vector of of length equal to the group size.
#' @param interact_cells character string with text and labels for cells in 
#' interactive plot.
#' @param interact_genes character string with text and labels for genes in 
#' interactive plots.
#' @inheritParams biMAP_plotter
#'
#' @returns
#' ggplot object
scatter_plot <- function(umap_cells, 
                         umap_genes,
                         color_by,
                         colors,
                         color_genes,
                         cell_size,
                         cell_alpha,
                         gene_size,
                         gene_alpha,
                         interact_genes,
                         interact_cells){
  
  
  if (isTRUE(color_genes)){
    color_by_genes <- color_by
    gene_colors <- colors
    
  } else {
    
    color_by_genes <- "type"
    
    if("not_in_meta_df" %in% names(colors)){
      cidx <- which(umap_cells[,color_by] == "not_in_meta_df")
      gidx <- which(umap_genes[,color_by] == "not_in_meta_df")
      
      if(length(cidx) == 0 & length(gidx) == nrow(umap_genes)){
        colors <- colors[-which(names(colors) == "not_in_meta_df")]
      }
    }
    
    if (is(color_genes, "character")){
      gene_colors <- c("gene" = color_genes)
    } else {
      gene_colors <- c("gene" = "#7393B3")
      
    }
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = umap_cells,
               mapping = ggplot2::aes_(x = ~x,
                                       y = ~y,
                                       color = as.name(color_by),
                                       text = quote(interact_cells)),
               alpha = cell_alpha,
               size = cell_size) +
    ggplot2::geom_point(data = umap_genes,
               mapping = ggplot2::aes_(x = ~x,
                                       y = ~y,
                                       fill = as.name(color_by_genes),
                                       text = quote(interact_genes)),
               color = "black",
               shape = 21,
               alpha = gene_alpha, 
               size = gene_size) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = gene_colors) +
    ggplot2::labs(x="Dim 1",
         y="Dim 2") +
    ggplot2::theme_bw()
  
  return(p)
  
}



#' internal helper to plot hex biMAP
#' 
#' @inheritParams biMAP_plotter
#' @inheritParams scatter_plot
#' 
#' @returns 
#' ggplot object
#' 
hex_plot <- function(umap_cells,
                     umap_genes,
                     color_by,
                     colors,
                     color_genes,
                     hex_n,
                     min_bin,
                     show_density,
                     cell_alpha,
                     gene_size,
                     gene_alpha,
                     interact_genes){
  
  xrange <- max(umap_cells$x)-min(umap_cells$x)
  yrange <- max(umap_cells$y)-min(umap_cells$y)
  
  bin_size <- c(xrange/hex_n, yrange/hex_n)
  
  if (isTRUE(color_genes)){
    color_by_genes <- color_by
  } else {
    color_by_genes <- "type"
    
    if("not_in_meta_df" %in% names(colors)){
      cidx <- which(umap_cells[,color_by] == "not_in_meta_df")
      gidx <- which(umap_genes[,color_by] == "not_in_meta_df")
      
      if(length(cidx) == 0 & length(gidx) == nrow(umap_genes)){
        colors <- colors[-which(names(colors) == "not_in_meta_df")]
      }
    }
    
    if (is(color_genes, "character")){
      colors<- c(colors, "gene" = color_genes)
    } else {
      colors<- c(colors, "gene" = "#7393B3")
      
    }
  }
  
  hexb <- hexbin::hexbin(umap_cells$x,
                         umap_cells$y,
                         xbins = hex_n,
                         xbnds = c(min(umap_cells$x),
                                   max(umap_cells$x)),
                         ybnds = c(min(umap_cells$y),
                                   max(umap_cells$y)),
                         IDs = TRUE)
  
  gghex <- data.frame(hexbin::hcell2xy(hexb),
                      count = hexb@count,
                      cell = hexb@cell,
                      xo = hexb@xcm,
                      yo = hexb@ycm,
                      hexclust = NA)
  
  for (i in seq_along(gghex$cell)){
    
    cell_id <- gghex$cell[i]
    hcnt <- gghex$count[i]
    
    orig_id <- which(hexb@cID == cell_id)
    umap_cells[orig_id,"hexbin"] <- cell_id
    
    gghex$hexclust[i] <- get_majority(umap_cells[orig_id, color_by])
    
  }
  
  
  hex_colors <- vector(mode = "character", length = length(gghex$cell))
  
  for (n in seq_along(gghex$cell)){
    hex_colors[n] <- mix_rgb(umap_cells,
                             colors = colors,
                             cell = gghex$cell[n],
                             color_by = color_by)
    
  }
  
  gghex$colors <- hex_colors
  
  if(min_bin > 0){
    gghex <- gghex[gghex$count >= min_bin,]
  }
  
  # names(hex_colors) <- as.character(gghex$cell)
  
  p <- ggplot2::ggplot()
  
  if (isTRUE(show_density)){
    
    p <- p + ggplot2::geom_hex(data = gghex,
                               mapping = ggplot2::aes(x = x,
                                                      y = y,
                                                      alpha = count),
                               fill = gghex$colors,
                               stat = "identity") + 
      ggplot2::scale_alpha(range = c(0.2, 1))
    
  } else {
    
    # gghex$std_count <- (gghex$count-min(gghex$count))/(max(gghex$count)-min(gghex$count))
    p <- p + ggplot2::geom_hex(data = gghex,
                               mapping = ggplot2::aes(x = x,
                                                      y = y),
                               fill = gghex$colors,
                               alpha = cell_alpha,
                               stat = "identity")
    
  }
  
  p <- p + ggplot2::geom_point(data = umap_genes,
                               mapping = ggplot2::aes_(x = ~x,
                                                       y = ~y,
                                                       fill = as.name(color_by_genes),
                                                       text = quote(interact_genes)
                               ),
                               color = "black",
                               shape = 21,
                               size = gene_size,
                               alpha = gene_alpha ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x="Dim 1",
                  y="Dim 2") +
    ggplot2::theme_bw()
  
  return(p)
}



#' internal helper to plot contour biMAP
#' 
#' @inheritParams biMAP_plotter
#' @inheritParams scatter_plot
#' @returns 
#' contour plot as ggplot object.
#' 
contour_plot <- function(umap_cells,
                         umap_genes,
                         color_by,
                         color_genes,
                         colors,
                         interact_genes,
                         gene_size,
                         gene_alpha,
                         contour_n){
  
  xrange <- max(umap_cells$x)-min(umap_cells$x)
  yrange <- max(umap_cells$y)-min(umap_cells$y)
  
  if (isTRUE(color_genes)){
    color_by_genes <- color_by
    
  } else {
    color_by_genes <- "type"
    
    if("not_in_metadata" %in% names(colors)){
      cidx <- which(umap_cells[,color_by] == "not_in_meta_df")
      gidx <- which(umap_genes[,color_by] == "not_in_meta_df")
      
      if(length(cidx) == 0 & length(gidx) == nrow(umap_genes)){
        colors <- colors[-which(names(colors) == "not_in_meta_df")]
      }
    }
    
    if (is(color_genes, "character")){
      gene_colors <- c("gene" = color_genes)
    } else {
      gene_colors <- c("gene" = "#7393B3")
    }
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_density_2d(data = umap_cells,
                             mapping = ggplot2::aes_(x = ~x,
                                                     y = ~y,
                                                     colour = as.name(color_by)),
                             contour_var = "ndensity",
                             breaks = seq(0, 1.0, length.out = contour_n),
    ) +
    ggplot2::geom_point(data = umap_genes,
                        mapping = ggplot2::aes_(x = ~x,
                                                y = ~y,
                                                fill = as.name(color_by_genes),
                                                text = quote(interact_genes)
                        ),
                        color = "black",
                        size = gene_size,
                        shape = 21,
                        # stroke = 0.25,
                        alpha = gene_alpha ) +
    ggplot2::scale_fill_manual(values = gene_colors) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x="Dim 1",
                  y="Dim 2") +
    ggplot2::theme_bw()
  
  
  return(p)
  
}

#' Plotting biMAP embedding
#' @description
#' TODO
#' @param obj A caclust object or SingleCellExperiment object  
#' @param caclust_meta_name character. The name of caclust object stored in metadata(SingleCellExperiment object)
#' @inheritParams plot_biMAP
#' @details
#' TODO
#' @return
#' A caclust object or SingleCellExperiment object.

#' @export
setGeneric("plot_biMAP", function(obj,
                                  caclust_meta_name = 'caclust',
                                  meta_df = NULL,
                                  color_by = "cluster",
                                  type = "scatter",
                                  cell_size = 1,
                                  gene_size = 3,
                                  hex_n = 40,
                                  min_bin = 2,
                                  contour_n = 5,
                                  cell_alpha = 0.5,
                                  gene_alpha = 1,
                                  show_density = FALSE,
                                  color_genes = FALSE,
                                  label_groups = TRUE,
                                  group_label_size=4,
                                  labels_per_group=1,
                                  label_marker_genes = FALSE,
                                  ...){
  standardGeneric("plot_biMAP")
})

#
#' @rdname plot_biMAP
#' @param color_by character which can be chosen from 'type', 'cluster' and column in the input meta_df
#' @param meta_df data.frame.
#' @inheritParams plot_biMAP
#' @export
setMethod(f = "plot_biMAP",
          signature(obj = "caclust"),
          function(obj, 
                   meta_df = NULL,
                   color_by = "cluster",
                   ...){
            
            p <- biMAP_plotter(caclust = obj,
                               color_by = color_by,
                               meta_df = meta_df,
                               ...)
            return(p)
            
          })

#
#' @rdname plot_biMAP
#' @param obj SingleCellExperiment object
#' @param caclust_meta_name character. Slot name of caclust object stored in meatadata(obj)
#' @param color_by character which can be chosen from 'type', 'cluster',column in the input meta_df and columns in colData of the obj (if meta_df == NULL)
#' @param meta_df data.frame.
#' @inheritParams plot_biMAP
#' @export
setMethod(f = "plot_biMAP",
          signature(obj = "SingleCellExperiment"),
          function(obj, 
                   caclust_meta_name = 'caclust',
                   meta_df = NULL,
                   color_by = "cluster",
                   type = 'scatter',
                   ...){
            
            if(isFALSE(caclust_meta_name %in% names(S4Vectors::metadata(obj)))){
              stop(paste('The aclust object with name', caclust_meta_name, 'is not found in metadata(sce), please try a different "biMAP_meta_name".'))
            }
            caclust <- S4Vectors::metadata(obj)[[caclust_meta_name]]
            
            if (is.null(meta_df) & (isFALSE(color_by %in%  colnames(caclust@bimap)))){
              meta_df = colData(obj)
            }
            
            if(isFALSE(color_by %in% c(colnames(meta_df), colnames(caclust@bimap)))){
              stop('color_by not found in either meta_df or obj')
            }
            
            p <- biMAP_plotter(caclust = caclust,
                               color_by = color_by,
                               meta_df = meta_df,
                               type = type,
                               ...)
            return(p)
            
          })




#' Plot biMAP with hexagonal bins
#' @family plot_*_biMAP
#' 
#' @description 
#' Bins cells into hexagons and colors them proportionally to their group label.
#' @param obj Object that stores biMAP coordinates. Can be either a "caclust"
#' object or of class "SingleCellExperiment".
#' @inheritParams biMAP_plotter
#' @param ... Further arguments
#' @returns 
#' hex-plot biMAP as a ggplot.
#' 
#' @export
setGeneric("plot_hex_biMAP", function(obj,
                                      meta_df = NULL,
                                      color_by = "cluster",
                                      gene_size = 2,
                                      cell_alpha = 0.8,
                                      gene_alpha = 0.8,
                                      hex_n = 40,
                                      min_bin = 2,
                                      show_density = FALSE,
                                      color_genes = FALSE,
                                      label_groups = TRUE,
                                      group_label_size=4,
                                      label_marker_genes = FALSE,
                                      ...) {
  standardGeneric("plot_hex_biMAP")
})

#' @rdname plot_hex_biMAP
#' @export
setMethod(f = "plot_hex_biMAP",
          signature(obj = "caclust"),
          function(obj,
                   meta_df = NULL,
                   color_by = "cluster",
                   gene_size = 2,
                   cell_alpha = 0.8,
                   gene_alpha = 0.8,
                   hex_n = 40,
                   min_bin = 2,
                   show_density = FALSE,
                   color_genes = FALSE,
                   label_groups = TRUE,
                   group_label_size=4,
                   label_marker_genes = FALSE,
                   ...){

  
  
  p <- biMAP_plotter(caclust = obj,
                     meta_df = meta_df,
                     color_by = color_by,
                     type = "hex",
                     cell_size = NA,
                     gene_size = gene_size,
                     hex_n = hex_n,
                     min_bin = min_bin,
                     contour_n = NA,
                     cell_alpha = cell_alpha,
                     gene_alpha = gene_alpha,
                     show_density = show_density,
                     color_genes = color_genes,
                     label_groups = label_groups,
                     group_label_size = group_label_size,
                     label_marker_genes = label_marker_genes)
  
  return(p)
})

#' @rdname plot_hex_biMAP
#' @export
setMethod(f = "plot_hex_biMAP",
          signature(obj = "SingleCellExperiment"),
          function(obj,
                   meta_df = NULL,
                   color_by = "cluster",
                   gene_size = 2,
                   cell_alpha = 0.8,
                   gene_alpha = 0.8,
                   hex_n = 40,
                   min_bin = 2,
                   show_density = FALSE,
                   color_genes = FALSE,
                   label_groups = TRUE,
                   group_label_size=4,
                   label_marker_genes = FALSE,
                   ...){
            
            
            
            p <- plot_biMAP( obj = obj,
                               meta_df = meta_df,
                               color_by = color_by,
                               type = "hex",
                               cell_size = NA,
                               gene_size = gene_size,
                               hex_n = hex_n,
                               min_bin = min_bin,
                               contour_n = NA,
                               cell_alpha = cell_alpha,
                               gene_alpha = gene_alpha,
                               show_density = show_density,
                               color_genes = color_genes,
                               label_groups = label_groups,
                               group_label_size = group_label_size,
                               label_marker_genes = label_marker_genes)
            
            return(p)
          })



#' Plots a biMAP scatter plot
#' @family plot_*_biMAP
#'
#' @description 
#' Plots cells and genes as points. Genes are plotted on top of the cells in
#' a larger size to make visual differentiation easier.
#' @param obj Object that stores biMAP coordinates. Can be either a "caclust"
#' object or of class "SingleCellExperiment".
#' @inheritParams biMAP_plotter
#' @param ... Further arguments.
#' @returns
#' biMAP scatter plot as ggplot object.
#'
#' @export
setGeneric("plot_scatter_biMAP", function(obj,
                                          meta_df = NULL,
                                          color_by = "cluster",
                                          cell_size = 1,
                                          gene_size = 2,
                                          cell_alpha = 0.7,
                                          gene_alpha = 0.8,
                                          color_genes = FALSE,
                                          label_groups = TRUE,
                                          group_label_size = 4,
                                          label_marker_genes = FALSE,
                                          ...) {
  standardGeneric("plot_scatter_biMAP")
})

#' @rdname plot_scatter_biMAP
#' @export
setMethod(f = "plot_scatter_biMAP",
          signature = (obj = "caclust"),
          function(obj,
                  meta_df = NULL,
                  color_by = "cluster",
                  cell_size = 1,
                  gene_size = 2,
                  cell_alpha = 0.7,
                  gene_alpha = 0.8,
                  color_genes = FALSE,
                  label_groups = TRUE,
                  group_label_size = 4,
                  label_marker_genes = FALSE,
                  ...){

  p <- biMAP_plotter(caclust = obj,
                      meta_df = meta_df,
                      color_by = color_by,
                      type = "scatter",
                      cell_size = cell_size,
                      gene_size = gene_size,
                      hex_n = NA,
                      min_bin = NA,
                      contour_n = NA,
                      cell_alpha = cell_alpha,
                      gene_alpha = gene_alpha,
                      show_density = NA,
                      color_genes = color_genes,
                      label_groups = label_groups,
                      group_label_size = group_label_size,
                      label_marker_genes = label_marker_genes)

  return(p)
})

#' @rdname plot_scatter_biMAP
#' @export
setMethod(f = "plot_scatter_biMAP",
          signature = (obj = "SingleCellExperiment"),
          function(obj,
                   meta_df = NULL,
                   color_by = "cluster",
                   cell_size = 1,
                   gene_size = 2,
                   cell_alpha = 0.7,
                   gene_alpha = 0.8,
                   color_genes = FALSE,
                   label_groups = TRUE,
                   group_label_size = 4,
                   label_marker_genes = FALSE,
                   ...){

            p <- plot_biMAP(obj = obj,
                               meta_df = meta_df,
                               color_by = color_by,
                               type = "scatter",
                               cell_size = cell_size,
                               gene_size = gene_size,
                               hex_n = NA,
                               min_bin = NA,
                               contour_n = NA,
                               cell_alpha = cell_alpha,
                               gene_alpha = gene_alpha,
                               show_density = NA,
                               color_genes = color_genes,
                               label_groups = label_groups,
                               group_label_size = group_label_size,
                               label_marker_genes = label_marker_genes)

            return(p)
            })


#' Plot contour biMAP
#' @family plot_*_biMAP
#'
#' @description 
#' Cell density is plotted as contour lines and genes are plotted as points
#' over the contours.
#' @param obj Object that stores biMAP coordinates. Can be either a "caclust"
#' object or of class "SingleCellExperiment".
#' @inheritParams biMAP_plotter
#' @param ... Further arguments.
#' @returns 
#' contour biMAP plot as ggplot object.
#' @export
setGeneric("plot_contour_biMAP", function(obj,
                                          meta_df = NULL,
                                          color_by = "cluster",
                                          gene_size = 2,
                                          gene_alpha = 0.8,
                                          contour_n = 5,
                                          color_genes = FALSE,
                                          label_groups = TRUE,
                                          group_label_size = 4,
                                          label_marker_genes = FALSE,
                                          ...) {
  standardGeneric("plot_contour_biMAP")
})

#' @rdname plot_contour_biMAP
#' @export
setMethod(f = "plot_contour_biMAP", 
          signature = (obj = "caclust"),
          function(obj,
                   meta_df = NULL,
                   color_by = "cluster",
                   gene_size = 2,
                   gene_alpha = 0.8,
                   contour_n = 5,
                   color_genes = FALSE,
                   label_groups = TRUE,
                   group_label_size = 4,
                   label_marker_genes = FALSE,
                   ...){
  
  p <- biMAP_plotter(caclust = obj,
                      meta_df = meta_df,
                      color_by = color_by,
                      type = "contour",
                      cell_size = NA,
                      gene_size = gene_size,
                      hex_n = NA,
                      min_bin = NA,
                      contour_n = contour_n,
                      cell_alpha = NA,
                      gene_alpha = gene_alpha,
                      show_density = FALSE,
                      color_genes = color_genes,
                      label_groups = label_groups,
                      group_label_size = group_label_size,
                      label_marker_genes = label_marker_genes)
  
  return(p)
  
})


#' @rdname plot_contour_biMAP
#' @export
setMethod(f = "plot_contour_biMAP", 
          signature = (obj = "SingleCellExperiment"),
          function(obj,
                   meta_df = NULL,
                   color_by = "cluster",
                   gene_size = 2,
                   gene_alpha = 0.8,
                   contour_n = 5,
                   color_genes = FALSE,
                   label_groups = TRUE,
                   group_label_size = 4,
                   label_marker_genes = FALSE,
                   ...){
            
            p <- plot_biMAP(obj = obj,
                               meta_df = meta_df,
                               color_by = color_by,
                               type = "contour",
                               cell_size = NA,
                               gene_size = gene_size,
                               hex_n = NA,
                               min_bin = NA,
                               contour_n = contour_n,
                               cell_alpha = NA,
                               gene_alpha = gene_alpha,
                               show_density = FALSE,
                               color_genes = color_genes,
                               label_groups = label_groups,
                               group_label_size = group_label_size,
                               label_marker_genes = label_marker_genes)
            
            return(p)
            
          })


#' plot biMAP with gene expression
#' @rdname plot_feature_biMAP
#' @param sce SinleCellExperiment object
#' @param caclust caclust object
#' @param feature character. Name of gene to visualize
#' @param color_cells_by character. Default: expression
#' @param assay character. Name of assay in SingleCellExperiment used for visualization.
#' @export
plot_feature_biMAP <- function(sce,
                               caclust, 
                              feature = NULL, 
                              color_cells_by="expression", 
                              assay = "logcounts"){
  
  stopifnot(is(caclust, "caclust"))
  umap_coords <- caclust@bimap
  
  stopifnot(is(sce, "SingleCellExperiment"))
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
    
    # stopifnot(isTRUE(feature %in% umap_coords$name))
    cnts <- SummarizedExperiment::assay(sce, assay)
    umap_coords$expression <- NA
    umap_coords[cell_idx,]$expression <- cnts[feature, umap_coords$name[cell_idx]]

    umap_coords[cell_idx,] <- umap_coords[cell_idx,][order(umap_coords[cell_idx,"expression"], decreasing = FALSE),]
  }
  
  
  ggplot()+
    geom_point(umap_coords[umap_coords$type == "gene",],
               mapping=aes_(~x, ~y, text = paste0(
                                       "Type: ", quote(type), "\n",
                                       "Name: ", quote(name), "\n",
                                       "Cluster: ", quote(cluster))), color ="#A9A9A9", alpha = 0.5) +  #grey
    geom_point(umap_coords[umap_coords$type == "cell",],
               mapping=aes_(~x, ~y, color = as.name(color_cells_by), text = paste0(
                                                    "Type: ", quote(type), "\n",
                                                    "Name: ", quote(name), "\n",
                                                    "Cluster: ", quote(cluster))),
               alpha = 0.8) +
    geom_point(data = na.omit(umap_coords[feature,c("name", "x","y")]),
               aes_(~x, ~y),
               color = "red") +
    geom_text_repel(data = na.omit(umap_coords[feature,c("name", "x","y")]),
                    aes_(~x, ~y, label= ~name),
                    color = "red") +
    viridis::scale_color_viridis(name=lgnd, discrete = isExpr) +
    labs(x="Dim 1",
         y="Dim 2")+
    theme_bw()
    
}


metadata_biMAP <- function(umap_coords, 
                           sce, 
                           color_cells_by, 
                           continous = FALSE){
  meta_df <- colData(sce)

  cell_idx <- which(umap_coords$type == "cell")

  umap_coords[,color_cells_by] <- NA
  umap_coords[cell_idx, color_cells_by] <- meta_df[umap_coords[cell_idx,]$name, color_cells_by] 



  p <- ggplot()+
    geom_point(umap_coords[umap_coords$type == "gene",],
               mapping=aes_(~x, ~y), color ="#A9A9A9", alpha = 0.5) +  #grey
    geom_point(umap_coords[umap_coords$type == "cell",],
               mapping=aes_(~x, ~y, color = as.name(color_cells_by)),
               alpha = 0.8)

    if(isTRUE(continous)){

      p <- p + viridis::scale_color_viridis(name=color_cells_by, discrete = FALSE)

    } 


    p<- p + 
      labs(x="Dim 1",
           y="Dim 2")+
      theme_bw()

    return(p)
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
                   rm.show = TRUE,
                   ...){
            
            if (!is(caobj,"cacomp")){
              stop("Not a CA object. Please run cacomp() first!")
            }
            
            ngenes <- nrow(caobj@std_coords_rows)
            ncells <- nrow(caobj@std_coords_cols)
            gene.idx <- which(rownames(caobj@prin_coords_rows) %in% names(gene_clusters(caclust)))
            cell.idx <- which(rownames(caobj@prin_coords_cols) %in% names(cell_clusters(caclust)))
            
            cells <- data.frame(clusters = rep('trimmed', ncells)) 
            genes <- data.frame(clusters = rep('trimmed', ngenes))
            
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
            
            if (isFALSE(rm.show)){
              rows <- rows[rows$cluster != 'trimmed',]
              cols <- cols[cols$cluster != 'trimmed',]
            }
            
            # rows <- rows[rownames(rows) %in% names(genes),]
            # cols <- cols[rownames(cols) %in% names(cells),]
            if (type == "ggplot"){
              
              # rows <- as.data.frame(rows)
              # cols <- as.data.frame(cols)
              #
              rnmx <- colnames(rows)[xdim]
              rnmy <- colnames(rows)[ydim]
              cnmx <- colnames(cols)[xdim]
              cnmy <- colnames(cols)[ydim]
              
              p <- ggplot2::ggplot()+
                ggplot2::geom_point(data=cols,
                                    ggplot2::aes_(x = as.name(cnmx), y = as.name(cnmy),
                                                  colour= cols$clusters),
                                    # colour = "#990000",
                                    shape = 4) +
                ggplot2::geom_point(data=rows,
                                    ggplot2::aes_(x = as.name(rnmx), y = as.name(rnmy),
                                                  colour= rows$clusters),
                                    # colour = "#0066FF",
					alpha = 1, 
                                    shape = 1) +
                ggplot2::theme_bw() #+ ggsci::scale_color_npg()
              
              if (!is.null(row_labels)){
                p <- p +
                  ggplot2::geom_point(data=rows[row_labels,],
                                      ggplot2::aes_(x = as.name(rnmx),
                                                    y = as.name(rnmy),
                                                    colour= rows$clusters[row_labels]),
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

#' BiMAP visualzation of expression of features
#' @rdname feature_biMAP
#' @param obj SinleCellExperiment object
#' @param caclust caclust object or NULL
#' @param biMAP_meta_name character. Name of biMAP coordiantes data.frame stored in metadata(sce obj)
#' @inheritParams plot_feature_biMAP
#' @export
setGeneric("feature_biMAP", function(obj,
                                     caclust = NULL,
                                     biMAP_meta_name = 'biMAP_SNNdist',
                                     feature = NULL, 
                                     color_cells_by="expression", 
                                     assay = "logcounts",
                                     ...){
  standardGeneric("feature_biMAP")
})

#
#' @rdname feature_biMAP
#' @param obj SinleCellExperiment object
#' @param caclust caclust with bimap slot
#' @inheritParams plot_feature_biMAP
#' @export
#' 
setMethod(f = "feature_biMAP",
          signature(obj = "SingleCellExperiment", caclust = "caclust"),
          function(obj, 
                   caclust,
                   feature = NULL, 
                   color_cells_by="expression", 
                   assay = "logcounts",
                   ...){
            
            p <- plot_feature_biMAP(sce = obj,
                                caclust = caclust,
                                feature = feature, 
                                color_cells_by=color_cells_by, 
                                assay =assay)
            return(p)
            
          })

#
#' @rdname feature_biMAP
#' @param obj SinleCellExperiment object
#' @param caclust caclust object with bimap slot
#' @param biMAP_meta_name character. Name of biMAP coordiantes data.frame stored in metadata(sce obj)
#' @inheritParams plot_feature_biMAP
#' @export
#' 
setMethod(f = "feature_biMAP",
          signature(obj = "SingleCellExperiment"),
          function(obj, 
                   caclust,
                   caclust_meta_name = 'caclust',
                   feature = NULL, 
                   color_cells_by="expression", 
                   assay = "logcounts",
                   ...){
            
            if(is.null(caclust)){
              if(caclust_meta_name %in% names(S4Vectors::metadata(obj))){
                caclust <- S4Vectors::metadata(obj)[[biMAP_meta_name]]
                # cat(paste0('Plotting by data.frame ', caclist_meta_name, ' from metadata(sce).\n'))
              }else{
              #   stop(paste('The biMAP coordinate data.frame with name', 
              #              biMAP_meta_name, 'is not found in metadata, please try a different "biMAP_meta_name".'))
              # }
                stop()
                }
            }
            
            p <- plot_feature_biMAP(sce = obj,
                                    caclust = caclust,
                                    feature = feature, 
                                    color_cells_by=color_cells_by, 
                                    assay =assay)
            return(p)
            
          })
