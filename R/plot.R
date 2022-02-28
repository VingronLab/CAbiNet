

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



