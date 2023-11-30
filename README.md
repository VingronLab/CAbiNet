
# CAbiNet
**Correspondence Analysis based Biclustering on Networks**

This package provides functions to for the visualization and biclustering of single-cell RNA-seq data.

A longer vignette explaining how to install and use the package can be found here:
https://vingronlab.github.io/CAbiNet/

## Installation

You can install the package with:

``` r

# please also install this version of APL from github:
devtools::install_github("VingronLab/APL", ref = "cabinet_freeze")
devtools::install_github("VingronLab/CAbiNet")
```

> ⚠️ A bug in ggplot2 versions >3.3.0 and <3.4.1 lead to incorrect behaviour  in `plot_hex_biMAP`. Please make sure you have ggplot2 3.4.1 or higher installed. 

## Quick start

Here we provide a very short example of how to use the package. We hope to provide a more detailed description of how to use CAbiNet to perform your analysis in the near future.

``` r
library(CAbiNet)
library(APL)
library(scRNAseq)


sce <- DarmanisBrainData()

# Here you might want to do some preprocessing.

# Correspondence Analysis
caobj = cacomp(sce,
               dims = 50,
               top = 1000, # number of genes with highest inertia to keep.
               python = TRUE)

# SNN graph & biclustering
cabic <- caclust(obj = caobj,
                 k = 10,
                 loops = FALSE,
                 SNN_prune = 1/15,
                 mode = "all",
                 select_genes = TRUE,
                 prune_overlap = TRUE,
                 overlap = 0.2,
                 calc_gene_cell_kNN = FALSE,
                 resolution = 1,
                 algorithm = 'leiden')

sce$cabinet <- cell_clusters(cabic)

cabic <- biMAP(cabic, k = 30)

# plot results
plot_biMAP(cabic, color_genes = TRUE)

# Interactive biMAP where you can mouse over the points to see their identities
plot_biMAP(cabic, color_by = "cluster",
           color_genes = TRUE,
           interactive = TRUE)

plot_scatter_biMAP(cabic,
                   gene_alpha = 0,
                   color_by = "cell.type",
                   meta_df = colData(sce))


```
