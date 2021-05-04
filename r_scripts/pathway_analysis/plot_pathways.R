#! /usr/bin/Rscript

library(ggplot2)
library(enrichplot) 
library(cowplot)


########################################################################################
################################ --------> plots <---------- ###########################
# Source: https://yulab-smu.github.io/clusterProfiler-book/chapter12.html
# Gene Concept Network
fig.pathways.REACTOME <- function(reactome_res, entrezid_log2fc, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(reactome_res)) 
    {
      showCategories = nrow(reactome_res)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(reactome_res$Description, showCategories)
  }

  par(mfrow = c(1, 1))
  p1 = enrichplot::cnetplot(reactome_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "REACTOME: Pathway Enrichment Analysis using DB Reactome") +
    scale_colour_gradient2(name = "fold change", low = "green", mid = "white", high = "red")
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}

#' Plot pathways in a dotplot
#' 
#' @param pathway_res : enrichobject
#' @param showCategories : str, int
#'   which or how many pathways to visualise; limited by max. number of enriched pathways
fig.pathway.dotplot <- function(pathway_res, showCategories, method) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(pathway_res)) 
    {
      showCategories = nrow(pathway_res)
    }
  } 
  p1 = enrichplot::dotplot(pathway_res, showCategory = showCategories, color='p.adjust') + 
    ggtitle(paste(method, ": Pathway Enrichment Analysis", sep = " "))
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}

#' Visualise pathways in a network
#' 
#' @param pathway_res : enrichobject
#' @param genelist_log2FC : list of log2FC values received from DGE analysis
#' @param nodelabel choose between: category, all, gene, none
fig.cnetplot <- function(pathway_res, genelist_log2FC, nodelabel) 
{
  p1 = enrichplot::cnetplot(pathway_res, categorySize = "pvalue", foldChange = genelist_log2FC, 
                            circular = TRUE, colorEdge = TRUE, node_label = nodelabel)
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  
}

#' Plot pathways in a dotplot
#' 
#' @param pathway_res : enrichobject
#' @param showCategories : int
#'   how many pathways to visualise; limited by max. number of enriched pathways
fig.dotplot <- function(pathway_res, showCategories) 
{
  if (showCategories > nrow(pathway_res)) 
  {
    showCategories = nrow(pathway_res)
  }
  p1 = enrichplot::dotplot(pathway_res, showCategory = showCategories) 
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}

