#! /usr/bin/Rscript

library(ggplot2)
library(enrichplot) 
library(cowplot)


########################################################################################
################################ --------> plots <---------- ###########################
# Source: https://yulab-smu.github.io/clusterProfiler-book/chapter12.html
# Gene Concept Network
fig.pathways.REACTOME <- function(reactome_res, entrezid_log2fc, showCategories, title,
                                  width, height, output.dir) 
{
  if (typeof(showCategories) != 'character') 
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
  
  if (all(entrezid_log2fc) < 0) 
  {
    low = 'red'
    mid = 'red'
    high = 'blue'
  } else {
    low = 'blue'
    mid = 'blue'
    high = 'red'
  }

  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = enrichplot::cnetplot(reactome_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + 
    ggtitle("REACTOME: Pathway Enrichment Analysis using DB Reactome") 
  # Color upper and lower border
  min.value <- floor( min(p1$data$color, na.rm = TRUE) )
  max.value <- ceiling( max(p1$data$color, na.rm = TRUE) )

  p1 <- p1 + scale_color_gradientn(
    name = "fold change", colours = c("blue", "red"), limits= c(min.value, max.value))
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
  
  p1 
}


#' Plot pathways in a dotplot
#' 
#' @param pathway_res : enrichobject
#' @param showCategories : str, int
#'   which or how many pathways to visualise; limited by max. number of enriched pathways
fig.pathway.dotplot <- function(pathway_res, showCategories, method, title,
                                width, height, output.dir) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(pathway_res)) 
    {
      showCategories = nrow(pathway_res)
    }
  } 
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = enrichplot::dotplot(pathway_res, showCategory = showCategories, color='p.adjust') + 
    ggtitle(paste(method, ": Pathway Enrichment Analysis", sep = " "))
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}
