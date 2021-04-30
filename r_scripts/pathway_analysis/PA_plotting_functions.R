
# Plot for pathway enrichemnt analysis:
library(ggplot2)
library(enrichplot) 
library(cowplot)


########################################################################################
################################ --------> plots <---------- ###########################
# Source: https://yulab-smu.github.io/clusterProfiler-book/chapter12.html
# Gene Concept Network
fig.pathways <- function(reactome_res, kegg_res, entrezid_log2fc, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(reactome_res)) 
    {
      showCategories = nrow(reactome_res)
    }
  } 
  p1 = enrichplot::cnetplot(reactome_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "REACTOME: Pathway Enrichment Analysis using DB Reactome")
  p2 = enrichplot::cnetplot(kegg_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "KEGG: Pathway Enrichment Analysis using DB KEGG")
  
  cowplot::plot_grid(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2])
}


fig.gsea_pathways <- function(gsea_res, entrezid_log2fc, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(gsea_res)) 
    {
      showCategories = nrow(gsea_res)
    }
  } 
  p1 = enrichplot::cnetplot(gsea_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "GSEA: Gene Set Enrichment Analysis")
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.go_enrichment <- function(clusterprofiler_res, entrezid_log2fc, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(clusterprofiler_res)) 
    {
      showCategories = nrow(clusterprofiler_res)
    }
  } 
  p1 = enrichplot::cnetplot(clusterprofiler_res, showCategory = showCategories,
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "ClusterProfiler: Cell Membrane GO-term Enrichment")
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.gsea_GO <- function(gsea.go_object, entrezid_log2fc, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(gsea.go_object)) 
    {
      showCategories = nrow(gsea.go_object)
    }
  } 
  p1 = enrichplot::cnetplot(gsea.go_object, showCategory = showCategories,
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "ClusterProfiler: GO-GSEA")
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.pathways.dose <- function(dose_res, entrezid_log2fc, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(dose_res)) 
    {
      showCategories = nrow(dose_res)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(dose_res$Description, showCategories)
  } 
  
  par(mfrow = c(1, 1))
  p1 = enrichplot::cnetplot(dose_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "DOSE: Disease Ontology Semantic and Enrichment")
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.pathways.KEGG <- function(kegg_res, entrezid_log2fc, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(kegg_res)) 
    {
      showCategories = nrow(kegg_res)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(kegg_res$Description, showCategories)
  }
  
  par(mfrow = c(1, 1))
  p1 = enrichplot::cnetplot(kegg_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + ggtitle(
                              "KEGG: Pathway Enrichment Analysis using DB KEGG")
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


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
                              "REACTOME: Pathway Enrichment Analysis using DB Reactome")
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.pathway.dotplots <- function(reactome_res, kegg_res, showCategories) 
{
  if (type(showCategories) != 'character') 
  {
    if (showCategories > nrow(reactome_res)) 
    {
      showCategories = nrow(reactome_res)
    }
  } 
  p1 = enrichplot::dotplot(reactome_res, showCategory = showCategories) + ggtitle(
    "REACTOME: Pathway Enrichment Analysis using DB Reactome")
  p2 = enrichplot::dotplot(kegg_res, showCategory = showCategories) + ggtitle(
    "KEGG: Pathway Enrichment Analysis using DB KEGG")
  
  cowplot::plot_grid(p1, p2, ncol = 2, nrow = 1)
}

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


fig.pubmed <- function(pathway_res, showCategories) 
{
  #' Display how many publications have been published in the last years according to the 
  #' found enriched pathways
  #' 
  #' @param pathway_res (Large enrichResult)
  
  if (showCategories > nrow(pathway_res)) 
  {
    showCategories = nrow(pathway_res)
  }
  
  # Plot published Publications until current year
  current_year <- c(as.integer(format(Sys.Date(), "%Y")))
  
  terms <- pathway_res$Description[1:showCategories]
  par(mfrow = c(1, 2))
  p <- pmcplot(terms, 2010:current_year)
  p2 <- pmcplot(terms, 2010:current_year, proportion = FALSE)
  cowplot::plot_grid(p, p2, ncol = 2)
}

fig.cnetplot <- function(pathway_res, genelist_log2FC, nodelabel) 
{
  #' 
  #' @param pathway_res
  #' @param genelist_log2FC
  #' @param nodelabel choose between: category, all, gene, none
  
  p1 = enrichplot::cnetplot(pathway_res, categorySize = "pvalue",
                            foldChange = genelist_log2FC, 
                            circular = TRUE, colorEdge = TRUE, node_label = nodelabel)
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  
}


fig.dotplot <- function(pathway_res, showCategories) 
{
  if (showCategories > nrow(pathway_res)) 
  {
    showCategories = nrow(pathway_res)
  }
  p1 = enrichplot::dotplot(pathway_res, showCategory = showCategories) 
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.emapplot <- function(pathway_res) 
{
  #' enrichment map 
  p1 <- emapplot(pathway_res, color = "pvalue", vertex.label.cex = 1.2, 
                 pie_scale = 1.5, layout = "kk") 
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.heatplot <- function(pathway_res, entrezid_log2FC, showCategories) 
{
  #' heatmap plot of enriched terms
  if (showCategories > nrow(pathway_res)) 
  {
    showCategories = nrow(pathway_res)
  }
  p1 <- enrichplot::heatplot(pathway_res, foldChange = entrezid_log2FC, 
                             showCategory = showCategories)
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.upsetplot <- function(pathway_res) 
{
  #' Upset plot
  p1 <- upsetplot(pathway_res)
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.ridgeplot <- function(pathway_res) 
{
  #' ridgeline plot for expression distribution of GSEA result
  p1 <- ridgeplot(pathway_res)
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
}


fig.geseaplot <- function(gsea_res, genesetid, by_scoring) 
{
  #' 
  #' @param gsea_res GSEA enrichment result
  #' @param genesetid which gene set shall be visualised
  #' @param by_scoring choose between all, runningScore, position
  p1 <- gseaplot(gsea_res, geneSetID = 1, title = gsea_res$Description[1], by = "all")
  cowplot::plot_grid(p1, ncol = 1)
}


fig.geseaplot2 <- function(gsea_res, show_genesetids) 
{
  #' color = c("#E495A5", "#86B875", "#7DB0DD"), title = gsea_res$Description[1]
  #' @param gsea_res GSEA enrichment result
  #' @param show_genesetid which gene sets or set shall be visualised
  p1 <- gseaplot2(gsea_res, geneSetID = 1:show_genesetids, pvalue_table = TRUE,
                  ES_geom = "dot", subplots = 1:3)
  cowplot::plot_grid(p1, ncol = 1)
}


fig.gesearank <- function(gsea_res, genesetid) 
{
  #'
  #' @param gsea_res GSEA enrichment result
  #' @param genesetid which enriched gene set shall be visualised
  gsearank(gsea_res, genesetid, title = gsea_res[genesetid, "Description"])
}


fig.gseatable <- function(df, ranked_gene_list, fgsea_res) 
{
  p1 <- plotGseaTable(df, ranked_gene_list, fgsea_res,  gseaParam = 0.5)
  cowplot::plot_grid(p1, ncol = 1)
}

fig.gseaenrichment <- function(df, ranked_gene_list) 
{
  p1 <- plotEnrichment(df, ranked_gene_list)
  cowplot::plot_grid(p1, ncol = 1)
}
