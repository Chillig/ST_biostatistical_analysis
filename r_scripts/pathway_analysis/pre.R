#! /usr/bin/Rscript

# R version 4.0.3 Patched (2020-10-23 r79366)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7

# install Packages: BiocManager::install("reactome.db"), BiocManager::install("pathview"), install.packages("cowplot"), install.packages("xlsx")
# BiocManager::install("enrichplot"), install.packages("ggupset"), BiocManager::install("ReactomePA"), install.packages("qdapTools"), BiocManager::install("fgsea"), BiocManager::install("goseq"), BiocManager::install("clusterProfiler"), install.packages("Hmisc")
########################################################################################
# Analysis included in that script:
# ReactomePA pathway enrichment analysis 

# The steps inclcude:
# I. Replace gene symbol by entrezid of gene 
# II. Define significant DEx genes, rank them by their log2FC and read out background genes
# III. Perform Pathway enrichment analysis
# IV. Visualise Pathways in Dotplots and as Networks

########################################################################################
# .libPaths() : "/Users/christina.hillig/anaconda3/envs/py37_R4/lib/R/library"

# remove all variables in global environment
rm(list = ls())


# libraries
library(docstring)

# Gene and Pathway database:
library(org.Hs.eg.db) # human organism

# variable or objects functions:
library(dplyr) 


#################################### -> functions <- ##################################
rename_genetoentrezid <- function(df.dge_results) 
{
  # I. replace gene symbol by entrezid of gene
  if ("entrezid" %nin% colnames(df.dge_results))
  {
    # works approach 
    Entrez_ID <- mapIds(org.Hs.eg.db, row.names(df.dge_results), 'ENTREZID', 'SYMBOL')
    Entrez_ID_df = list2df(Entrez_ID, col1 = "entrezid", col2 = "gene_symbol")
    Entrez_ID_df = Entrez_ID_df[!is.na(Entrez_ID_df$entrezid), ]
    
    gseaDat <- filter(df.dge_results, !is.na(Entrez_ID))
    df.dge_results <- filter(df.dge_results, !is.na(Entrez_ID))
    df.dge_results$entrezid <- Entrez_ID_df$entrezid
  }
  
  return(df.dge_results)
}

do_rank_genes <- function(df.dge_results) 
{
  # II. Rank all genes based on their fold change
  ranked_genes <- df.dge_results$log2fc
  names(ranked_genes) <- df.dge_results$gene_symbol
  df.dge_results <- sort(df.dge_results, decreasing = T)
  
  return(ranked_genes)
}


get_significantgenes <- function(df.dge_results, p_value, lfc_factor, op) 
{
  # II. Get significantly differentially expressed genes
  df.sig <- df.dge_results[df.dge_results$pval < p_value & 
                             !is.na(df.dge_results$pval) & 
                             op(df.dge_results$log2fc, lfc_factor), ]
  df.sig <- df.sig$entrezid
  df.sig <- as.character(na.exclude(df.sig))       
  
  return(df.sig)
}


