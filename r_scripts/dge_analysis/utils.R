# Helper / utils Functions

# variable or objects functions:
library(dplyr)    # alternatively, this also loads %>%

# database
library(org.Hs.eg.db) 

# source required R-script
source("../scripts/gene_lists.R")

# ------------------------- Add information ------------------------- 
map_genesymbol_entrezid_genename <- function(gene.symbols) 
{
  #'
  #' @description 
  #' @param gene.symbols
  #' @return entrezid and gene name of genes
  
  genename <- vector()
  entrezid <- vector()
  for (sym_c in 1:length(gene.symbols)) 
  {
    # Map gene symbols to entrezid
    entrezid[sym_c] <- as.character(mget(gene.symbols[sym_c], org.Hs.eg.db::org.Hs.egSYMBOL2EG,
                                         ifnotfound = NA))
    # Map entrezid to gene name
    genename[sym_c] <- as.character(mget(entrezid[sym_c], org.Hs.eg.db::org.Hs.egGENENAME,
                                         ifnotfound = NA))
  }

  
  return(list(entrezid, genename))
}


do_add.column <- function(df, entrezid, gene.name) 
{
  # add genenames and entrezid to dataframe
  df <- tibble::add_column(df, entrezid = entrezid, .after = "gene_symbol")
  df <- tibble::add_column(df, gene_name = gene.name, .after = "gene_symbol")
  
  return(df)
}


add.hkg <- function(dge.df.object) 
{
  # house keeping genes
  hkg <- get_hkg(genes = dge.df.object$gene_symbol)
  house_keeping_genes <- hkg[[1]]
  if (length(hkg[[2]]) > 0) {
    index_hkg <- match(house_keeping_genes, dge.df.object$gene_symbol)
    # remove NA/s
    house_keeping_genes <- house_keeping_genes[!index_hkg %in% NA]
    index_hkg <- index_hkg[!index_hkg %in% NA]
    
    df_degres_hkg <- dge.df.object[index_hkg, ]
    df_degres_hkg$hkg <- rep("y", length.out = length(rownames(df_degres_hkg)))
    
    # index of non hkg
    non_hkg_index <- rownames(
      dge.df.object[!dge.df.object$gene_symbol %in% house_keeping_genes, ])
    
    df_degres_non_hkg <- dge.df.object[!dge.df.object$gene_symbol %in% house_keeping_genes, ]
    df_degres_non_hkg$hkg <- rep("n", length.out = length(rownames(df_degres_non_hkg)))
    
    dge.df.object <- rbind(df_degres_non_hkg, df_degres_hkg)
    
    
  } else 
  {
    dge.df.object$hkg <- rep("n", length.out = length(rownames(dge.df.object)))
    }
  
  
  return(dge.df.object)
}


# ------------------------- Up and Down regulated Data frames------------------------- 
save_updown.df <- function(df.res, pval.cut, lfc.cut, output.path, dge.method) 
{
  #' @description  Save up and down regulated genes in csv and html file
  #' @param df.res: dataframe containing DGE Analysis 
  #' @param pval.cut: p-value threshold
  #' @param lfc.cut: Log2FC threshold
  #' @param output.path: path wehre to save the files
  #' @param dge.method: which DGE Analysis method has been used
  
  try(
    {
      # Get up and down regulated genes
      up_deg <- df.res[df.res$log2fc > lfc.cut & df.res$pval < pval.cut,]
      up_deg$DEx <- rep("up", length.out = length(rownames(up_deg)))
      down_deg <- df.res[df.res$log2fc < -lfc.cut & df.res$pval < pval.cut,]
      down_deg$DEx <- rep("down", length.out = length(rownames(down_deg)))
      combined_up_down_genes <- rbind(up_deg, down_deg)
      
      # highlight FDR < 0.05 cells
      dt_combined_genes <- datatable(combined_up_down_genes, escape = FALSE, 
                                     editable = list(target = "all")) %>%
        formatStyle(columns = "padj", 
                    background = styleInterval(c(0., 0.05), c("white", "red", "white"))) 
      
      # save up and down regulated genes + their results form DGE analysis in .csv file 
      write.table(combined_up_down_genes, 
                  file.path(output.path, 
                            paste(sample, str_comparison,"DEG", dge.method, "R.csv", 
                                  sep = "_")), 
                  row.names = FALSE, sep = ",", 
                  col.names = colnames(combined_up_down_genes), append = FALSE) 
      
      # save data table highlight padj values below 0.05
      DT::saveWidget(dt_combined_genes, 
                     file.path(output.path, 
                               paste("w_replicate", sample, 
                                     str_comparison, "DEG", dge_method, "R.html",  sep = "_")))
    }, silent = TRUE)
}


get_df.mean_updown_hkg <- function(df, log2fc.cut, pval.cut, fdr.cut, cutby)
{
  #' @description Define the up/down regulated genes 
  #' @param df: data frame containing results of DGE Analysis
  #' @param log2fc.cut: log2FC threshold
  #' @param pval.cut: p-value threshold
  #' @param fdr.cut: FDR / p-adjusted value threshold
  #' @param cutby: Apply  threshold such that DEx genes are determined by FDR or p-value
  
  df_dge.res <- df[order(df$padj), ]
  df.log2fc <- df_dge.res$log2fc
  df.padj <- df_dge.res$padj
  df.pval <- df_dge.res$pval
  
  # rename first column if it is not "gene_symbol"
  if (is.null(df_dge.res$gene_symbol)) {
    df_dge.res <- rename(df_dge.res, c("gene_symbol" = colnames(df_dge.res)[1])) 
  } 
  
  # Check if a FDR cut or a p-value threshold shall be applied
  if (cutby == 'fdr') 
  {
    # Up
    up_deg <- df_dge.res[df.log2fc > log2fc.cut & df.padj < fdr.cut, ]
    up_deg$DEx <- rep("up", length.out = length(rownames(up_deg)))
    
    # down
    down_deg <- df_dge.res[df.log2fc < (-log2fc.cut) & df.padj < fdr.cut, ]
    down_deg$DEx <- rep("down", length.out = length(rownames(down_deg)))
    
    # intermediate
    inter_deg <- df_dge.res[!(abs(df.log2fc) > log2fc.cut & df.padj < fdr.cut), ]
    inter_deg$DEx <- rep("intermediate", length.out = length(rownames(inter_deg)))
  } else {
    # Up
    up_deg <- df_dge.res[df.log2fc > log2fc.cut & df.pval < pval.cut, ]
    up_deg$DEx <- rep("up", length.out = length(rownames(up_deg)))
    
    # down
    down_deg <- df_dge.res[df.log2fc < (-log2fc.cut) & df.pval < pval.cut, ]
    down_deg$DEx <- rep("down", length.out = length(rownames(down_deg)))
    
    # intermediate
    inter_deg <- df_dge.res[!(abs(df.log2fc) > log2fc.cut & df.pval < pval.cut), ]
    inter_deg$DEx <- rep("intermediate", length.out = length(rownames(inter_deg)))
  }

  # Merge data frames with up, down and intermediate genes 
  combined_updowninter_genes <- rbind(up_deg, down_deg, inter_deg)
  
  return(combined_updowninter_genes)
}
