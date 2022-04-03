#! /usr/bin/Rscript

# R version 4.0.3 Patched (2020-10-23 r79366)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7

# install Packages: BiocManager::install("reactome.db"), BiocManager::install("pathview"), install.packages("cowplot"), install.packages("xlsx")
# BiocManager::install("enrichplot"), install.packages("ggupset"), BiocManager::install("ReactomePA"), install.packages("qdapTools"), BiocManager::install("fgsea"), BiocManager::install("goseq"), BiocManager::install("clusterProfiler"), install.packages("Hmisc")
########################################################################################
# Analysis included in that script:
# ReactomePA pathway enrichment analysis 

# The steps inlcude:
# I. Replace gene symbol by entrezid of gene
# II. Define significant DEx genes and background genes
# III. Perform Pathway enrichement analysis
# IV. Visualise Pathways in Dotplots and as Networks

########################################################################################
# .libPaths() : "/Users/christina.hillig/anaconda3/envs/py37_R4/lib/R/library"

# remove all variables in global environment
rm(list = ls())


# libraries
library(docstring)

# GO-term Analysis:
library(clusterProfiler)

# for pathway enrichemnt analysis:
library(ReactomePA)

# Gene and Pathway database:
library(org.Hs.eg.db) # human organism
# library(GO.db)
# library(reactome.db)

# variable or objects functions:
library(dplyr) 
# library(Hmisc)
library(pathview)
library(openxlsx)
library("xlsx")
library("qdapTools")

# Plot functions
source("PA_plotting_functions.R")
source("ImmunePublication_pathways.R")
source('utils.R')


#################################### -> functions <- ##################################
#' Read in .csv of .xlsx file 
#' 
#' Gets the full path and reads in csv file including the header
#' @note The .csv file entries should be separated by ","
#' @param path_name_file Path to .csv file
#' @return Dataframe
load_files <- function(path_name_file)
{
  if (tail(strsplit(path_name_file, split = "[.]")[[1]], n = 1) == "xlsx") 
  {
    dge_df <- xlsx::read.xlsx(path_name_file, sheetIndex = 1, header = TRUE)
  } else {
    dge_df <- xlsx::read.csv(path_name_file, header = TRUE, sep = ";")
    
    if (ncol(dge_df) == 1) {dge_df = read.csv(path_name_file, header = TRUE, sep = ",")}
  }
  
  return(dge_df)
}


main = function(sample, date_file, replicate_type, dge_approach, minGSSize,
                lfc_factor,  fdr_value, p_value, pval_cut, correction_method, plot_signaturecytokine)
{
  print("-------------------------->  Pathway enrichment Analysis  <--------------------------")
  
  ############################### 1. Initialization phase
  sub_folders <- list.dirs(
    path = file.path(getwd(), "Input", paste(dge_approach, "output", sep = "_"),
                     date_file, sample), full.names = TRUE, recursive = FALSE)
  
  # output paths
  dir.create(file.path(getwd(), "output", "Fig3C__PA_analysis"),  recursive = TRUE,
             showWarnings = FALSE)
  dir.create(file.path(getwd(), "Output", "Fig3C__PA_analysis", Sys.Date()),  recursive = TRUE,
             showWarnings = FALSE)
  save_dir <- file.path(getwd(), "Output", "Fig3C__PA_analysis", Sys.Date(), 
                        paste(dge_approach, "output", sep = "_"))
  dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)
  
  
  ############################### 2. START ANALYSIS
  for (sub_folder in sub_folders) 
  {
    # get all DGE .csv files in subfolders
    all_filenames = list.files(path = sub_folder, pattern = c("*.csv|*.xlsx"), 
                               recursive = TRUE)
    for (dge_filename in all_filenames[!grepl("metaData*", all_filenames)])
    {
      print(dge_filename)
      
      # output PA Analysis
      name_comparison = tail(strsplit(dge_filename, split = "[0-9]\\__", 
                                      perl = TRUE)[[1]], n = 1)
      name_comparison <- strsplit(name_comparison, "_DGE|_DEG")[[1]][1]
      results_save_dir <- file.path(save_dir, sample,
                                    tail(strsplit(sub_folder, .Platform$file.sep)[[1]], 
                                         n = 1), 
                                    strsplit(strsplit(dge_filename, .Platform$file.sep)[[1]], 
                                             "[.]")[[1]][1])
      dir.create(results_save_dir, recursive = TRUE, showWarnings = FALSE)
      
      # path to input csv files
      dge_results_load_dir <- file.path(sub_folder,  dge_filename)
      
      # Load dge list with gene names and p-value (and/or log2FC)
      dge_list_results <- load_files(path_name_file = dge_results_load_dir)
      
      # Remove redundant inforamtion:
      # # -remove duplicated rows
      dge_list_results <- dge_list_results[!duplicated(dge_list_results$gene_symbol), ]
      rownames(dge_list_results) <- dge_list_results$gene_symbol
      dge_list_results <- dge_list_results[, c(colnames(dge_list_results) != "X"),
                                           drop = FALSE]
      # # -remove NA columns
      dge_list_results <- dge_list_results[, colSums(is.na(dge_list_results)) < 
                                             nrow(dge_list_results)]
      # # -remove NA rows
      dge_list_results <- na.omit(dge_list_results)
      
      # If cytokine shall be considered in analysis 
      if (plot_signaturecytokine == TRUE) 
      {
        signature_gene = strsplit(name_comparison, .Platform$file.sep)[[1]][1]
      } else 
      {
        # remove signature cytokine for plots
        dge_list_results = dge_list_results[dge_list_results$gene_symbol != signature_gene, ]
      }
      
      ########################### ---> Preparation <--- ##########################
      # I. replace gene symbol by entrezid of gene
      if ("entrezid" %nin% colnames(dge_list_results))
      {
        # works approach 
        Entrez_ID <- mapIds(org.Hs.eg.db, row.names(dge_list_results), 'ENTREZID', 'SYMBOL')
        Entrez_ID_df = list2df(Entrez_ID, col1 = "entrezid", col2 = "gene_symbol")
        Entrez_ID_df = Entrez_ID_df[!is.na(Entrez_ID_df$entrezid), ]
        
        gseaDat <- filter(dge_list_results, !is.na(Entrez_ID))
        dge_list_results <- filter(dge_list_results, !is.na(Entrez_ID))
        dge_list_results$entrezid <- Entrez_ID_df$entrezid
      }
      
      # I. Rank all genes based on their fold change
      # Option 2: Determining the DE genes using edgeR 
      if ('log2fc' %in% colnames(dge_list_results)) 
      {
        ranks <- dge_list_results$log2fc 
      } else 
      {
        ranks <- dge_list_results$log2FC 
        names(dge_list_results)[names(dge_list_results) == 'log2FC'] <- 'log2fc'
      }
      
      names(ranks) <- dge_list_results$gene_symbol #dge_list_results$entrezid
      ranked_genes <- sort(ranks, decreasing = T)
      
      # II. Define significant DEx genes and background genes
      # Option 2: Determining the DE genes using edgeR 
      # II.a) sort genes into groups belonging either to reference or test (control) condition
      df.ref <- dge_list_results[
        dge_list_results$pval < p_value & !is.na(dge_list_results$pval) & 
          dge_list_results$log2fc < -lfc_factor, ]
      de.ref <- df.ref$entrezid
      de.ref <- as.character(na.exclude(de.ref))        
      ranked_genes.ref <- df.ref$log2fc
      names(ranked_genes.ref) <- df.ref$gene_symbol
      
      df.ctrl <- dge_list_results[
        dge_list_results$pval < p_value & !is.na(dge_list_results$pval) & 
          dge_list_results$log2fc > lfc_factor, ]
      de.ctrl <- df.ctrl$entrezid
      de.ctrl <- as.character(na.exclude(de.ctrl))
      ranked_genes.ctrl <- df.ctrl$log2fc
      names(ranked_genes.ctrl) <- df.ctrl$gene_symbol
      
      # II.b) get differentially expressed genes either up or down regulated 
      char_columns <- 2
      dge_list_results[ , char_columns] <- as.data.frame(sapply(dge_list_results$log2fc, 
                                                                as.numeric))
      de.common <- dge_list_results$entrezid[dge_list_results$pval < p_value 
                                             & !is.na(dge_list_results$pval) & 
                                               abs(dge_list_results$log2fc) > lfc_factor]
      de.common <- as.character(na.exclude(de.common))

      # II.c) Background genes are all genes from our (sup-) data set
      bg_genes <- as.character(dge_list_results$entrezid)
      
      # transform entrez_id to factor
      dge_list_results$entrezid = as.factor(dge_list_results$entrezid)

      #############################################################################
      #################### ---> Pathway Enrichment Analysis <--- ##################
      #############################################################################
      # III. ReactomePA Pathway enrichment analysis of a gene set: Database REACTOME
      # Note: Used to determine occuring protein receptors in dataset
      # Given vector of genes, function returns enriched pathways with FDR control.
      # III.a) Find enriched Pathways for reference condition
      reactome_object.ref <- enrichPathway(gene = de.ref, # a vector of entrez gene id
                                           universe = bg_genes,
                                           organism = 'human', 
                                           qvalueCutoff = fdr_value, 
                                           pvalueCutoff = pval_cut, 
                                           pAdjustMethod = correction_method, 
                                           minGSSize = minGSSize,
                                           maxGSSize = 500,
                                           readable = T)
      
      # III.b) Find enriched Pathways for control condition
      reactome_object.ctrl <- enrichPathway(gene = de.ctrl, # a vector of entrez gene id
                                            universe = bg_genes,
                                            organism = 'human', 
                                            qvalueCutoff = fdr_value, 
                                            pvalueCutoff = pval_cut, 
                                            pAdjustMethod = correction_method, 
                                            minGSSize = minGSSize,
                                            maxGSSize = 500,
                                            readable = T)

      
      ################### ---> convert gene ID to Symbol <--- ################### 
      # Pathway Enrichment
      reactome.ctrl = setreadable_pa(paenrich_object = reactome_object.ctrl)
      reactome.ref = setreadable_pa(paenrich_object = reactome_object.ref) 
      
      
      #############################################################################
      ################### ---> Save results to csv file <--- ###################### 
      #############################################################################
      # Attention: 
      # ctrl (= negative log2FC) and ref (= positive log2FC) are switched for Immune publication
      enrichobject_to_df(paenrich_object = reactome.ctrl, condition = 'Cytoneg', 
                         pa_database = 'REACTOME', output_path = results_save_dir) 
      enrichobject_to_df(paenrich_object = reactome.ref, condition = 'Cytopos', 
                         pa_database = 'REACTOME', output_path = results_save_dir) 

      #############################################################################
      ############################### ---> PLOTS <--- ############################# 
      #############################################################################
      # IV.a) Plot variables
      # select pathways or Enriched gene sets manually
      publication_pas = pathwaysofinterest()
      if (sample == 'single_cell') 
      {
        pas_publication = grep(
          paste('sc', signature_gene, sep = "_"), keys(publication_pas), value = TRUE)
      } else 
      {
        pas_publication = grep(
          paste('st', signature_gene, sep = "_"), keys(publication_pas), value = TRUE)
      }
      show_categories = publication_pas[[pas_publication]] #  3
      show_dotplot_categories = 15
      
      width_img = 16
      height_img = 8

      ######### ---> Save Pathway Enrichment Analysis Plots and Files <--- #########
      # IV.b) Reference Condition
      # If a gene is associated with two or more enriched PAs 
      # but less than those are shown than this results in a bug 
      # -> the log2fc of that gene is not correctly shown
      if (!is.null(nrow(reactome.ref)))
      {
        if (nrow(reactome.ref) > 1 & any(show_categories %in% reactome.ref$Description)) 
        {
          # Cnetplots to visualise enriched pathways
          pdf(file = file.path(results_save_dir, "Cytopos_REACTOME_Pathway_Enrichment_Analysis.pdf"),
              width = width_img, height = height_img)
          print(fig.pathways.REACTOME(reactome_res = reactome.ref, 
                                      entrezid_log2fc = ranked_genes.ref,
                                      showCategories = show_categories))
          dev.off()
          
          # Dotplot to visualise enriched pathways
          pdf(file = file.path(results_save_dir, "Cytopos_REACTOME_dotplot.pdf"),
              width = height_img, height = height_img)
          print(fig.pathway.dotplot(pathway_res = reactome.ref,
                                    showCategories = show_dotplot_categories, 
                                    method = 'REACTOME'))
          dev.off()
          
        }
      }
      
      # IV.c) Control Condition
      if (!is.null(nrow(reactome.ctrl)))
      {
        if (nrow(reactome.ctrl) > 1 & any(show_categories %in% reactome.ctrl$Description)) 
        {
          # Cnetplots to visualise enriched pathways
          pdf(file = file.path(results_save_dir, "Cytoneg_REACTOME_Pathway_Enrichment_Analysis.pdf"),
              width = width_img, height = height_img)
          print(fig.pathways.REACTOME(reactome_res = reactome.ctrl, 
                                      entrezid_log2fc = ranked_genes.ctrl,
                                      showCategories = show_categories))
          dev.off()
          
          # Dotplot to visualise enriched pathways
          pdf(file = file.path(results_save_dir, "Cytoneg_REACTOME_dotplot.pdf"),
              width = height_img, height = height_img)
          print(fig.pathway.dotplot(pathway_res = reactome.ctrl,
                                    showCategories = show_dotplot_categories, 
                                    method = 'REACTOME'))
          dev.off()
        }
      }
    }
  }
  
  sessionInfo()
}

date_file <- "2021-02-01"
replicate_type = "biological"
sample = "spatial"
dge_approach = "T-cell_countmatrix"

lfc_factor = 1
fdr_value = 0.05
pval_cut = 0.05
p_value = 0.05
minGSSize = 10

# FDR and BH are more conservative than bonferroni
test_method =  "BH" 

plot_cytokine = TRUE

# TODO: 
# 1. save used params in excel file
# 2. save read in DEG files used in analysis in a .txt file


main(date_file = date_file, sample = sample, dge_approach = dge_approach, 
     replicate_type = replicate_type,
     lfc_factor = lfc_factor, fdr_value = fdr_value, p_value = p_value, pval_cut = pval_cut, 
     correction_method = test_method, minGSSize=minGSSize, plot_signaturecytokine=plot_cytokine)
