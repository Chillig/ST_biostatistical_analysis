#! /usr/bin/Rscript

# R version 4.0.3 Patched (2020-10-23 r79366)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7

# install Packages: BiocManager::install("reactome.db"), BiocManager::install("pathview"), install.packages("cowplot"), install.packages("xlsx")
# BiocManager::install("enrichplot"), install.packages("ggupset"), BiocManager::install("ReactomePA"), install.packages("qdapTools"), BiocManager::install("fgsea"), BiocManager::install("goseq"), BiocManager::install("clusterProfiler"), install.packages("Hmisc"), BiocManager::install("KEGG.db"), BiocManager::install("rWikiPathways")
########################################################################################
# GSEA analysis
# Gene Set Enrichment Analysis GSEA tests whether a set of genes of interest is enriched
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

# The analysis is performed by:
# I. ranking all genes in the data set
# II. identifying the rank positions of all members of the gene set in the ranked data set
# III. calculating an enrichment score (ES) that represents the difference between the observed rankings 
# and that which would be expected assuming a random rank distribution.

# Analysis included in that script:
# 1. GSEA analysis (Sergushichev 2016) --> get ES 
# 2. GO enrichment analysis  (Young et al. 2010)
# 3. KEGG pathway enrichment analysis (Yu et al. 2012)
########################################################################################
# .libPaths() : "/Users/christina.hillig/anaconda3/envs/py37_R4/lib/R/library"

# remove all variables in global environment
rm(list = ls())
library(docstring)

# libraries
# library("rWikiPathways")
# Gene set enrichment analysis:
# library(fgsea)

# GO-term Analysis:
# library(goseq)
library(clusterProfiler)

# for pathway enrichemnt analysis:
# library(DOSE) 
library(ReactomePA)

# Gene, GO-term and Pathway database:
library(org.Hs.eg.db) # human organism
# library(GO.db)
# library(reactome.db)

# Plot for pathway enrichemnt analysis:
# library(ggplot2)
# library(enrichplot) 
# library(cowplot)

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
load_files <- function(path_name_file)
{
  #' Read in .csv file 
  #' 
  #' Gets the full path and reads in csv file including the header
  #' @note The .csv file entries should be separated by ","
  #' @param path_name_file Path to .csv file
  #' @return Dataframe
  
  if (tail(strsplit(path_name_file, split = "[.]")[[1]], n = 1) == "xlsx") 
  {
    dge_df <- read.xlsx(path_name_file, sheetIndex = 1, header = TRUE)
  } else {
    dge_df <- read.csv(path_name_file, header = TRUE, sep = ";")
    
    if (ncol(dge_df)==1) {dge_df = read.csv(path_name_file, header = TRUE, sep = ",")}
  }
  
  return(dge_df)
}

insertRow <- function(existingDF, newrow, r) 
{
  existingDF[seq(r + 1, nrow(existingDF) + 1),] <- existingDF[seq(r, nrow(existingDF)),]
  existingDF[r,] <- newrow
}


get_gene_entrez <- function(dge_matrix)
{
  # get gene names
  gene_names <- row.names(dge_matrix)
  
  EG_IDs = mget(gene_names, revmap(org.Hs.egSYMBOL), ifnotfound = NA)
  dge_matrix <- dge_matrix[-which(is.na(EG_IDs)), ] 
  
  return(EG_IDs)
}

convert_gene_keggid <- function(dge_matrix)
{
  # get gene names
  gene_names <- row.names(dge_matrix)
  ##Get the Entrez gene IDs associated with those symbols
  # EG_IDs = unlist(mget(x=gene, envir=org.Hs.egALIAS2EG, ifnotfound=NA))
  
  EG_IDs = mget(gene_names, revmap(org.Hs.egSYMBOL), ifnotfound = NA)
  dge_matrix <- dge_matrix[-which(is.na(EG_IDs)), ] 
  EG_IDs <- EG_IDs[-which(is.na(EG_IDs))] 
  
  ##Then get the KEGG IDs associated with those entrez genes.
  KEGG_IDs = mget(as.character(EG_IDs), org.Hs.egPATH, ifnotfound = NA)
  
  
  ## Map KEGG_IDs to DEG matrix
  # change gene names to EG_IDs
  row.names(dge_matrix) <- unlist(EG_IDs, use.names = FALSE)
  dge_matrix <- dge_matrix[-which(is.na(KEGG_IDs)), ]
  KEGG_IDs <- KEGG_IDs[-which(is.na(KEGG_IDs))]
  
  # assign KEGG_IDs to EG_IDS in matrix
  # TODO how to treat rows with same same?
  counter <- 1
  for (i in seq_along(KEGG_IDs))
  {
    if (names(KEGG_IDs[i]) %in% row.names(dge_matrix)) {
      if (length(KEGG_IDs[i][[1]]) > 1)
      {
        index_id <- which(row.names(dge_matrix) == names(KEGG_IDs[i]))
        for (id_names in KEGG_IDs[i][[1]])
        {
          # copy and add row
          # check for dulicated row names
          if (id_names %in%  row.names(dge_matrix))
          {
            counter <- counter + 1
          }else{
            rownames(dge_matrix)[index_id] <- id_names
            dge_matrix[seq(index_id + 1, nrow(dge_matrix) + 1), ] <- dge_matrix[
              seq(index_id, nrow(dge_matrix)), ]
            dge_matrix[index_id, ] <- dge_matrix[index_id, ]
            counter <- counter + 1
            # insertRow(dge_matrix, dge_matrix[index_id, ], index_id)
          }
        }
      }else
      {
        index_id <- which(row.names(dge_matrix)  == names(KEGG_IDs[i]))
        # copy and add row
        if (KEGG_IDs[i][[1]] %in%  row.names(dge_matrix))
        {
          counter <- counter + 1
        }else
        {
          rownames(dge_matrix)[index_id] <- KEGG_IDs[i][[1]]
          counter <- counter + 1}
      }
    }
  }
  
  
  return(dge_matrix)
}


############## initialization phase
main = function(sample, date_file, replicate_type, dge_approach, minGSSize,
                lfc_factor,  fdr_value, p_value, pval_cut, correction_method) 
{
  print("-------------------------->  Pathway enrichment Analysis  <--------------------------")
  
  ############################### 1. Initialization phase
  sub_folders <- list.dirs(path = file.path(getwd(), "Input",
                                            paste(dge_approach, "output", sep = "_"),
                                            date_file, sample), 
                           full.names = TRUE, recursive = FALSE)
  
  # output paths
  dir.create(file.path(getwd(), "output", "Fig3C__PA_analysis"),  recursive = TRUE,
             showWarnings = FALSE)
  dir.create(file.path(getwd(), "Output", "Fig3C__PA_analysis", Sys.Date()),  recursive = TRUE,
             showWarnings = FALSE)
  save_dir <- file.path(getwd(), "Output", "Fig3C__PA_analysis", Sys.Date(), 
                        paste(dge_approach, "output", sep = "_"))
  dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)
  
  
  ######### START ANALYSIS
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
      # remove duplicated rows ..
      dge_list_results <- dge_list_results[!duplicated(dge_list_results$gene_symbol), ]
      rownames(dge_list_results) <- dge_list_results$gene_symbol
      dge_list_results <- dge_list_results[, c(colnames(dge_list_results) != "X"),
                                           drop = FALSE]
      # remove NA columns
      dge_list_results <- dge_list_results[, colSums(is.na(dge_list_results)) < 
                                             nrow(dge_list_results)]
      # remove NA rows
      dge_list_results <- na.omit(dge_list_results)
      
      # TDecide if cytokine shall be considered in analysis 
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
      # II. Identify the rank positions
      # barplot(ranked_genes)
      
      # III. Define significant DEx genes and background genes
      # Option 2: Determining the DE genes using edgeR 
      # 1.a) sort genes into groups belonging either to reference or test (control) condition
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
      
      # 1.b) get differentially expressed genes either up or down regulated 
      char_columns <- 2
      dge_list_results[ , char_columns] <- as.data.frame(sapply(dge_list_results$log2fc, 
                                                                as.numeric))
      de.common <- dge_list_results$entrezid[dge_list_results$pval < p_value 
                                             & !is.na(dge_list_results$pval) & 
                                               abs(dge_list_results$log2fc) > lfc_factor]
      de.common <- as.character(na.exclude(de.common))

      # b.) Background genes are all genes from our (sup-) data set
      bg_genes <- as.character(dge_list_results$entrezid)
      
      # transform entrez_id to factor
      dge_list_results$entrezid = as.factor(dge_list_results$entrezid)

      #############################################################################
      #################### ---> Pathway Enrichment Analysis <--- ##################
      #############################################################################
      # 1. ReactomePA Pathway enrichment analysis of a gene set: Database REACTOME
      # Note: Used to determine occuring protein receptors in dataset
      # Given vector of genes, function returns enriched pathways with FDR control.
      # 1.1 Find enriched Pathways for reference condition
      reactome_object.ref <- enrichPathway(gene = de.ref, # a vector of entrez gene id
                                           universe = bg_genes,
                                           organism = 'human', 
                                           qvalueCutoff = fdr_value, 
                                           pvalueCutoff = pval_cut, 
                                           pAdjustMethod = correction_method, 
                                           minGSSize = minGSSize,
                                           maxGSSize = 500,
                                           readable = T)
      
      # 1.2 Find enriched Pathways for control condition
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
      reactome.ctrl = setreadable_pa(paenrich_object=reactome_object.ctrl)
      reactome.ref = setreadable_pa(paenrich_object=reactome_object.ref) 
      
      
      #############################################################################
      ################### ---> Save results to csv file <--- ###################### 
      #############################################################################
      # Attention: 
      # ctrl (= negative log2FC) and ref (= positive log2FC) are switched for Immune publication
      enrichobject_to_df(paenrich_object=reactome.ctrl, condition='Cytoneg', 
                         pa_database='REACTOME', output_path=results_save_dir) 
      enrichobject_to_df(paenrich_object=reactome.ref, condition='Cytopos', 
                         pa_database='REACTOME', output_path=results_save_dir) 

      #############################################################################
      ############################### ---> PLOTS <--- ############################# 
      #############################################################################
      # Plot variables
      # select pathways or Enriched gene sets manually
      publication_pas = pathwaysofinterest()
      if (sample == 'single_cell') 
      {
        pas_publication = grep(
          paste('sc', signature_gene, sep="_"), keys(publication_pas), value=TRUE)
      } else 
      {
        pas_publication = grep(
          paste('st', signature_gene, sep="_"), keys(publication_pas), value=TRUE)
      }
      show_categories = publication_pas[[pas_publication]] #  3
      # show_categories = 4
      show_dotplot_categories = 15
      
      width_img = 16
      height_img = 8

      ######### ---> Save Pathway Enrichment Analysis Plots and Files <--- #########
      # 1. Reference Condition
      # If a gene is associated with two or more enriched PAs 
      # but less than those are shown than this results in a bug 
      # -> the log2fc of that gene is not correctly shown
      if (!is.null(nrow(reactome.ref)))
      {
        if (nrow(reactome.ref) > 1 & any(show_categories %in% reactome.ref$Description)) 
        {
          # Cnetplots to visualise enriched pathways
          pdf(file = file.path(results_save_dir,
                               "Cytopos_REACTOME_Pathway_Enrichment_Analysis.pdf"),
              width = width_img, height = height_img)
          print(fig.pathways.REACTOME(reactome_res = reactome.ref, 
                                      entrezid_log2fc = ranked_genes.ref,
                                      showCategories = show_categories))
          dev.off()
          
          # Dotplot to visualise enriched pathways
          pdf(file = file.path(results_save_dir, "Cytopos_REACTOME_dotplot.pdf"),
              width = height_img, height = height_img)
          print(fig.pathway.dotplot(pathway_res=reactome.ref,
                                    showCategories=show_dotplot_categories, 
                                    method='REACTOME'))
          dev.off()
          
        }
      }
      
      # # 2. Control Condition
      if (!is.null(nrow(reactome.ctrl)))
      {
        if (nrow(reactome.ctrl) > 1 & any(show_categories %in% reactome.ctrl$Description)) 
        {
          # Cnetplots to visualise enriched pathways
          pdf(file = file.path(results_save_dir,
                               "Cytoneg_REACTOME_Pathway_Enrichment_Analysis.pdf"),
              width = width_img, height = height_img)
          print(fig.pathways.REACTOME(reactome_res = reactome.ctrl, 
                                      entrezid_log2fc = ranked_genes.ctrl,
                                      showCategories = show_categories))
          dev.off()
          
          # Dotplot to visualise enriched pathways
          pdf(file = file.path(results_save_dir, "Cytoneg_REACTOME_dotplot.pdf"),
              width = height_img, height = height_img)
          print(fig.pathway.dotplot(pathway_res=reactome.ctrl,
                                    showCategories=show_dotplot_categories, 
                                    method='REACTOME'))
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

main(date_file = date_file, sample = sample, dge_approach = dge_approach, 
     replicate_type = replicate_type,
     lfc_factor = lfc_factor, fdr_value = fdr_value, p_value = p_value, pval_cut = pval_cut, 
     correction_method = test_method, minGSSize=minGSSize)
