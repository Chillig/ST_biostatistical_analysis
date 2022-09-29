# Figure 8
rm(list = ls())
set.seed(42)  # important: as permutations and random initialization are used
options(java.parameters = "-Xmx4g")
library(fgsea)
library(ggplot2)
library('xlsx')

# 1. Load patient Bulk DEGs of:
# - Pso, L vs Pso NL
# - LP & LE, L vs LP & LE, NL
# - AE, L vs AE, NL
# 2. Rank List by log2FC or signed p.adj value
# 3. Gene Set: 
# - load DEGs from 3D skin models 
# - filter DEGs by up/down regulated genes
# - create df with column names up / down for each disease
# 4. Perform fgsea 

# 1. Create output directory
save_dir <- file.path(
  "/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/output/SFig7EFG",
  "Enrichmentplot_Keras_STDEGs",
  Sys.Date())
dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)

comparisons <- c('IL17A', 'IL13', 'IFNG')
for (val in seq(1, length(comparisons))) 
{
  # 1. Load Keras DEGs
  df.keras_degs <- readxl::read_excel(
    '/Users/christina.hillig/Documents/Projects/IGSSE-TUM_Projects/ST_ncISD/S_Eyerich_ncISD__ImmunePublication/Kera arrays recombinant cytokines__CH.xlsx', 
    sheet = comparisons[val], col_names = TRUE)
  # add singed p-value column
  df.keras_degs$signed.pval <- sign(df.keras_degs$foldchange) * df.keras_degs$pvalue
  # 2. Rank List
  df.ranked_degs <- df.keras_degs[order(df.keras_degs$signed.pval), ]
  # names vector
  ranked.genes <- df.ranked_degs$foldchange
  names(ranked.genes) <- df.ranked_degs$GenSymbol
  
  # 3. Gene Set
  df.degs <- readxl::read_excel(path = file.path(paste0('/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/output/reviewers/Whole_T_cell_matrix__cdr_project_patient_annotation_cyto/', comparisons[val], 
                                                        '/', comparisons[val], '_glmGamPoi_DGE.xlsx')), 
                                col_names = TRUE, sheet = 1)
  # 3.1 Get up/down regulated genes
  up.cytopos_spots <- df.degs[df.degs$log2fc > 1 & df.degs$padj < 0.05, ]$gene_symbol
  down.cytopos_spots <- df.degs[df.degs$log2fc < -1 & df.degs$padj < 0.05, ]$gene_symbol
  # 3.2 Create new named list 
  genesets <- list(up.cytopos_spots, down.cytopos_spots)
  names(genesets) <- c(paste0(comparisons[val], "_up"), paste0(comparisons[val], "_down"))
  
  # Enrichmentplots
  pdf(file = file.path(
    save_dir,
    paste('Enrichmentplot', paste0(comparisons[val], '__up', ".pdf"))), 
    width = 6, height =4)
  p.up <- plotEnrichment(genesets[[paste0(comparisons[val], "_up")]], ranked.genes) + 
    ggplot2::labs(title= paste0("Up-regulated genes in ",comparisons[val], "+ spots"))
  print(p.up)
  dev.off() 
  
  # pdf(file = file.path(
  #   save_dir,
  #   paste('Enrichmentplot', paste0(comparisons[val], '_down', ".pdf"))), 
  #   width = 6, height =4)
  # p.down <- plotEnrichment(genesets[[paste0(comparisons[val], "_down")]], ranked.genes) + 
  #   ggplot2::labs(title=paste0("Down-regulated genes in ",comparisons[val], "+ spots"))
  # print(p.down)
  # dev.off() 
}


