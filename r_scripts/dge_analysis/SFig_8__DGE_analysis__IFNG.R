
library(rhdf5)
library(glmGamPoi)
library('tibble')
library("org.Hs.eg.db")


#' Perform DGE Analysis using Bioconductor package 'glmGamPoi' designed for single-cell data sets
#'
#' @family DGE Analysis functions
#' @title DGE Analysis using glmGamPoi
#' @description Calculate DGE genes from single-cell or ST data set. 
#' @param countmatrix counts
#' @param design_function design function e.g.: ~ condition
#' @param metaData metaData containing additional information of each cell or spot
#' @param sizefactor size factors calculated on the whole data set
#' @return data frame containing the DGE Analysis result
#' @references
#' @export
run_glmgampoi <- function(count.matrix, design, metaData, multitest.method, 
                          sizefactor="deconvolution") 
{
  message("glmGamPoi")
  
  # Create a 'glmGamPoi' object
  gampoi_dds = glmGamPoi::glm_gp(as.matrix(count.matrix), design = design, on_disk = FALSE,
                                 size_factors = sizefactor)
  
  # Test for Differential Expression
  # Conduct a quasi-likelihood ratio test for a Gamma-Poisson fit. 
  # set: 
  # contrast = 'conditionreference'
  # -> find differentially expressed genes between conditions
  # pseudobulk_by = rownames(metaData) 
  # -> assuming that each spot is an individual ST samples as the spots 
  # are separated by distance of 100 um
  
  result_gampoi= glmGamPoi::test_de(gampoi_dds, 
                                    contrast = tail(colnames(gampoi_dds$model_matrix), n=1), 
                                    pseudobulk_by = rownames(metaData),
                                    full_design = gampoi_dds$model_matrix,
                                    pval_adjust_method = multitest.method )
  
  # Store result in data frame with pre-determined columns
  df = data.frame(gene_symbol = result_gampoi$name,
                  pval = result_gampoi$pval,
                  padj = result_gampoi$adj_pval,
                  log2fc = result_gampoi$lfc)
  
  return(df)
}


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


# ------------------------- Prepare Data ------------------------- 
prepare_data_dgeanalysis <- function(matrix, metadata, design.function) 
{
  # check if columns of Count Matrix are equal to rows of MetaData
  bool_gene_matrix_metadata <- all(rownames(metadata) %in% 
                                     colnames(matrix))
  if (bool_gene_matrix_metadata == FALSE)
  {
    matrix <- matrix[, rownames(metadata)]
    print(all(rownames(metadata) == colnames(matrix)))
  }
  
  # 3. create design matrix
  # add cellular detection rate (cdr)
  cdr <- scale(colMeans(matrix > 0))
  metadata['cdr'] = cdr
  design.matrix = model.matrix(design.function, metadata)
  
  # convert metaData column to factors
  metadata['patient'] <- lapply(metadata['patient'], as.factor)
  metadata['DISEASE'] <- lapply(metadata['DISEASE'], as.factor)
  metadata['biopsy_type'] <- lapply(metadata['biopsy_type'], as.factor)
  metadata['sample'] <- lapply(metadata['sample'], as.factor)
  metadata['specimen'] <- lapply(metadata['specimen'], as.factor)
  metadata['tissue_layer'] <- lapply(metadata['tissue_layer'], as.factor)
  
  return(list(matrix, metadata, design.matrix))
}


adata.dir <- '/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/output/reviewers/Figure_S8/2022-09-16/SDC_adata.h5'

# Read out observables
adata.obs <- rhdf5::h5read(adata.dir, "/obs/")
adata.vars <- rhdf5::h5read(adata.dir, "/var/")
# load log counts
# adata.X <- h5read(adata.dir, "/X/data")
# adata.X <- matrix(unlist(adata.X), nrow = 16685, ncol = 8526)
adata.counts <-rhdf5:: h5read(adata.dir, "/layers/counts")
row.names(adata.counts) <- adata.vars$gene_name
# make barcode names unique
unique.barcodes <- make.names(adata.obs$`_index`, unique = TRUE)
colnames(adata.counts) <- unique.barcodes

# Read out metaData with sizeFactors, patient, project, tissue layers
df.metaData <- as.data.frame(adata.obs[c('sample', 'project', 'patient', 'DISEASE', 'biopsy_type', 
                           'tissue_layer', 'size_factors', 'specimen', 'IFNG_in_sdcc_r4')])
row.names(df.metaData) <- unique.barcodes
# Design
design.func = ~cdr + project + patient + tissue_layer + IFNG_in_sdcc_r4

# Remove nn spots label = 2
mask.nn <- df.metaData$IFNG_in_sdcc_r4 == 2
df.metaData <- df.metaData[!mask.nn, ]
adata.counts <- adata.counts[ , !mask.nn]


# Prepare count matrix and metaData
matrix_metadata_design.matrix <- prepare_data_dgeanalysis(matrix=adata.counts, metadata=df.metaData,
                                                          design.function=design.func) 
df.metaData <- matrix_metadata_design.matrix[[2]]

# Run DGE analysis
df.dge <- run_glmgampoi(count.matrix=matrix_metadata_design.matrix[[1]], 
                        design=matrix_metadata_design.matrix[[3]],
                        metaData=df.metaData, multitest.method='BH', sizefactor=df.metaData$size_factors) 

# add entrezID and description
entrezid.gene.name <- map_genesymbol_entrezid_genename(
  gene.symbols=df.dge$gene_symbol) 
df.dge <- do_add.column(df=df.dge, entrezid=entrezid.gene.name[[1]], gene.name=entrezid.gene.name[[2]]) 

print(df.dge[df.dge$gene_symbol == 'IFNG',  ])

write.csv(df.dge, file=file.path("/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/output/reviewers/Figure_S8/glmGamPoi_DEGs", '2022-09-16',
                                 "IFNG", "DEGs_IFNG_optimal_radius4.csv"))
write.xlsx(df.dge, file=file.path("/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/output/reviewers/Figure_S8/glmGamPoi_DEGs", '2022-09-16',
                                  "IFNG", "DEGs_IFNG_optimal_radius4.xlsx"))

df.dge[order(df.dge$padj),]


# References
citation("glmGamPoi")
