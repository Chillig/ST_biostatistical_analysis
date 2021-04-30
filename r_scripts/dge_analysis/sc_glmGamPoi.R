citation("glmGamPoi")
library(glmGamPoi)


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
