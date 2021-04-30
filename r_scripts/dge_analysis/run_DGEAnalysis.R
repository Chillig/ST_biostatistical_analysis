source('../scripts/data_preparation.R')
source('../scripts/sc_glmGamPoi.R')

run_dge.analysis <- function(matrix, metadata, num_samples_donor, exp_design.function,
                             multitest.method, dge_method, fdr.cut) 
{
  matrix.metadata.design = prepare_data_dgeanalysis(
    matrix = matrix, metadata = metadata, num_samples_donor  = num_samples_donor,
    design.function = exp_design.function)
  
  expr.matrix <-  matrix.metadata.design[[1]]
  meta.Data <- matrix.metadata.design[[2]]
  design.function <- matrix.metadata.design[[3]]
  size.factor <- metadata$Sizefactor
  
  message("Choosen method: glmGamPoi")
  dge_result <- run_glmgampoi(count.matrix = expr.matrix, 
                              design = design.function, 
                              metaData = meta.Data, sizefactor=size.factor,
                              multitest.method = multitest.method) 
  
  #   sort summary by adjusted p-val
  dge_result <- dge_result[order(dge_result$padj), ]
  
  #   remove NA's from result
  dge_result <- na.omit(dge_result)
  
  return(list(dge_result, design.function, meta.Data))
}
