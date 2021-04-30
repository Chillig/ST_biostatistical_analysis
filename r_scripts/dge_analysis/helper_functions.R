# ------------------------- Helpfer functions -------------------------
delete_zero_conditions <- function(countmatrix)
{
  #' @description  Remove all conditions which have only zero entries
  #'
  #' @param countmatrix
  #' @return cleaned countmatrix
  
  message('Remove conditions with only zeros')
  
  # get zero columns
  zero_columns <- colSums(countmatrix) [colSums(countmatrix) == 0]
  conditions <- colnames(countmatrix)
  
  vec_cond = vector()
  sample_name_short = vector()
  for (cond in 1:length(conditions))
  {
    temp_condname <- strsplit(conditions[cond], "[.]")[[1]]
    vec_cond[cond] <- temp_condname[3]
    sample_name_short[cond] <- temp_condname[1]
  }
  
  names_conditions <- names(zero_columns)
  if (length(names_conditions) > 0)
  {
    vec_condition <- vector()
    for (cond in 1:length(names_conditions))
    {
      vec_condition[cond] <- strsplit(names_conditions[cond], "[.]")[[1]][2]
    }
    
    num_samples <- length(unique(sample_name_short))
    occurence_conditions <- table(vec_condition)
    
    index_zero_conditions <- vector()
    for (o_c in 1:length(occurence_conditions))
    {
      if (occurence_conditions[o_c] == num_samples)
      {
        index_zero_conditions <- c(index_zero_conditions, occurence_conditions[o_c])
      }
    }
    
    # remove zero conditions
    x <- names(index_zero_conditions)
    for (i_c in 1:length(x))
    {
      countmatrix <- countmatrix[, !stringr::str_detect(vec_cond, x[i_c])]
      vec_cond <- vec_cond[!stringr::str_detect(vec_cond, x[i_c])]
    }
  }
  
  return(countmatrix)
}
