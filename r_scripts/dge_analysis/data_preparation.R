
prepare_data <- function(countmatrix, condmetadata)
{
  #'
  #' @param countmatrix
  #' @param condmetadata
  #' @return pairwise_countmatrix, pairwise_metadata, cell_name
  
  cells <- colnames(countmatrix)
  vec_cond = vector()
  sample_name_short = vector()
  condition_names = vector()
  for (cond in 1:length(cells))
  {
    temp_conditionnames <- strsplit(cells[cond], "[.]")[[1]]
    condition <- tail(temp_conditionnames, n = 1)
    if ("_" %in% condition) {
      vec_cond[cond] <- strsplit(condition, "_")[[1]][2]
    } else {
      vec_cond[cond] <- condition
    }
    sample_name_short[cond] <- temp_conditionnames[1] 
    condition_names[cond] <- tail(temp_conditionnames, n = 2)[1]
  }
  # get number of conditions
  num_levels <- length(unique(vec_cond))
  name_cond <- unique(vec_cond)
  # check if values 1 and 2 in conditions: 
  # produces Warning message: NAs introduced by coercion 
  if (is.na(as.integer(name_cond[1]))) {
    vec_cond <- as.integer(x = factor(x = vec_cond))
    
    # rename rownames of count matrix
    colcountmatrix <- vector()
    for ( cond in 1:length(cells)) 
    {
      colcountmatrix[cond] <- paste(paste(head(strsplit(cells[cond], "[.]")[[1]], -1), 
                                          collapse = "."), 
                                    vec_cond[cond], sep = "_")
    }
    colnames(countmatrix) <- colcountmatrix
  }
  
  # get combinations of comparisons
  vec_comb <- gtools::combinations(num_levels, r = 2, v = seq(1, num_levels, by = 1))
  
  # Create count_matrix for pair-wise comparison
  pairwise_countmatrix <- list()
  pairwise_metadata <- list()
  for (i_c in 1:nrow(vec_comb))
  {
    bool_index_conditions <- vec_cond %in% vec_comb[i_c, ]
    pairwise_countmatrix[[i_c]] <- countmatrix[ , bool_index_conditions]
    gene_names <- row.names(countmatrix)
    rownames(pairwise_countmatrix[[i_c]]) <- make.names(gene_names, unique = TRUE)
    
    # read out meta data for pairwise comparison
    pairwise_metadata[[i_c]] <- condmetadata[bool_index_conditions, ]
  }
  
  cell_name <- vector()
  for (cond in 1:length(cells)) 
  {
    cell_name[cond] <- tail(strsplit(cells[cond], "[.]")[[1]], n=1)
  }
  
  cell_name <- factor(cell_name, levels=name_cond)
  names(cell_name) <- colnames(countmatrix)  
  pairwise_info <- list(pairwise_countmatrix, pairwise_metadata, cell_name)
  return(pairwise_info)
}


prepare_data_dgeanalysis <- function(matrix, metadata, num_samples_donor, design.function) 
{
  # 2. Prepare count matrix and metaData
  # Add Pseudo-count value of 1 to count matrix 
  # --> no removing of any broken samples necessary 
  # --> only in edgeR: DESeq2 need raw counts
  exprmatrix.groups <- prep_GEx_matrix(expr_matrix = matrix)
  metadata <- prep_metaData(metadata = metadata,
                            sample_id = colnames(exprmatrix.groups[[1]]),
                            groups = exprmatrix.groups[[2]])
  
  # check if columns of Count Matrix are equal to rows of MetaData
  bool_gene_matrix_metadata <- all(rownames(metadata) %in% 
                                     colnames(exprmatrix.groups[[1]]))
  if (bool_gene_matrix_metadata == FALSE)
  {
    exprmatrix.groups[[1]] <- exprmatrix.groups[[1]][, rownames(metadata)]
    print(all(rownames(metadata) == colnames(exprmatrix.groups[[1]])))
  }
  
  # 3. create design matrix
  # TAdd cellular detection rate (cdr)
  cdr <- scale(colMeans(exprmatrix.groups[[1]] > 0))
  metadata['cdr'] = cdr
  design.matrix = model.matrix(design.function, metadata)
  
  # convert metaData column to factors
  metadata['patient'] <- lapply(metadata['patient'], as.factor)
  metadata['disease'] <- lapply(metadata['disease'], as.factor)
  metadata['biopsy_type'] <- lapply(metadata['biopsy_type'], as.factor)
  metadata['sample'] <- lapply(metadata['sample'], as.factor)
  metadata['batch'] <- lapply(metadata['batch'], as.factor)
  metadata['annotation'] <- lapply(metadata['annotation'], as.factor)
  
  return(list(exprmatrix.groups[[1]], metadata, design.matrix))
}


prep_GEx_matrix <- function(expr_matrix) 
{
  # Prepare Expression matrix 
  # define control and replicate samples and create condition and label vectors
  sample_names <- colnames(expr_matrix)
  vec_cond = vector()
  sample_name_short = vector()
  for (cond in 1:length(sample_names))
  {
    temp_conditionnames <- strsplit(sample_names[cond], "[.]")[[1]]
    vec_cond[cond] <- tail(temp_conditionnames, n = 1)
    sample_name_short[cond] <- paste(temp_conditionnames[1], temp_conditionnames[2],
                                     tail(strsplit(vec_cond[cond], "_")[[1]], n = 1), sep = "_")  
  }
  
  # create condition comparisons
  # replace levels by labeling samples with "reference" and "control"
  str_group = vector()
  comp <- unique(vec_cond)
  freq_conditions = table(vec_cond)
  for (sa in vec_cond)
  {
    if (sa == comp[which.max(freq_conditions)[[1]]])
    {
      # set group with the most spots as reference
      str_group <- c(str_group, "reference") #" reference" "Cytokine"
    }
    else
    {
      str_group <- c(str_group,  "control") # "control" "Others"
    }
  }
  colnames(expr_matrix) <- sample_name_short
  
  return(list(expr_matrix, str_group))
}


prep_metaData <- function(metadata, sample_id, groups) 
{
  # Prepare MetaData 
  # MetaData dataframe with columns id, condition, label (needed to create design matrix)
  metaData <- subset(metadata, select=-c(sample_id, condition))
  metaData <- tibble::add_column(metaData, condition = groups, .after = "barcode")
  metaData <- tibble::add_column(metaData, sample = metadata$sample_id, .after = "barcode")
  metaData['annotation'] = metadata['label']
  metaData['label'] = metadata['condition']
  rownames(metaData) <- sample_id
  
  return(metaData)
}
