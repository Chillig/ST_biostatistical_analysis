source('../scripts/helper_functions.R')

# ------------------------- Load files -------------------------
load_metadata <- function(metadata.path)
{
  #' Load MetaData as data frame
  #'
  #' @note The .csv file entries should be separated by "," or ";"
  #' @param metadata.path Path to .csv or .xlsx MetaData file
  #' @return data frame containing the information about the MetaData
  
  message("-- Load metaData --")
  
  # read in metaData from csv file
  if (tail(strsplit(metadata.path, split = "[.]")[[1]], n = 1) == "xlsx") 
  {
    df_metaData <- xlsx::read.xlsx(metadata.path, sheetIndex = 1, header = TRUE) 
  } else {
    fl <- readLines(metadata.path, n = 1)
    if (grepl(";", fl)) 
    {
      df_metaData <-  read.csv2(metadata.path, header = TRUE, sep = ";")
    }else 
    {
      df_metaData <- read.csv(metadata.path, header = TRUE, sep = ",")
    }
  }
  
  # remove columns with only indexes as entries
  df_metaData <- df_metaData[, c(colnames(df_metaData) != "X"), drop = FALSE]
  
  return(df_metaData)
}


load_countmatrix <- function(path_file.name) 
{
  #' Read in .csv or .xlsx file 
  #' 
  #' Gets the full path and reads in csv file including the header
  #' @note The .csv file entries should be separated by "," or ";"
  #' @param path_file.name Path to .csv or .xlsx file
  #' @return Dataframe and Genes
  
  message("-- Load count matrix --")
  
  if (tail(strsplit(path_file.name, split = "[.]")[[1]], n = 1) == "xlsx") 
  {
    df_cond <- xlsx::read.xlsx(path_file.name, sheetIndex = 1, header = TRUE)
  } else {
    fl <- readLines(path_file.name, n = 1)
    if (grepl(";", fl)) 
    {
      df_cond <-  read.csv2(path_file.name, header = TRUE, sep = ";", row.names = 'geneNames')
    }else 
    {
      df_cond <- read.csv(path_file.name, header = TRUE, sep = ",", row.names = 'geneNames')
    }
  }
  
  gene_names <- make.unique(row.names(df_cond))
  rownames(df_cond) <- gene_names
  
  # remove columns with only zeros as entries
  df_cond <- df_cond[, c(colnames(df_cond) != "X"), drop = FALSE]
  
  if (ncol(df_cond) > 0)
  {
    df_cond <- delete_zero_conditions(df_cond)
  }
  
  return(df_cond)
}
