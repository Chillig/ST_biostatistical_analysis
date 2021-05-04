library("xlsx")
# ------------------------- Load files -------------------------
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
    dge_df <- read.csv(path_name_file, header = TRUE, sep = ";")
    
    if (ncol(dge_df) == 1) {dge_df = read.csv(path_name_file, header = TRUE, sep = ",")}
  }
  
  return(dge_df)
}
