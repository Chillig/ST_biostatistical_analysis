# Preparation

#' Get input directory
#' 
#' @param input.folder
#' @param date.file
#' @param dataset.type
#' @param seq.technique
#' @param comparison
#' @param genename
#' @param design.function
get_inputdir <- function(input.folder, date.file, dataset.type, seq.technique, comparison,
                         design.function, genename) 
{
  input.dir <- file.path(input.folder, date.file, dataset.type, 
                         seq.technique,paste(seq.technique, design.function, sep = "_"), genename,
                         comparison)
  return(input.dir)
}


#' Create output directory
#' 
#' @param output.folder
#' @param dataset.type
#' @param seq.technique
#' @param genename
#' @return output path
get_outputdir <- function(output.folder, dataset.type, seq.technique, genename) 
{
  # output paths
  # output_dir <- file.path("..", "..", "output", "Figure_3K", dataset.type,  Sys.Date())
  output_dir <- file.path(output.folder, "output", "Figure_3K", Sys.Date(), dataset.type, genename)
  dir.create(output_dir, showWarnings = FALSE)
  output_path <- file.path(output_dir, seq.technique)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  return(output_path)
}

