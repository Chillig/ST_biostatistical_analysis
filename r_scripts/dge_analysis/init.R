# Preparation

get_inputdir <- function(input.dir, date.file, dataset.type, seq.technique, comparison) 
{
  input.dir <- file.path(input.dir, paste(date.file, dataset.type, sep = '__'), 
                         seq.technique, comparison)
  return(input.dir)
}

get_outputdir <- function(dataset.type, dge.method, seq.technique) 
{
  #'
  #' @param method
  #' @param project
  #' @return output path
  
  # output paths
  output_dir <- file.path(getwd(), "Output", dataset.type, paste(dge.method, "output", sep = "_"), 
                          Sys.Date())
  dir.create(output_dir, showWarnings = FALSE)
  output_path <- file.path(output_dir, seq.technique)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  return(output_path)
}

