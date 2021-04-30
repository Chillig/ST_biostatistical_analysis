library('hash')
library(dplyr)

get_hkg <- function(genes)
{
  #' 
  #' @description Get Housekeeping genes and find them in provided gene list
  #' @param genes gene list
  #' @return hkg and index of hkg in provided gene list
  # house keeping genes: https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
  house_keeping_genes = c('ACTB', 'GAPDH', 'PGK1', 'PPIA', 'RPLP0', 
                          'ARBP', 'B2M', 'YWHAZ', 'SDHA', 
                          'TFRC', 'GUSB', 'HMBS', 'HPRT1', 'TBP')
  df_gene_names <- data.frame(genes)
  ind_house_keep_genes <- row.names.data.frame(df_gene_names)[genes %in% house_keeping_genes]
  ind_house_keep_genes <- strtoi(ind_house_keep_genes)
  
  hkg <- list(house_keeping_genes, ind_house_keep_genes)
  
  return(hkg)
}


get_cyto.responders <- function() 
{
  #' @note Responder genes have been derived using the Kerasarray 
  #' (Keratinocytes which are only localised in the EPIDERMIS)
  #' @description Get Cytokines and their responder genes
  #' @param genes gene list
  #' @return cytokines and responder gene list
  # Responder and Cytokine genes
  cyto.responders = c("IL17A", "SPRR2B", "LCN2", "IL19", "CXCL6", "SPRR2E", "IGF2",
                      "SPRR2A", "SPRR2F", "SPRR2D", "DEFB4A",
                      "IFNG", "UBD", "CXCL9", "CXCL10", "CXCL11", "IL32", "ICAM1",
                      "BATF2", "GBP5", "GBP4", "CCL8", "IRF1", "IRF1", "RSAD2", "MUC1",
                      "IL13", "CCL26", "NTRK1", "HSD3B1", "BET3L", "SERPINB13", "CH25H",
                      "EML5", "P2RY1", "LOXL4", "FAM26D", "THBS1")
  
  return(cyto.responders)
}


gsg.epidermis <- function() 
{
  goldenstandard_genes = hash::hash()
  # Source LCE: https://www.genenames.org/data/genegroup/#!/group/627
  goldenstandard_genes[["uE"]] <- c("FLG", "HRN", 
                                    "LCE1A", "LCE1B", "LCE1C", "LCE1D", "LCE1E", "LCE1F",
                                    "LCE2A", "LCE2B", "LCE2C", "LCE2D",
                                    "LCE3A", "LCE3B", "LCE3C", "LCE3D", "LCE3E",
                                    "LCE4A", "LCE5A", "LCE6A")
  goldenstandard_genes[["mE"]] = c("FLG", "HRN", 
                                   "LCE1A", "LCE1B", "LCE1C", "LCE1D", "LCE1E", "LCE1F",
                                   "LCE2A", "LCE2B", "LCE2C", "LCE2D",
                                   "LCE3A", "LCE3B", "LCE3C", "LCE3D", "LCE3E",
                                   "LCE4A", "LCE5A", "LCE6A")
  # Source S100: https://www.genenames.org/data/genegroup/#!/group/459
  goldenstandard_genes[["bE"]] <- c("KRT10", "DEFB4", 
                                    "S100A1", "S100A2", "S100A3", "S100A4", "S100A5",
                                    "S100A6", "S100A7", "S100A7A", "S100A7L2", "S100A7P1",
                                    "S100A7P2", "S100A8", "S100A9",
                                    "S100A10", "S100A11", "S100A12", "S100A13", "S100A14",
                                    "S100A15A", "S100A16",
                                    "S100B", "S100G", "S100P", "S100Z")
  
  return(goldenstandard_genes)
}
