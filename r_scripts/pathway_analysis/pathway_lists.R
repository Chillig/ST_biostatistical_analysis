#! /usr/bin/Rscript
library('hash')

#' Pathways which shall be shown in the paper
#' 
#' @return hash dictionary containing pathways for each DGE Analysis
pathwaysofinterest <- function() 
{
  pa_publication = hash::hash()
  # Single cell enriched Pathways
  pa_publication['sc_IFNG_pa'] <- c(
    'G alpha (i) signalling events', 'Signaling by Interleukins', 
    'Metabolism of steroids', 'Integrin cell surface interactions')
  pa_publication['sc_IL17A_pa'] <- c(
    'Signaling by Retinoic Acid', 'Interleukin-17 signaling', 
    'Fatty acid metabolism', 'Signaling by Interleukins')
  
  # Spatial transcriptomics enriched Pathways
  pa_publication['st_IL17A_pa'] <- c(
    'Signaling by Interleukins', 'Chemokine receptors bind chemokines', 
    'Antimicrobial peptides', 'Keratinization')
  pa_publication['st_IFNG_pa'] <- c(
    'Chemokine receptors bind chemokines',  'Signaling by Interleukins', 
    'Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell',
    'Interferon Signaling')
  pa_publication['st_IL13_pa'] <- c(
    'Chemokine receptors bind chemokines', 'Signaling by Interleukins', 
    'Interleukin-4 and Interleukin-13 signaling', 
    'Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell')
  
  return(pa_publication)
}
