library("xlsx")
library('DOSE')
library('org.Hs.eg.db')


filter_data <- function(df.data, signature_gene,plot_signaturecytokine ) 
{
  # Remove redundant information:
  # # -remove duplicated rows
  df.data <- df.data[!duplicated(df.data$gene_symbol), ]
  rownames(df.data) <- df.data$gene_symbol
  df.data <- df.data[, c(colnames(df.data) != "X"), drop = FALSE]
  # # -remove NA columns
  df.data <- df.data[, colSums(is.na(df.data)) < nrow(df.data)]
  # # -remove NA rows
  df.data <- na.omit(df.data)
  
  # If cytokine shall be considered in analysis 
  if (plot_signaturecytokine == FALSE) 
  {
    # remove signature cytokine for plots
    df.data = df.data[df.data$gene_symbol != signature_gene, ]
  } 
  
  return(df.data)
}


rename_genetoentrezid <- function(df.dge_results) 
{
  # I. replace gene symbol by entrezid of gene
  if ("entrezid" %nin% colnames(df.dge_results))
  {
    # works approach 
    Entrez_ID <- mapIds(org.Hs.eg.db, row.names(df.dge_results), 'ENTREZID', 'SYMBOL')
    Entrez_ID_df = list2df(Entrez_ID, col1 = "entrezid", col2 = "gene_symbol")
    Entrez_ID_df = Entrez_ID_df[!is.na(Entrez_ID_df$entrezid), ]
    
    gseaDat <- filter(df.dge_results, !is.na(Entrez_ID))
    df.dge_results <- filter(df.dge_results, !is.na(Entrez_ID))
    df.dge_results$entrezid <- Entrez_ID_df$entrezid
  }
  
  return(df.dge_results)
}

do_rank_genes <- function(df.dge_results) 
{
  # II. Rank all genes based on their fold change
  ranked_genes <- df.dge_results$log2fc
  names(ranked_genes) <- df.dge_results$gene_symbol
  ranked_genes <- sort(ranked_genes, decreasing = T)
  
  return(ranked_genes)
}


get_significantgenes <- function(df.dge_results, p_value, lfc_factor, op) 
{
  # II. Get significantly differentially expressed genes
  df.sig <- df.dge_results[df.dge_results$pval < p_value & 
                             !is.na(df.dge_results$pval) & 
                             op(df.dge_results$log2fc, lfc_factor), ]
  degenes.sig <- df.sig$entrezid
  degenes.sig <- as.character(na.exclude(degenes.sig))       
  
  return(list(df.sig, degenes.sig))
}


setreadable_pa <- function(paenrich_object) 
{
  if (!is.null(paenrich_object) && dim(paenrich_object)[1] > 0 && dim(paenrich_object)[2] > 0)
  {
    paenrich_object <- DOSE::setReadable(paenrich_object, 'org.Hs.eg.db', 'ENTREZID')
  }else{
    paenrich_object <- paenrich_object
  }
  return(paenrich_object)
}


save_enrichobject_as_csv <- function(paenrich_object, condition, pa_database, output_path) 
{
  df_enrich = paenrich_object@result 
  xlsx_data =  paste(condition, pa_database, 'Pathway_Enrichment_Analysis', '.xlsx', sep = "_")
  csv_data = paste(condition, pa_database, 'Pathway_Enrichment_Analysis', '.csv', sep = "_")
  if (nrow(df_enrich) > 1) 
  {
    write.xlsx(df_enrich, file.path(output_path, xlsx_data), row.names = FALSE)
    write.csv(df_enrich, file.path(output_path, csv_data), row.names = FALSE)
  }
}

