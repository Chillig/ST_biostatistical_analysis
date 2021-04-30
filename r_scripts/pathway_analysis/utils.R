library("xlsx")
library('DOSE')
library('org.Hs.eg.db')


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


save_data_files <- function(paenrich_df, condition, pa_database, output_path) 
{
  xlsx_data =  paste(condition, pa_database, 'Pathway_Enrichment_Analysis', '.xlsx', sep = "_")
  csv_data = paste(condition, pa_database, 'Pathway_Enrichment_Analysis', '.csv', sep = "_")
  if (nrow(paenrich_df) > 1) 
  {
    write.xlsx(paenrich_df, file.path(output_path, xlsx_data), row.names = FALSE)
    write.csv(paenrich_df, file.path(output_path, csv_data), row.names = FALSE)
  }
}


enrichobject_to_df <- function(paenrich_object, condition, pa_database, output_path) 
{
  df_enrich = paenrich_object@result 
  save_data_files(paenrich_df=df_enrich, condition=condition, pa_database=pa_database, 
                  output_path=output_path)
}

gseaobject_to_df <- function(gsea_object_res, condition, pa_database, output_path)
{
  # exclude core_enrichment and leading_edge
  df_gsea = gsea_object_res@result[1:(length(gsea_object_res@result) - 2)]
  df_gsea$core_enrichment = gsea_object_res@result$core_enrichment
  save_data_files(paenrich_df=df_gsea, condition=condition, pa_database=pa_database, 
                  output_path=output_path)
}
