#' Creates a maf table of non-silent snv's with some variants removed via cesR filtering. To be used as the maf file for the cancer epistasis analysis.

library(data.table)

construct_maf = function(cesa_maf, original_maf, preloaded_maf, variants_per_gene=NULL, save_results = TRUE){
  #' Assigning a single gene to each variant
  cesa_maf = cesa_maf[lengths(genes) == 1 & !is.na(genes) & is.na(top_gene), top_gene := as.character(genes)]
  cesa_maf_variants = cesa_maf[variant_type == 'snv', .(Unique_Patient_Identifier, variant_id)]
  
  preload_maf_variants = preloaded_maf[variant_type == 'snv' & Variant_Classification != 'Silent' & is.na(problem) & !is.na(variant_id), .(Unique_Patient_Identifier, variant_id)]
  
  #' Taking the intersection of these loaded and preloaded variants as the ones to keep
  variants_to_keep = merge(cesa_maf_variants, preload_maf_variants) 
  #' Adding back samples that have no SNVs or only silent SNVs (and were therefore excluded by the filtering steps)
  variants_to_keep = rbind(variants_to_keep, 
                           data.frame(Unique_Patient_Identifier=setdiff(cesa_maf$Unique_Patient_Identifier, 
                                                                        variants_to_keep$Unique_Patient_Identifier)),
                           fill = T)
  #' Note: variants_per_gene should be a subset of variants_to_keep, with the only difference being that 
  #' variants_per_gene is restricted to the genes in gene_mut_rate_df (from genes_list.txt)
  
  #' Keeping desired variants in MAF and removing extraneous columns
  filtered_maf = merge(cesa_maf, variants_to_keep, by = c('Unique_Patient_Identifier','variant_id'), all.y=T) # must specify all.x or all.y = T because variant_id is NA for the samples added back
  filtered_maf = filtered_maf[,-c('prelift_chr', 'prelift_start', 'liftover_strand_flip','variant_type', 'genes', 'top_consequence')]

  #' Adding source information to all of the TGS samples
  source_panel_info = as.data.table(original_maf)
  source_panel_info = source_panel_info[!duplicated(Tumor_Sample_Barcode),.(Tumor_Sample_Barcode, Source, Panel)]
  source_panel_info[Panel == '', Panel := NA]
  
  filtered_maf = merge(filtered_maf, source_panel_info, by.x = 'Unique_Patient_Identifier', by.y = 'Tumor_Sample_Barcode', all.x = T)
  
  #' Renaming columns for consistency between Python and R
  colnames(filtered_maf) = c('Sample ID', 'Mutation', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Gene', 'Source', 'Panel')
  
  if(save_results){
    fwrite(filtered_maf, paste0(location_output, "cesR_maf_for_epistasis.txt"))
  }
  
  return(filtered_maf)
}
