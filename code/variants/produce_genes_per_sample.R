library(data.table)

produce_genes_per_sample = function(maf, gene_list, save_results = TRUE){
  maf = maf[,.(`Sample ID`,Source,Panel,Gene)]
  #some genes are mutated more than once but we only consider whether a gene is mutated or not
  maf = unique(maf)
  #remove all genes that we do not want to consider (also removes samples whose mutations only occur in those genes)
  subset_maf = maf[Gene %in% toupper(gene_list)]

  # adding the removed samples back
  subset_maf = rbind(subset_maf, 
                     dplyr::setdiff(maf[,.(`Sample ID`, Source, Panel)], subset_maf[,.(`Sample ID`, Source, Panel)]),
                     fill = T)
  subset_maf[is.na(subset_maf)] <- ''

  #' Pivot to make a table where each row is a sample and the columns are genes
  #' with the value being a binary representation of whether the sample has the
  #' mutation or not
  samples_genes = dcast(subset_maf, `Sample ID` + Source + Panel ~ Gene, fill = 0)
  samples_genes$V1 = NULL
  tmp = samples_genes[,-c('Sample ID','Source','Panel')]
  tmp[tmp != 0] <- 1
  samples_genes = cbind(samples_genes[,.(`Sample ID`,Source,Panel)], tmp)

  if(save_results){
    fwrite(samples_genes, paste0(location_output, 'genes_per_sample.txt'))
  }

  return(samples_genes)
}