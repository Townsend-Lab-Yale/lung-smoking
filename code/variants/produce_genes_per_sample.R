library(data.table)

produce_genes_per_sample = function(maf, gene_list, save_results = TRUE){
  maf = maf[,.(`Sample ID`,Source,Panel,Gene)]
  #some genes are mutated more than once but we only consider whether a gene is mutated or not
  maf = unique(maf)
  #remove all genes that we do not want to consider (and samples whose mutations only occur in those genes)
  maf = maf[Gene %in% toupper(gene_list)]

  #' Pivot to make a table where each row is a sample and the columns are genes
  #' with the value being a binary representation of whether the sample has the
  #' mutation or not
  samples_genes = dcast(maf, `Sample ID` + Source + Panel ~ Gene, fill = 0)
  temp = samples_genes[,-c('Sample ID','Source','Panel')]
  temp[temp != 0] <- 1
  samples_genes = cbind(samples_genes[,.(`Sample ID`,Source,Panel)], temp)

  if(save_results){
    fwrite(samples_genes, paste0(location_output, 'genes_per_sample.txt'))
  }

  return(samples_genes)
}
