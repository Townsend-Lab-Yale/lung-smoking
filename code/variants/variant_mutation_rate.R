# pick gene
# get rate
# collect all variants of interest for gene, reverse complementing as necessary
# find freq of that variant type (trinuc_mut) for the gene
# find #times that context is observed
# var_mut_rate = gene_rate * freq / #times-context-observed


variant_mutation_rate = function(gene, cadd_variants, cesa, trinuc_proportion_matrix, cadd_threshold = 20){
  gene_x = gene
  # unlist because some codons have multiple variants
  variants = c(unlist(cesa$variants[gene == gene_x & variant_type == 'aac',constituent_snvs]), cesa$variants[gene == gene_x & variant_type == 'snv',variant_id])
  
  #' 
  #' ASSUMPTION: Relies on CADD providing scores for ALL variants of interest
  #' 
  variants = variants[variants %in% cadd_variants[gene == gene_x & PHRED > cadd_threshold, variant_id]]
  
  temp1 = str_split_fixed(variants, pattern = ':', n=2)
  temp2 = str_split_fixed(temp1[,2], pattern= '_', n=2)
  variants = as.data.table(cbind(temp1[,1],temp2))
  colnames(variants) = c('chr','pos','mut')
  variants$pos = strtoi(variants$pos)
  
  #get trinuc context for variants
  
  #' 
  #' ASSUMPTION: all sequences taken from POSITIVE strand
  #' 
  variants$trinuc = as.vector(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, names = str_c('chr', variants$chr), start = variants$pos-1, width=3, strand = '+')))
  
  #reverse complementing variants as necessary
  middle = str_c(substr(variants$trinuc,2,2),collapse="")
  revcomp_ind = c(str_locate_all(middle,'G|A')[[1]][,1])
  variants[revcomp_ind, ':=' (trinuc = as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, trinuc]))), 
                              mut = paste0(as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, substr(mut,1,1)]))),
                                           '>',
                                           as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, substr(mut,3,3)])))))]
  
  #generating COSMIC type format for trinuc mutations
  variants$trinuc_mut = paste0(substr(variants$trinuc,1,1),'[', variants$mut,']',substr(variants$trinuc,3,3))
  
  #gets trinuc composition for each variant
  
  #' Previous implementation: has different trinucleotide mutation context proportions, resulting in differences in calculated mutation rates
  #' trinuc_mut_index_table = colnames(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix)
  #' ind = match(variants$trinuc_mut, trinuc_mut_index_table)
  #' variants$trinuc_mut_prop= as.vector(unlist(lapply(ind, function(y) cancereffectsizeR:::.ces_ref_data[[cesa@ref_key]][["gene_trinuc_comp"]][[gene_x]][y])))
  
  variants$trinuc_mut_prop = as.vector(trinuc_proportion_matrix[variants$trinuc_mut])
  variants$gene_rate = cesa$gene_rates[gene == gene_x,rate]
  
  #' This only works with the singular table returned by compute tri-nt contexts, if we were to vectorize this, would need to make it dataframe-workable
  tri_nt_contexts = compute_trinucleotide_contexts(gene_x)
  variants$context_freq = tri_nt_contexts[match(variants$trinuc, names(tri_nt_contexts))]
  
  #calculate mutation rate for variants
  variants[,mut_rate := gene_rate * trinuc_mut_prop / context_freq]
  
  return(variants$mut_rate)
}
