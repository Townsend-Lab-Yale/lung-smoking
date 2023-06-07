library(data.table)
library(cancereffectsizeR)
library(Biostrings)
library(stringr)


# pick gene
# get rate
# collect all variants of interest for gene, reverse complementing as necessary
# find freq of that variant type (trinuc_mut) for the gene
# find #times that context is observed
# var_mut_rate = gene_rate * freq / #times-context-observed

variant_mutation_rate = function(gene, cesa, trinuc_proportion_matrix, samples=NA){
  gene_x = gene
  print(gene_x)
  # certain variants have an NA variant_id, but none of these are in genes in our 1290 genes of interest.
  maf = cesa$maf
  
  if(any(is.na(samples))){stop('Please enter a list of samples. If you want to include all samples, use cesa$samples$Unique_Patient_Identifier')}
  
  variants = unique(cesa$maf[Unique_Patient_Identifier %in% samples & top_gene == gene_x & 
                               variant_type == 'snv' & !is.na(variant_id), 
                             variant_id])
  
  if(length(variants) == 0){return(NA)}

  
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
  variants$trinuc_mut_prop = as.vector(trinuc_proportion_matrix[variants$trinuc_mut])
  
  if(colnames(cesa$gene_rates)[2] == "rate"){
    variants$gene_rate = cesa$gene_rates[gene == gene_x, rate]
  } else if(colnames(cesa$gene_rates)[2] == "rate_grp_1"){
    variants$gene_rate = cesa$gene_rates[gene == gene_x, rate_grp_1]
  } else{stop('Gene mutation rate column not recognized')}
  
  
  #' This only works with the singular table returned by compute tri-nt contexts, 
  #' if we were to vectorize this, we would need to make it dataframe-workable
  tri_nt_contexts = compute_trinucleotide_contexts(gene_x)
  variants$context_freq = tri_nt_contexts[variants$trinuc]
  
  #calculate mutation rate for variants
  variants[,mut_rate := gene_rate * trinuc_mut_prop / context_freq]
  
  return(variants$mut_rate)
}
