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

variant_mutation_rate = function(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, variants=NULL, maf_df=NULL, samples=NULL){
  gene_x = gene
  # print(gene_x)
  
  # certain variants have an NA variant_id, but none of these are in genes in our 1290 genes of interest.
  
  if(is.null(variants)){
    print('No variant list provided. Using only the variants observed within the gene from the MAF dataframe provided.')
    if(is.null(maf_df)){stop("MAF dataframe must be provided if variant list is not provided")}
    if(is.null(samples)){stop("List of samples must be provided if variant list is not provided.
                              If you want to include all samples, use cesa$samples$Unique_Patient_Identifier")}
    variants = unique(maf_df[Unique_Patient_Identifier %in% samples & top_gene == gene_x & 
                               variant_type == 'snv' & !is.na(variant_id), 
                             variant_id])
    if(length(variants) == 0){warning('There are no variants present in the dataset corresponding to this gene'); return(NA)}
  } else if(length(variants) == 0){warning('The list of variants provided is empty.'); return(NA)}
  
  variants = unique(variants)
  
  tmp1 = str_split_fixed(variants, pattern = ':', n=2)
  tmp2 = str_split_fixed(tmp1[,2], pattern= '_', n=2)
  variants = as.data.table(cbind(tmp1[,1],tmp2))
  colnames(variants) = c('chr','pos','mut')
  variants[,pos := strtoi(pos)]
  
  #get trinuc context for variants
  
  #' 
  #' ASSUMPTION: all sequences taken from POSITIVE strand
  #' 
  variants[, trinuc := as.vector(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, names = str_c('chr', chr), start = pos-1, width=3, strand = '+')))]
  
  #reverse complementing variants as necessary
  middle = str_c(substr(variants$trinuc,2,2),collapse="")
  revcomp_ind = c(str_locate_all(middle,'G|A')[[1]][,1])
  variants[revcomp_ind, ':=' (trinuc = as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, trinuc]))), 
                              mut = paste0(as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, substr(mut,1,1)]))),
                                           '>',
                                           as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, substr(mut,3,3)])))))]
  
  #generating COSMIC type format for trinuc mutations
  variants[,trinuc_mut := paste0(substr(trinuc,1,1),'[', mut,']',substr(trinuc,3,3))]
  
  #gets trinuc composition for each variant
  variants[,trinuc_mut_prop := as.vector(trinuc_proportion_matrix[trinuc_mut])]
  
  variants[,gene_rate := gene_syn_mutation_rate]
  
  tri_nt_contexts = table(variants$trinuc)
  variants[, context_freq := tri_nt_contexts[trinuc]]
  
  #' #' This only works with the singular table returned by compute tri-nt contexts, 
  #' #' if we were to vectorize this, we would need to make it dataframe-workable
  #' tri_nt_contexts = compute_trinucleotide_contexts(gene_x)
  #' variants[,context_freq := tri_nt_contexts[trinuc]] # should use data.table notation
  
  #calculate mutation rate for variants
  variants[,mut_rate := gene_rate * trinuc_mut_prop / context_freq]
  
  return(variants$mut_rate)
}

gene_mutation_rate = function(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, variants=NULL, maf_df=NULL, samples=NULL){
  return(sum(variant_mutation_rate(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, variants, maf_df, samples)))
}
