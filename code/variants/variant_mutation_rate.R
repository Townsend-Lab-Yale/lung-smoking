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

variant_mutation_rate = function(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, ref_genome_version, variants=NULL, maf_df=NULL, samples=NULL){
  gene_x = gene

  if(ref_genome_version == "hg19"){
    reference_genome = BSgenome.Hsapiens.UCSC.hg19}
  else if(ref_genome_version == "hg38"){
    reference_genome = BSgenome.Hsapiens.UCSC.hg38}
  else{stop('Reference genome version must be one of `hg19` or `hg38`')}

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
    if(length(variants) == 0){warning('There are no variants present in the dataset corresponding to this gene'); return(NULL)}
  } else if(length(variants) == 0){warning('The list of variants provided is empty.'); return(NULL)}
  
  unique_variants = unique(variants)
  
  tmp1 = str_split_fixed(unique_variants, pattern = ':', n=2)
  tmp2 = str_split_fixed(tmp1[,2], pattern= '_', n=2)
  variant_mut_rate_df = as.data.table(cbind(tmp1[,1],tmp2))
  colnames(variant_mut_rate_df) = c('chr','pos','mut')
  variant_mut_rate_df[,pos := strtoi(pos)]
  
  #get trinuc context for variants
  
  #' 
  #' ASSUMPTION: all sequences taken from POSITIVE strand
  #' 
  variant_mut_rate_df[, trinuc := as.vector(as.character(getSeq(reference_genome, names = str_c('chr', chr), start = pos-1, width=3, strand = '+')))]
  
  #reverse complementing variants as necessary
  middle = str_c(substr(variant_mut_rate_df$trinuc,2,2),collapse="")
  revcomp_ind = c(str_locate_all(middle,'G|A')[[1]][,1])
  variant_mut_rate_df[revcomp_ind, ':=' (trinuc = as.character(reverseComplement(DNAStringSet(variant_mut_rate_df[revcomp_ind, trinuc]))), 
                              mut = paste0(as.character(reverseComplement(DNAStringSet(variant_mut_rate_df[revcomp_ind, substr(mut,1,1)]))),
                                           '>',
                                           as.character(reverseComplement(DNAStringSet(variant_mut_rate_df[revcomp_ind, substr(mut,3,3)])))))]
  
  #generating COSMIC type format for trinuc mutations
  variant_mut_rate_df[,trinuc_mut := paste0(substr(trinuc,1,1),'[', mut,']',substr(trinuc,3,3))]
  
  #gets trinuc composition for each variant
  variant_mut_rate_df[,trinuc_mut_prop := as.vector(trinuc_proportion_matrix[trinuc_mut])]
  
  variant_mut_rate_df[,gene_rate := gene_syn_mutation_rate]
  
  tri_nt_contexts = table(variant_mut_rate_df$trinuc)
  variant_mut_rate_df[, context_freq := tri_nt_contexts[trinuc]]
  
  #' #' This only works with the singular table returned by compute tri-nt contexts, 
  #' #' if we were to vectorize this, we would need to make it dataframe-workable
  #' tri_nt_contexts = compute_trinucleotide_contexts(gene_x)
  #' variant_mut_rate_df[,context_freq := tri_nt_contexts[trinuc]] # should use data.table notation
  
  #calculate mutation rate for variants
  variant_mut_rate_df[,mut_rate := gene_rate * trinuc_mut_prop / context_freq]
  
  return(variant_mut_rate_df)
}

gene_mutation_rate = function(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, ref_genome_version, variants=NULL, maf_df=NULL, samples=NULL){
  return(sum(variant_mutation_rate(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, ref_genome_version, variants, maf_df, samples)$mut_rate))
  
  
  # variant_mut_rate_df = variant_mutation_rate(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, variants, maf_df, samples)
  # if(!is.null(variant_mut_rate_df)){
  #   return(sum(variant_mut_rate_df$mut_rate))}
  # else{return(NA)}
}
