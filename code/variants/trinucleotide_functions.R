#' Gets proportion of each trinucleotide mutation (N[N>N]N) in the dataset
#' Any trinucleotide mutations not in the output table have a proportion of 0
#' Result is heavily dependent on dataset as the count for each trinucleotide mutation
#'   comes from the set of variants in the dataset. Thus biases which influence
#'   which variants are represented in the dataset will also influence the proportions.

compute_trinucleotide_mut_proportions = function(gene, cesa){
  gene_x = gene
  variants = c(unlist(cesa$variants[gene == gene_x & variant_type == 'aac',constituent_snvs]), cesa$variants[gene_x == gene & variant_type == 'snv',variant_id])
  
  temp1 = str_split_fixed(variants, pattern = ':', n=2)
  temp2 = str_split_fixed(temp1[,2], pattern= '_', n=2)
  variants = as.data.table(cbind(temp1[,1],temp2))
  colnames(variants) = c('chr','pos','mut')
  variants$pos = strtoi(variants$pos)
  
  variants$trinuc = as.vector(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, names = str_c('chr', variants$chr), start = variants$pos-1, width=3, strand = '+')))
  
  #reverse complementing variants as necessary
  middle = str_c(substr(variants$trinuc,2,2),collapse="")
  revcomp_ind = c(str_locate_all(middle,'G|A')[[1]][,1])
  variants[revcomp_ind, ':=' (trinuc = as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, trinuc]))), mut = paste0(as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, substr(mut,1,1)]))),'>',as.character(reverseComplement(DNAStringSet(variants[revcomp_ind, substr(mut,3,3)])))))]
  
  #generating COSMIC type format for trinuc mutations
  variants$trinuc_mut = paste0(substr(variants$trinuc,1,1),'[', variants$mut,']',substr(variants$trinuc,3,3))
  
  trinuc_mut_proportions = table(variants$trinuc_mut) / nrow(variants)
  
  return(trinuc_mut_proportions)
}



#' Calculates the 32 trinucleotide contexts for every site in any gene of interest in the hg19 genome.
#' 32 instead of 64 because any purine site is equivalent to its corresponding pyrimidine site 
#'  because of the way mutations to those sites would be interpreted (see COSMIC)
#' For purine sites, their context numbers are added to their reverse complement

compute_trinucleotide_contexts <- function(genes) {
  coding_seqs = sapply(genes, function(x) DNAString(paste0(as.character(cancereffectsizeR:::.ces_ref_data$ces.refset.hg19$RefCDS[[x]]$seq_cds1up[1:2]), as.character(cancereffectsizeR:::.ces_ref_data$ces.refset.hg19$RefCDS[[x]]$seq_cds1down))))
  tri_nt_contexts = sapply(coding_seqs, function(x) Biostrings::trinucleotideFrequency(x))
  rows_to_remove = c()
  for(i in 1:nrow(tri_nt_contexts)){
    if(substr(rownames(tri_nt_contexts)[i],2,2) %in% c('G','A')){
      rows_to_remove = c(rows_to_remove, i)
      tri_nt_contexts[as.character(reverseComplement(DNAString(rownames(tri_nt_contexts)[i]))),] =
        tri_nt_contexts[as.character(reverseComplement(DNAString(rownames(tri_nt_contexts)[i]))),] +
        tri_nt_contexts[i,]
    }
  }
  tri_nt_contexts = tri_nt_contexts[-rows_to_remove,]
  return(tri_nt_contexts)
}
