library(data.table)
library(cancereffectsizeR)
library(Biostrings)
library(stringr)


order_of_trinuc_muts = 
   c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T",
  "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T",
  "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T",
  "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",
  "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T",
  "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
  "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T",
  "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",
  "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
  "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
  "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
  "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
  "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T",
  "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T",
  "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
  "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
  "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T",
  "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
  "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T",
  "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
  "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
  "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
  "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
  "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")


#' Gets proportion of each trinucleotide mutation (N[N>N]N) in the dataset
#' Any trinucleotide mutations not in the output table have a proportion of 0
#' Result is heavily dependent on dataset as the count for each trinucleotide mutation
#'   comes from the set of variants in the dataset. Thus biases which influence
#'   which variants are represented in the dataset will also influence the proportions.

compute_trinuc_mut_proportions <- function(variants){
  temp1 = str_split_fixed(variants, pattern = ':', n=2)
  temp2 = str_split_fixed(temp1[,2], pattern= '_', n=2)
  variants = as.data.table(cbind(temp1[,1],temp2))
  colnames(variants) = c('chr','pos','mut')
  variants$pos = strtoi(variants$pos)
  ## One sample did not have a position (strangely)
  variants = variants[!is.na(variants$pos)]
  variants$trinuc = as.vector(as.character(getSeq(
    BSgenome.Hsapiens.UCSC.hg19,
    names=str_c('chr', variants$chr),
    start=variants$pos-1,
    width=3,
    strand = '+')))
  ## reverse complementing variants as necessary
  middle = str_c(substr(variants$trinuc,2,2),collapse="")
  revcomp_ind = c(str_locate_all(middle,'G|A')[[1]][,1])
  variants[revcomp_ind,
           ':=' (trinuc=as.character(reverseComplement(DNAStringSet(
             variants[revcomp_ind, trinuc]))),
             mut=paste0(as.character(reverseComplement(DNAStringSet(
               variants[revcomp_ind, substr(mut,1,1)]))),
               '>',
               as.character(reverseComplement(DNAStringSet(
                 variants[revcomp_ind, substr(mut,3,3)])))))]
  ## generating COSMIC type format for trinuc mutations
  variants$trinuc_mut = paste0(substr(variants$trinuc, 1, 1),
                               '[', variants$mut, ']',
                               substr(variants$trinuc, 3, 3))
  ## concatenated table contains the names of trinuc_muts not
  ## represented in the dataset (subtract 1 because table is
  ## auto-iniated with 1s) divides number of observations by number of
  ## variants to get proportions
  trinuc_mut_proportions = (
    c(table(variants$trinuc_mut),
      table(setdiff(order_of_trinuc_muts, variants$trinuc_mut)) - 1)
    / nrow(variants))
  ## TODO: order table before returning value
  return(trinuc_mut_proportions)
}


#' Computes proportions of each trinucleotide mutation in the CESA MAF file across the entire genome
genome_trinuc_mut_proportions <- function(cesa){
  variants = unique(cesa$maf[!is.na(variant_id) & variant_type == 'snv', variant_id])
  return(compute_trinuc_mut_proportions(variants))
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