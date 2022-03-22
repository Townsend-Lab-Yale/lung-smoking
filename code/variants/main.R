library(data.table)
library(cancereffectsizeR)
library(Biostrings)
library(stringr)

source('../cadd/cadd.R')
source('trinucleotide_functions.R')
source('variant_mutation_rate.R')

# TO-DO: make into a function that accepts a gene list and returns a dataframe

# EXAMPLE USAGE
print('loading data')
cadd_table = fread('../../temp/cadd_scores.tsv')
cesa = load_cesa('../../data/pan-dataset_samples_cesa.rds')
maf_file = fread("../../output/merged_luad_maf.txt")


print('creating maf file with variant classification column and CADD scores')
#unsure if I need to preload each dataset separately
maf_file = preload_maf(maf_file, refset = ces.refset.hg19, chain_file = "../../data/hg38ToHg19.over.chain", sample_col = 'Sample ID', keep_extra_columns = T)
var_classification = maf_file[,.(variant_id, Variant_Classification)]
var_classification = var_classification[!duplicated(var_classification$variant_id)]
# MAF to be used in cancer epistasis analysis
cesa = integrate_cadd_scores(cesa, cadd_table)
cesa@maf = merge(cesa$maf, var_classification, by='variant_id', all.x = T)


#' now that we can integrate cadd scores into maf and want to avoid using variants table, 
#' will deprecate integrate_cadd_score in favor of merge_maf_cadd_scores
#' will need to modify variant_mutation_rate function accordingly
gene_df = as.data.table(c('KRAS','TP53','LRP1B','BRAF','STK11'))
colnames(gene_df) = 'gene'
print('computing gene mutation rates')
genome_trinucleotide_mutation_proportions = compute_genome_trinucleotide_mut_proportions(cesa)
gene_df$mutation_rates = lapply(gene_df$gene, function(x) sum(variant_mutation_rate(x, cesa, genome_trinucleotide_mutation_proportions)))
