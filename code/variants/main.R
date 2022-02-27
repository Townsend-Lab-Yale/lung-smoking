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
cadd_table = fread('../../temp/GRCh37-v1.6_47f52a54c6b1adfe91e9afd010da4d32.tsv')
cesa = load_cesa('../../data/pan-dataset_samples_cesa.rds')
print('integrating cadd')
cadd_scores = integrate_cadd_scores(cesa, cadd_table)

gene_df = as.data.table(c('KRAS','TP53','LRP1B','BRAF','STK11'))
colnames(gene_df) = 'gene'
print('computing gene mutation rates')
gene_df$mutation_rates = lapply(gene_df$gene, function(x) sum(variant_mutation_rate(x, cadd_scores, cesa)))