#'
#' Creates the MAF file and mutation rates necessary for input into cancer epistasis analysis
#'

.libPaths(c("./.Rlibs", .libPaths()))

#' Load in relevant functions
#source('../cadd/cadd.R')
source('trinucleotide_functions.R')
source('variant_mutation_rate.R')

#~~~~~~~~~~~~#

#' NOTE: Different users will need to change this
location_data = '../../data/'
location_output = '../../output/'

#' Create CESA object for mutation rate calculation and MAF construction
#' Output location: 'data/pan_data_cesa_for_cancer_epistasis.rds'
source('create_cesa_for_epistasis.R')

#' List of genes for which to calculate variant-level mutation rates
gene_df = fread(paste0(location_output,"genes_list.txt"), header = F)
colnames(gene_df) = 'gene'
gene_df[,gene := toupper(gene)]

#' Calculate trinucleotide mutation proportions across the whole genome in preparation for calculating variant mutation rates
genome_trinucleotide_mutation_proportions = compute_genome_trinucleotide_mut_proportions(cesa)

#' Calculate variant level mutation rates
#' Then calculate the gene level mutation rate as the sum of all the variant level mutation rates within the gene
gene_df[, mutation_rates := lapply(gene, function(x) sum(variant_mutation_rate(x, cesa, genome_trinucleotide_mutation_proportions)))]
fwrite(gene_df, paste0(location_output,"variant_based_mutation_rates.txt"))

#' Produce MAF for cancer epistasis analysis
#' Output location: 'data/cesR_maf_for_epistasis_analysis.txt'
source('maf_construction.R')

#' Additionally, create genes per sample table for new compute_samples functionality
#' Output location: 'data/genes_per_sample.txt'
source('produce_samples_per_gene.R')



#~~~deprecated~~~#

# EXAMPLE USAGE WITH CADD
#print('loading data')
#cadd_table = fread('../../temp/cadd_scores.tsv')
#cesa = load_cesa('../../data/pan_data_cesa_for_cancer_epistasis.rds')
#maf_file = fread("../../output/merged_luad_maf.txt")

#print('creating maf file with variant classification column and CADD scores')
#unsure if I need to preload each dataset separately
##maf_file = preload_maf(maf_file, refset = ces.refset.hg19, chain_file = "../../data/hg38ToHg19.over.chain", sample_col = 'Sample ID', keep_extra_columns = T)
##var_classification = maf_file[,.(variant_id, Variant_Classification)]
##var_classification = var_classification[!duplicated(var_classification$variant_id)]
# MAF to be used in cancer epistasis analysis
##cesa = integrate_cadd_scores(cesa, cadd_table)
##cesa@maf = merge(cesa$maf, var_classification, by='variant_id', all.x = T)

#gene_df = fread('../../data/genes_list.txt', header = F)
#colnames(gene_df) = 'gene'
#gene_df[,gene := toupper(gene)]

#print('computing gene mutation rates')
#genome_trinucleotide_mutation_proportions = compute_genome_trinucleotide_mut_proportions(cesa)
#gene_df[, mutation_rates := lapply(gene, function(x) sum(variant_mutation_rate(x, cesa, genome_trinucleotide_mutation_proportions)))]
#fwrite(gene_df, '../../output/variant_based_mutation_rates.txt')
