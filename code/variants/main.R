#'
#' Creates the MAF file and mutation rates necessary for input into cancer epistasis analysis. It assumes that the working directory is set to variants
#'

.libPaths(c("./.Rlibs", .libPaths()))

#' Load in relevant functions
#source('../cadd/cadd.R')
source('trinucleotide_functions.R')
source('variant_mutation_rate.R')
source('produce_genes_per_sample.R')
source('maf_construction.R')

#~~~~~~~~~~~~#

location_data = '../../data/'
location_output = '../../output/'

save_results = TRUE

#' Create CESA object for mutation rate calculation and MAF construction
#' Output locations: 'data/pan_data_cesa_for_cancer_epistasis.rds'
#'   'data/(panel_)(non)smoking_sample_ids.txt'
#'   'output/(non)smoking_(w_panel_)mutation_rates.txt'
source('create_cesa_for_epistasis.R')

#' List of genes for which to calculate variant-level mutation rates
gene_df = fread(paste0(location_data, "genes_list.txt"), header = F)
colnames(gene_df) = 'gene'
gene_df[,gene := toupper(gene)]
#' For now, only calculating mutation rates for the limited selection of genes that we are calculating fluxes for
gene_df = gene_df[1:103,]

#' Calculate trinucleotide mutation proportions across the whole genome in preparation for calculating variant mutation rates
genome_trinucleotide_mutation_proportions = genome_trinuc_mut_proportions(cesa$maf)
smoking_genome_trinucleotide_mutation_proportions = genome_trinuc_mut_proportions(cesa$maf[Unique_Patient_Identifier %in% 
                                                                                             c(smoking_samples, panel_smoking_samples)])
nonsmoking_genome_trinucleotide_mutation_proportions = genome_trinuc_mut_proportions(cesa$maf[Unique_Patient_Identifier %in% 
                                                                                             c(nonsmoking_samples, panel_nonsmoking_samples)])

#' Calculate variant level mutation rates
#' Then calculate the gene level mutation rate as the sum of all the variant level mutation rates within the gene
gene_df[, pan_data_mutation_rates := lapply(gene, function(x) sum(variant_mutation_rate(x, 
                                                                                        cesa, 
                                                                                        genome_trinucleotide_mutation_proportions, 
                                                                                        samples = cesa$samples$Unique_Patient_Identifier)))]
gene_df[, smoking_mutation_rates := lapply(gene, function(x) sum(variant_mutation_rate(x, 
                                                                                       cesa_smoking_w_panel,
                                                                                       smoking_genome_trinucleotide_mutation_proportions, 
                                                                                       samples = c(smoking_samples, panel_smoking_samples))))]
gene_df[, nonsmoking_mutation_rates := lapply(gene, function(x) sum(variant_mutation_rate(x, 
                                                                                          cesa_nonsmoking_w_panel,
                                                                                          nonsmoking_genome_trinucleotide_mutation_proportions, 
                                                                                          samples = c(nonsmoking_samples, panel_nonsmoking_samples))))]

fwrite(gene_df, paste0(location_output,"variant_based_mutation_rates.txt"))

#' Produce MAF for cancer epistasis analysis
#' Output location: 'data/cesR_maf_for_epistasis_analysis.txt'
preloaded_maf = rbind(rbindlist(Broad_maf), FMAD_maf, rbindlist(Genie_maf), MSK2015_maf, rbindlist(MSK2017_maf), rbindlist(MSK2018_maf), OncoSG_maf, TCGA_maf, TracerX_maf, TSP_maf, NCI_maf,fill = T)
final_maf = construct_maf(cesa, maf_file, preloaded_maf, save_results)

#' Additionally, create genes per sample table for new compute_samples functionality
#' Output location: 'output/genes_per_sample.txt'
samples_genes = produce_genes_per_sample(final_maf, gene_df$gene, save_results)
