#'
#' Creates the MAF file and mutation rates necessary for input into cancer epistasis analysis. It assumes that the working directory is set to variants
#'

.libPaths(c("./.Rlibs", .libPaths()))

#' Load in relevant functions
source('trinucleotide_functions.R')
source('variant_mutation_rate.R')
source('produce_genes_per_sample.R')
source('maf_construction.R')

#~~~~~~~~~~~~#

location_data = '../../data/'
location_output = '../../output/'

save_results = TRUE

#' Create CESA object for mutation rate calculation and MAF construction
#' Output locations: 
#'   'data/pan_data_cesa.rds'
#'   'data/(panel_)(non)smoking_sample_ids.txt'
#'   'output/(non)smoking_mutation_rates.txt'
print('Creating CESA object for mutation rate calculation and MAF cleaning')
source('create_cesa_for_epistasis.R')

#' List of genes for which to calculate variant-level mutation rates
subset_genes = c("TP53",    "KRAS",  "EGFR",  "BRAF",   "CTNNB1",
                "KEAP1",   "STK11", "ATM",   "PIK3CA", "RBM10",
                "SMARCA4", "SMAD4", "ALK",   "ARID1A", "APC",
                "MET",     "RB1",   "SETD2", "BRCA2",  "MGA",
                "GNAS")
gene_mut_rate_df = data.table(gene = subset_genes)
gene_mut_rate_df[,gene := toupper(gene)]

#' Calculate trinucleotide mutation proportions across the whole genome in preparation for calculating variant mutation rates

cesa = load_cesa(paste0(location_data,"pan_data_cesa.rds"))
smoking_samples = fread(paste0(location_data, 'smoking_sample_ids.txt'), header=F)[[1]]
nonsmoking_samples = fread(paste0(location_data, 'nonsmoking_sample_ids.txt'), header=F)[[1]]
panel_smoking_samples = fread(paste0(location_data, 'panel_smoking_sample_ids.txt'), header=F)[[1]]
panel_nonsmoking_samples = fread(paste0(location_data, 'panel_nonsmoking_sample_ids.txt'), header=F)[[1]]

ref_genome_version = 'hg19'

print('Calculating gene-level mutation rates from aggregate of mutation rates of individual variants')
genome_trinucleotide_mutation_proportions = genome_trinuc_mut_proportions(cesa$maf, ref_genome_version)
smoking_genome_trinucleotide_mutation_proportions = genome_trinuc_mut_proportions(cesa$maf[Unique_Patient_Identifier %in%
                                                                                             c(smoking_samples, panel_smoking_samples)], 
                                                                                    ref_genome_version)
nonsmoking_genome_trinucleotide_mutation_proportions = genome_trinuc_mut_proportions(cesa$maf[Unique_Patient_Identifier %in%
                                                                                             c(nonsmoking_samples, panel_nonsmoking_samples)],
                                                                                        ref_genome_version)

#' Store dataframe of all non-silent variants present within the gene across all samples (so that the same
#'   set of variants are used to calculate the gene mutation rate for all genes)
cesa_maf = copy(cesa$maf)
cesa_maf = cesa_maf[lengths(genes) == 1 & !is.na(genes) & is.na(top_gene), top_gene := as.character(genes)]
preloaded_maf = fread(paste0(location_data, 'all_preloaded_mafs.txt'))
# filtering out silent variants
variants_per_gene = merge(cesa_maf, preloaded_maf[is.na(problem) & Variant_Classification != 'Silent', 
                              .(Unique_Patient_Identifier, variant_id)], 
                          all=F)
# excluding all non-SNVs
variants_per_gene = variants_per_gene[top_gene %in% gene_mut_rate_df$gene & variant_type == 'snv' & !is.na(variant_id), .(top_gene, variant_id)]


#' Calculate variant level mutation rates
#' Then calculate the gene level mutation rate as the sum of all the variant level mutation rates within the gene

pan_data_gene_rates = fread(paste0(location_output, 'pan_data_mutation_rates.txt'))
gene_mut_rate_df[, pan_data := 
              sapply(gene, function(x)
                gene_mutation_rate(x,
                                   pan_data_gene_rates[gene==x, rate_grp_1],
                                   genome_trinucleotide_mutation_proportions,
                                   ref_genome_version,
                                   variants = variants_per_gene[top_gene == x, variant_id]))]
smoking_w_panel_gene_rates = fread(paste0(location_output, 'smoking_mutation_rates.txt'))
gene_mut_rate_df[, smoking := 
            sapply(gene, function(x) 
                gene_mutation_rate(x,
                                   smoking_w_panel_gene_rates[gene==x, rate_grp_1],
                                   smoking_genome_trinucleotide_mutation_proportions,
                                   ref_genome_version,
                                   variants = variants_per_gene[top_gene == x, variant_id]))]

nonsmoking_w_panel_gene_rates = fread(paste0(location_output, 'nonsmoking_mutation_rates.txt'))
gene_mut_rate_df[, nonsmoking := 
            sapply(gene, function(x) 
                gene_mutation_rate(x,
                                   nonsmoking_w_panel_gene_rates[gene==x, rate_grp_1],
                                   nonsmoking_genome_trinucleotide_mutation_proportions,
                                   ref_genome_version,
                                   variants = variants_per_gene[top_gene == x, variant_id]))]
print(unique(warnings()))

if(save_results){fwrite(gene_mut_rate_df, paste0(location_output,"variant_based_mutation_rates.txt"))}

#' Final MAF filtering steps prior to epistasis analysis
#' Output location: 'output/cesR_maf_for_epistasis_analysis.txt'
print("Filtering MAF")
maf_file = read.csv(paste0(location_output,"merged_luad_maf.txt"))
maf_file = maf_file %>% dplyr::rename('Tumor_Sample_Barcode' = 'Sample.ID',
                                      'Tumor_Allele' = 'Tumor_Seq_Allele2')

final_maf = construct_maf(cesa$maf, maf_file, preloaded_maf, variants_per_gene, save_results)

#' Additionally, create genes per sample table for new compute_samples functionality
#' Output location: 'output/genes_per_sample.txt'
print("Creating table of samples and mutations")
samples_genes = produce_genes_per_sample(final_maf, gene_mut_rate_df$gene, save_results)
