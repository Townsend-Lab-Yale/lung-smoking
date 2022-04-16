#' Creates a maf table of non-silent snv's with some variants removed via cesR filtering. To be used as the maf file for the cancer epistasis analysis.

#cesa constructed as in mutation_rates.R
cesa = load_cesa('../../data/pan_data_cesa_for_cancer_epistasis.rds')
cesa@maf = cesa$maf[lengths(genes) == 1 & is.na(top_gene), top_gene := genes]

preload_maf_variants = rbind(Broad_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              FMAD_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              Genie_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              MSK2015_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              MSK2017_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              MSK2018_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              OncoSG_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              TCGA_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              TracerX_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                              TSP_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)])
cesa_maf_variants = cesa$maf[variant_type == 'snv', .(Unique_Patient_Identifier, variant_id)]
variants_to_keep = merge(cesa_maf_variants, preload_maf_variants) 

filtered_maf = merge(cesa$maf, variants_to_keep, by = c('Unique_Patient_Identifier','variant_id'))
filtered_maf = filtered_maf[,-c('prelift_chr', 'prelift_start', 'liftover_strand_flip','variant_type', 'genes', 'top_consequence')]

source_panel_info = as.data.table(maf_file)
source_panel_info = source_panel_info[!duplicated(Tumor_Sample_Barcode),.(Tumor_Sample_Barcode, Source, Panel)]
source_panel_info[Panel == '', Panel := NA]

filtered_maf = merge(filtered_maf, source_panel_info, by.x = 'Unique_Patient_Identifier', by.y = 'Tumor_Sample_Barcode', all.x = T)

colnames(filtered_maf) = c('Sample ID', 'Mutation', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Gene', 'Source', 'Panel')

fwrite(filtered_maf, '../output/cesR_maf_for_epistasis_analysis.txt')
