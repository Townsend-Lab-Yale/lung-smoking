#' Creates a maf table of non-silent snv's with some variants removed via cesR filtering. To be used as the maf file for the cancer epistasis analysis.

library(cancereffectsizeR)
library(data.table)


#' Assigning a single gene to each variant
cesa@maf = cesa$maf[lengths(genes) == 1 & is.na(top_gene), top_gene := genes]

#' Generating a list of variants when the maf is preloaded and when the maf is loaded
#' Taking the intersection of these variants as the ones to keep
preload_maf_variants = rbind(Broad_maf$WES[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             Broad_maf$WGS[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             FMAD_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             rbindlist(Genie_maf)[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             MSK2015_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             MSK2017_maf$IMPACT341[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             MSK2017_maf$IMPACT410[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             MSK2018_maf$IMPACT341[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             MSK2018_maf$IMPACT410[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],  
                             MSK2018_maf$IMPACT468[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             OncoSG_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             TCGA_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             TracerX_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)],
                             TSP_maf[variant_type == 'snv' & is.na(problem) & Variant_Classification != 'Silent', .(Unique_Patient_Identifier, variant_id)])
cesa_maf_variants = cesa$maf[variant_type == 'snv', .(Unique_Patient_Identifier, variant_id)]
variants_to_keep = merge(cesa_maf_variants, preload_maf_variants) 

#' Keeping desired variants and removing extraneous columns
filtered_maf = merge(cesa$maf, variants_to_keep, by = c('Unique_Patient_Identifier','variant_id'))
filtered_maf = filtered_maf[,-c('prelift_chr', 'prelift_start', 'liftover_strand_flip','variant_type', 'genes', 'top_consequence')]

#' Adding source information to all of the TGS samples
source_panel_info = as.data.table(maf_file)
source_panel_info = source_panel_info[!duplicated(Tumor_Sample_Barcode),.(Tumor_Sample_Barcode, Source, Panel)]
source_panel_info[Panel == '', Panel := NA]

filtered_maf = merge(filtered_maf, source_panel_info, by.x = 'Unique_Patient_Identifier', by.y = 'Tumor_Sample_Barcode', all.x = T)

#' Renaming columns for consistency between Python and R
colnames(filtered_maf) = c('Sample ID', 'Mutation', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Gene', 'Source', 'Panel')

#' Output
fwrite(filtered_maf, paste0(location_data, "cesR_maf_for_epistasis_analysis.txt"))