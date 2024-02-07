library(cancereffectsizeR)
library(ces.refset.hg19)
library(data.table)
library(dplyr)
library(rtracklayer)
library(stringr)


#' READ IN ALL MAF FILES
maf_file <- read.csv(paste0(location_output,"merged_luad_maf.txt"))
colnames(maf_file)[2] <- 'Tumor_Sample_Barcode'
colnames(maf_file)[7] <- 'Tumor_Allele'
maf_list <- split(maf_file, maf_file$Source) 

liftover_file = paste0(location_data, "hg38ToHg19.over.chain")

#' PRELOAD ALL MAF FILES
Broad_maf <- cancereffectsizeR::preload_maf(maf = maf_list$Broad, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
Broad_maf <- Broad_maf[is.na(problem)]
Broad_maf <- Broad_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]
Broad_maf <- Broad_maf[!Unique_Patient_Identifier %in% c("LUAD-B01169","LUAD-D01382")]

FMAD_maf <- cancereffectsizeR::preload_maf(maf = maf_list$`FM-AD`, ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
FMAD_maf <- FMAD_maf[is.na(problem)]

Genie_maf <- cancereffectsizeR::preload_maf(maf = maf_list$Genie, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
Genie_maf <- Genie_maf[is.na(problem)]

MSK2015_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2015, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
MSK2015_maf <- MSK2015_maf[is.na(problem)]
MSK2015_maf <- MSK2015_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

MSK2017_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2017, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
MSK2017_maf <- MSK2017_maf[is.na(problem)]

MSK2018_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2018, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
MSK2018_maf <- MSK2018_maf[is.na(problem)]

OncoSG_maf <- cancereffectsizeR::preload_maf(maf = maf_list$OncoSG, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
OncoSG_maf <- OncoSG_maf[is.na(problem)]
OncoSG_maf <- OncoSG_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TCGA_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TCGA, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
TCGA_maf <- TCGA_maf[is.na(problem)]
TCGA_maf <- TCGA_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TracerX_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TracerX, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
TracerX_maf <- TracerX_maf[is.na(problem)]
TracerX_maf <- TracerX_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TSP_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TSP, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
TSP_maf <- TSP_maf[is.na(problem)]
TSP_maf <- TSP_maf[!Unique_Patient_Identifier %in% c("luad_tsp_16929", "luad_tsp_16901", "luad_tsp_16875","luad_tsp_16915")]

NCI_maf <- cancereffectsizeR::preload_maf(maf = maf_list$NCI, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
NCI_maf <- NCI_maf[is.na(problem)]
NCI_maf <- NCI_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

#' EXTRACTING NECESSARY COVERAGE INFORMATION

#READING IN INFORMATION NECESSARY FOR LOADING MAF FILES (EXOME/GENOME AND COVERAGE INTERVALS FOR TGS)
broad_exome_or_genome <- fread(paste0(location_data,"luad_broad/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Platform')]
broad_exome_or_genome$Platform <- sapply(broad_exome_or_genome$Platform, function(x){if(str_detect(x, 'WGS')){return('WGS')} else return('WES')})
msk2017_panels_used <- fread(paste0(location_data,"lung_msk_2017/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Gene Panel')]
msk2018_panels_used <- fread(paste0(location_data,"nsclc_pd1_msk_2018/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Gene Panel')]
genie_panels_used <- fread(paste0(location_data,"genie_9/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Sequence Assay ID')]

#READING IN GENES INCLUDED IN EACH PANEL AND CREATING GRANGES OBJECT TO PASS INTO COVERED_REGIONS PARAMETER OF LOAD_MAF
#once the granges are exported once, these functions don't need to be run anymore

gene_granges <- rtracklayer::import(paste0(location_data, "gencode.v38lift37.basic.annotation.gtf"))

location_bed <- paste0(location_data,'bed_files/')
location_gene_panels <- paste0(location_data,'gene_panels/')

if(!file.exists(paste0(location_bed,"fmad_targets.bed"))){
  fmad_genes <- unique(fread(paste0(location_gene_panels,"foundation_one.txt"))$Hugo_Symbol)
  fmad_granges <- gene_granges[gene_granges$gene_name %in% fmad_genes, ]
  fmad_granges <- fmad_granges[fmad_granges$type %in% c('CDS','stop_codon'),]
  fmad_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = fmad_granges)
  export(fmad_gr_clean, paste0(location_bed,"fmad_targets.bed"))
}

if(!file.exists(paste0(location_bed,"tsp_targets.bed"))){
  tsp_genes <- unique(fread(paste0(location_gene_panels,"tsp.txt"))$Hugo_Symbol)
  tsp_granges <- gene_granges[gene_granges$gene_name %in% tsp_genes, ]
  tsp_granges <- tsp_granges[tsp_granges$type %in% c('CDS','stop_codon'),]
  tsp_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = tsp_granges)
  export(tsp_gr_clean, paste0(location_bed,"tsp_targets.bed"))
}

if(!file.exists(paste0(location_bed,"msk341_targets.bed"))){
  msk_341_genes <- unique(fread(paste0(location_gene_panels,"msk341.txt"))$Hugo_Symbol)
  msk341_granges <- gene_granges[gene_granges$gene_name %in% msk_341_genes, ]
  msk341_granges <- msk341_granges[msk341_granges$type %in% c('CDS','stop_codon'),]
  msk341_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = msk341_granges)
  export(msk341_gr_clean, paste0(location_bed,"msk341_targets.bed"))
}

if(!file.exists(paste0(location_bed,"msk410_targets.bed"))){
  msk_410_genes <- unique(fread(paste0(location_gene_panels,"msk410.txt"))$Hugo_Symbol)
  msk410_granges <- gene_granges[gene_granges$gene_name %in% msk_410_genes, ]
  msk410_granges <- msk410_granges[msk410_granges$type %in% c('CDS','stop_codon'),]
  msk410_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = msk410_granges)
  export(msk410_gr_clean, paste0(location_bed,"msk410_targets.bed"))
}

if(!file.exists(paste0(location_bed,"msk468_targets.bed"))){
  msk_468_genes <- unique(fread(paste0(location_gene_panels,"msk468.txt"))$Hugo_Symbol)
  msk468_granges <- gene_granges[gene_granges$gene_name %in% msk_468_genes, ]
  msk468_granges <- msk468_granges[msk468_granges$type %in% c('CDS','stop_codon'),]
  msk468_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = msk468_granges)
  export(msk468_gr_clean, paste0(location_bed,"msk468_targets.bed"))
}

#this should be run everytime so that we don't have to save a bunch of grange files
genie_panel_genes <- fread(paste0(location_data,"genie_9/genomic_information.txt"))[,c('Chromosome', 'Start_Position', 'End_Position', 'Hugo_Symbol', 'Feature_Type', 'SEQ_ASSAY_ID')]
genie_granges_list <- makeGRangesListFromDataFrame(genie_panel_genes, ignore.strand = T, seqnames.field = 'Chromosome', start.field = 'Start_Position', end.field = 'End_Position', split.field = 'SEQ_ASSAY_ID')
seqlevels(genie_granges_list, pruning.mode = "fine") <- c(1:22,'X','Y')

#SOME PANELS IN GENIE DON'T COVER TP53 OR KRAS SO THEY MUST BE REMOVED

genie_panel_genes_2 <- genie_panel_genes[,c('Hugo_Symbol','SEQ_ASSAY_ID')]
genie_panel_genes_2 <- genie_panel_genes_2[!duplicated(genie_panel_genes_2)]
genie_panel_genes_list <- split(genie_panel_genes_2$Hugo_Symbol, genie_panel_genes_2$SEQ_ASSAY_ID)

Genie_maf <- merge(Genie_maf, genie_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')

# panels_to_remove <- c()
# for(panel in names(genie_panel_genes_list)){
#   if(!('TP53' %in% genie_panel_genes_list[[panel]]) | !('KRAS' %in% genie_panel_genes_list[[panel]])){
#     panels_to_remove <- c(panels_to_remove, panel)
#   }
# }
# 
# Genie_maf <- Genie_maf[!'Sequence Assay ID' %in% panels_to_remove]

#SPLITTING MAF FILES WHEN DIFFERENT SEQUENCING METHODS ARE USED
Broad_maf <- merge(Broad_maf, broad_exome_or_genome, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
Broad_maf <- split(Broad_maf, Broad_maf$Platform)

MSK2017_maf <- merge(MSK2017_maf, msk2017_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
MSK2017_maf <- split(MSK2017_maf, MSK2017_maf$`Gene Panel`)

MSK2018_maf <- merge(MSK2018_maf, msk2018_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
MSK2018_maf <- split(MSK2018_maf, MSK2018_maf$`Gene Panel`)

Genie_maf <- merge(Genie_maf, genie_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
Genie_maf <- split(Genie_maf, Genie_maf$`Sequence Assay ID.x`)


preloaded_maf = rbind(rbindlist(Broad_maf), FMAD_maf, rbindlist(Genie_maf), MSK2015_maf, rbindlist(MSK2017_maf), rbindlist(MSK2018_maf), OncoSG_maf, TCGA_maf, TracerX_maf, TSP_maf, NCI_maf,fill = T)

#' LOADING IN MAF FILES INTO CESA OBJECT

cesa <- CESAnalysis(ces.refset.hg19)
#consider using covered regions padding when variants are outside the intervals
cesa <- load_maf(cesa, maf = Broad_maf$WES)
cesa <- load_maf(cesa, maf = MSK2015_maf)
cesa <- load_maf(cesa, maf = OncoSG_maf)
cesa <- load_maf(cesa, maf = TCGA_maf)
cesa <- load_maf(cesa, maf = TracerX_maf)
cesa <- load_maf(cesa, maf = FMAD_maf, coverage = 'targeted', covered_regions = paste0(location_bed,"fmad_targets.bed"), covered_regions_name = 'fmad_regions', covered_regions_padding = 100) #padding based on 23 variants having distance from interval between 10 and 100.

for(i in 1:length(Genie_maf)){
  cesa <- load_maf(cesa, maf = Genie_maf[i][[1]], coverage = 'targeted', covered_regions = genie_granges_list[names(Genie_maf)[i]][[1]], covered_regions_name = paste0(names(Genie_maf)[i], '_regions'), covered_regions_padding = 100)
}

cesa <- load_maf(cesa, maf = MSK2017_maf$IMPACT341, coverage = 'targeted', covered_regions = paste0(location_bed,"msk341_targets.bed"), covered_regions_name = 'msk341_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2017_maf$IMPACT410, coverage = 'targeted', covered_regions = paste0(location_bed,"msk410_targets.bed"), covered_regions_name = 'msk410_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2018_maf$IMPACT341, coverage = 'targeted', covered_regions = paste0(location_bed,"msk341_targets.bed"), covered_regions_name = 'msk341_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2018_maf$IMPACT410, coverage = 'targeted', covered_regions = paste0(location_bed,"msk410_targets.bed"), covered_regions_name = 'msk410_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2018_maf$IMPACT468, coverage = 'targeted', covered_regions = paste0(location_bed,"msk468_targets.bed"), covered_regions_name = 'msk468_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = TSP_maf, coverage = 'targeted', covered_regions = paste0(location_bed,"tsp_targets.bed"), covered_regions_name = 'tsp_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = Broad_maf$WGS, coverage = 'genome')
cesa <- load_maf(cesa, maf = NCI_maf, coverage = 'genome')


#' CALCULATING MUTATIONS RATES (NECESSARY FOR VARIANT_MUTATION_RATE)
cesa <- gene_mutation_rates(cesa, covariates = "lung")

#' CALCULATING MUTATION RATES FOR SMOKERS AND NONSMOKERS INDEPENDENTLY
#' First have to use mutational signature convolution on WES/WGS and clinical
#'   data for TGS to split the samples into smokers and nonsmokers
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "LUAD", treatment_naive = F)

#' Only using exome or genome samples, as mutational signatures can't be calculated for TGS
cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                              signature_exclusions = signature_exclusions, 
                              samples = cesa$samples[coverage %in% c('exome','genome'), Unique_Patient_Identifier],
                              cores=4
)

#SUBSETTING TO SAMPLES WITH UNBLENDED SIGNATURE WEIGHTS (SEE MESSAGES WITH JEFF MANDELL) AND WITH GREATER THAN 50 SNVS
bio_weights <- cesa$mutational_signatures$biological_weights
bio_weights_unblended <- bio_weights[bio_weights$group_avg_blended == F]
snv_counts <- cesa$maf[variant_type == 'snv', .N, by = "Unique_Patient_Identifier"]
good_samples <- snv_counts[N > 50, Unique_Patient_Identifier]
# NSLC_NCI patients will be added to the nonsmoking_samples list
NSLC_NCI_patients = unique(maf_list$NCI$Tumor_Sample_Barcode)
good_samples <- good_samples[! good_samples %in% NSLC_NCI_patients]
good_sample_weights <- bio_weights_unblended[Unique_Patient_Identifier %in% good_samples,]
#SMOKING SAMPLES ARE ANY SAMPLES WITH >0 SIGNATURE WEIGHTS.
smoking_samples <- good_sample_weights[SBS4 > 0, Unique_Patient_Identifier]
nonsmoking_samples <- good_sample_weights[SBS4 == 0, Unique_Patient_Identifier]

# We are confident that these patients are never-smokers, and the publication indicated that they had
# low smoking signature despite some having a history of secondary smoking.
nonsmoking_samples <- c(nonsmoking_samples,
                        NSLC_NCI_patients)


#INCLUDING PANEL DATA
maf_clinical = fread(paste0(location_output, 'merged_final.txt'))
panel_smoking_samples = unique(maf_clinical[Source %in% c('MSK2017','MSK2018')][Smoker == T, `Sample ID`])
panel_nonsmoking_samples = unique(maf_clinical[Source %in% c('MSK2017','MSK2018')][Smoker == F, `Sample ID`])

#CALCULATING MUTATION RATES FOR SMOKERS AND NONSMOKERS
## mutation rates of only WES/WGS
cesa_smoking = clear_gene_rates(cesa)
cesa_smoking = gene_mutation_rates(cesa_smoking, covariates = 'lung', 
                                   samples = cesa_smoking$samples[Unique_Patient_Identifier %in% smoking_samples, Unique_Patient_Identifier])

cesa_nonsmoking = clear_gene_rates(cesa)
cesa_nonsmoking = gene_mutation_rates(cesa_nonsmoking, covariates = 'lung', 
                                   samples = cesa_nonsmoking$samples[Unique_Patient_Identifier %in% nonsmoking_samples, Unique_Patient_Identifier])

# mutation rates of WES/WGS/TGS
cesa_smoking_w_panel = clear_gene_rates(cesa)
cesa_smoking_w_panel = gene_mutation_rates(cesa_smoking_w_panel, covariates = 'lung',
                                           samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])

cesa_nonsmoking_w_panel = clear_gene_rates(cesa)
cesa_nonsmoking_w_panel = gene_mutation_rates(cesa_nonsmoking_w_panel, covariates = 'lung',
                                              samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])


if(save_results){
  save_cesa(cesa, paste0(location_data,"pan_data_cesa_for_cancer_epistasis.rds"))
  
  fwrite(preloaded_maf, paste0(location_data, 'all_preloaded_mafs.txt'))
  
  fwrite(list(smoking_samples), paste0(location_data, 'smoking_sample_ids.txt'))
  fwrite(list(nonsmoking_samples), paste0(location_data, 'nonsmoking_sample_ids.txt'))
  
  fwrite(list(panel_smoking_samples), paste0(location_data, 'panel_smoking_sample_ids.txt'))
  fwrite(list(panel_nonsmoking_samples), paste0(location_data, 'panel_nonsmoking_sample_ids.txt'))
  
  pan_data_gene_rates = cesa$gene_rates
  colnames(pan_data_gene_rates)[2] <- "rate_grp_1"
  fwrite(pan_data_gene_rates, paste0(location_output, 'pan_data_mutation_rates.txt'))
  fwrite(cesa_smoking$gene_rates, paste0(location_output, 'smoking_mutation_rates.txt'))
  fwrite(cesa_nonsmoking$gene_rates, paste0(location_output, 'nonsmoking_mutation_rates.txt'))
  fwrite(cesa_smoking_w_panel$gene_rates, paste0(location_output, 'smoking_w_panel_mutation_rates.txt'))
  fwrite(cesa_nonsmoking_w_panel$gene_rates, paste0(location_output, 'nonsmoking_w_panel_mutation_rates.txt'))
}
