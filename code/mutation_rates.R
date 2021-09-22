library(cancereffectsizeR)
library(ces.refset.hg19)
library(data.table)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(stringr)

#READ IN ALL MAF FILES
maf_file <- read.csv("../output/merged_luad_maf.txt")
colnames(maf_file)[2] <- 'Tumor_Sample_Barcode'
colnames(maf_file)[7] <- 'Tumor_Allele'
maf_list <- split(maf_file, maf_file$Source)

#couldn't make this into a loop or an sapply type function so for now here's a code wall
Broad_maf <- cancereffectsizeR::preload_maf(maf = maf_list$Broad, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
Broad_maf <- Broad_maf[is.na(problem)]
Broad_maf <- Broad_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]
Broad_maf <- Broad_maf[!Unique_Patient_Identifier %in% c("LUAD-B01169","LUAD-D01382")]

FMAD_maf <- cancereffectsizeR::preload_maf(maf = maf_list$`FM-AD`, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
FMAD_maf <- FMAD_maf[is.na(problem)]
FMAD_maf <- FMAD_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

Genie_maf <- cancereffectsizeR::preload_maf(maf = maf_list$Genie, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
Genie_maf <- Genie_maf[is.na(problem)]
Genie_maf <- Genie_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

MSK2015_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2015, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
MSK2015_maf <- MSK2015_maf[is.na(problem)]
MSK2015_maf <- MSK2015_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

MSK2017_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2017, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
MSK2017_maf <- MSK2017_maf[is.na(problem)]
MSK2017_maf <- MSK2017_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

MSK2018_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2018, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
MSK2018_maf <- MSK2018_maf[is.na(problem)]
MSK2018_maf <- MSK2018_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

OncoSG_maf <- cancereffectsizeR::preload_maf(maf = maf_list$OncoSG, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
OncoSG_maf <- OncoSG_maf[is.na(problem)]
OncoSG_maf <- OncoSG_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TCGA_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TCGA, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
TCGA_maf <- TCGA_maf[is.na(problem)]
TCGA_maf <- TCGA_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TracerX_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TracerX, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
TracerX_maf <- TracerX_maf[is.na(problem)]
TracerX_maf <- TracerX_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TSP_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TSP, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")
TSP_maf <- TSP_maf[is.na(problem)]
TSP_maf <- TSP_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]
TSP_maf <- TSP_maf[!Unique_Patient_Identifier %in% c("luad_tsp_16929", "luad_tsp_16901", "luad_tsp_16875","luad_tsp_16915")]

#maf_list <- sapply(maf_list, function(x){preload_maf(maf = x, refset = ces.refset.hg19, chain_file = "../data/hg38ToHg19.over.chain")})


#SIGNATURES NOT RELEVANT TO LUAD
signatures_to_remove <- suggest_cosmic_signatures_to_remove(cancer_type = "LUAD", treatment_naive = F)


#READING IN INFORMATION NECESSARY FOR LOADING MAF FILES (EXOME/GENOME AND COVERAGE INTERVALS FOR TGS)
broad_exome_or_genome <- fread('../data/luad_broad/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Platform')]
broad_exome_or_genome$Platform <- sapply(broad_exome_or_genome$Platform, function(x){if(str_detect(x, 'WGS')){return('WGS')} else return('WES')})
msk2017_panels_used <- fread('../data/lung_msk_2017/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Gene Panel')]
msk2018_panels_used <- fread('../data/nsclc_pd1_msk_2018/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Gene Panel')]
genie_panels_used <- fread('../data/genie_9/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Sequence Assay ID')]

tsp_genes <- fread('../data/gene_panels/tsp.txt')
foundation_one_genes <- fread('../data/gene_panels/foundation_one.txt')
msk341_genes <- fread('../data/gene_panels/msk341.txt')
msk410_genes <- fread('../data/gene_panels/msk410.txt')
msk468_genes <- fread('../data/gene_panels/msk468.txt')
genie_genes <- fread('../data/gene_panels/genie_panel_genes.txt')


#READING IN GENES INCLUDED IN EACH PANEL AND CREATING GRANGES OBJECT TO PASS INTO COVERED_REGIONS PARAMETER OF LOAD_MAF
#once the granges are exported once, these functions don't need to be run anymore

gene_granges <- rtracklayer::import('../data/gencode.v38lift37.basic.annotation.gtf')

if(!file.exists('../data/fmad_targets.bed')){
  fmad_genes <- unique(fread('../gene_panels/fm-ad_genes.txt', sep = '\n', header = F))$V1
  fmad_granges <- gene_granges[gene_granges$gene_name %in% fmad_genes, ]
  fmad_granges <- fmad_granges[fmad_granges$type %in% c('CDS','stop_codon'),]
  fmad_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(), gr = fmad_granges)
  export(fmad_gr_clean, '../data/fmad_targets.bed')
} 

if(!file.exists('../data/tsp_targets.bed')){
  tsp_genes <- unique(fread('../gene_panels/tsp.txt', sep = '\n', header = F))$V1
  tsp_granges <- gene_granges[gene_granges$gene_name %in% tsp_genes, ]
  tsp_granges <- tsp_granges[tsp_granges$type %in% c('CDS','stop_codon'),]
  tsp_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(), gr = tsp_granges)
  export(tsp_gr_clean, '../data/tsp_targets.bed')
}

if(!file.exists('../data/msk341_targets.bed')){
  msk_341_genes <- unique(fread('../gene_panels/msk341.txt', sep = '\n', header = F)[V1 != 'nan'])$V1
  msk341_granges <- gene_granges[gene_granges$gene_name %in% msk_341_genes, ]
  msk341_granges <- msk341_granges[msk341_granges$type %in% c('CDS','stop_codon'),]
  msk341_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(), gr = msk341_granges)
  export(msk341_gr_clean, '../data/msk341_targets.bed')
}

if(!file.exists('../data/msk410_targets.bed')){
  msk_410_genes <- unique(fread('../gene_panels/msk410.txt', sep = '\n', header = F)[V1 != 'nan'])$V1
  msk410_granges <- gene_granges[gene_granges$gene_name %in% msk_410_genes, ]
  msk410_granges <- msk410_granges[msk410_granges$type %in% c('CDS','stop_codon'),]
  msk410_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(), gr = msk410_granges)
  export(msk410_gr_clean, '../data/msk410_targets.bed')
}

if(!file.exists('../data/msk468_targets.bed')){
  msk_468_genes <- unique(fread('../gene_panels/msk468.txt', sep = '\n', header = F)[V1 != 'nan'])$V1
  msk468_granges <- gene_granges[gene_granges$gene_name %in% msk_468_genes, ]
  msk468_granges <- msk468_granges[msk468_granges$type %in% c('CDS','stop_codon'),]
  msk468_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(), gr = msk468_granges)
  export(msk468_gr_clean, '../data/msk468_targets.bed')
}

#this should be run everytime so that we don't have to save a bunch of grange files
genie_panel_genes <- fread('../data/genie_9/genomic_information.txt')[,c('Chromosome', 'Start_Position', 'End_Position', 'Hugo_Symbol', 'Feature_Type', 'SEQ_ASSAY_ID')]
genie_granges_list <- makeGRangesListFromDataFrame(genie_panel_genes, ignore.strand = T, seqnames.field = 'Chromosome', start.field = 'Start_Position', end.field = 'End_Position', split.field = 'SEQ_ASSAY_ID')
seqlevels(genie_granges_list, pruning.mode = "fine") <- c(1:22,'X','Y')

#SOME PANELS IN GENIE DON'T COVER TP53 OR KRAS SO THEY MUST BE REMOVED

genie_panel_genes_2 <- genie_panel_genes[,c('Hugo_Symbol','SEQ_ASSAY_ID')]
genie_panel_genes_2 <- genie_panel_genes_2[!duplicated(genie_panel_genes_2)]
genie_panel_genes_list <- split(genie_panel_genes_2$Hugo_Symbol, genie_panel_genes_2$SEQ_ASSAY_ID)

Genie_maf <- merge(Genie_maf, genie_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')

panels_to_remove <- c()
for(panel in names(genie_panel_genes_list)){
  if(!('TP53' %in% genie_panel_genes_list[[panel]]) | !('KRAS' %in% genie_panel_genes_list[[panel]])){
    panels_to_remove <- c(panels_to_remove, panel)
  }
}

Genie_maf <- Genie_maf[!`Sequence Assay ID` %in% panels_to_remove]

#LIST OF ALL GENES INCLUDED IN ANALYSIS
genes_list <- fread('../data/genes_list.txt', header = F)$V1


#SPLITTING MAF FILES WHEN DIFFERENT SEQUENCING METHODS ARE USED
Broad_maf <- merge(Broad_maf, broad_exome_or_genome, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
Broad_maf <- split(Broad_maf, Broad_maf$Platform)

MSK2017_maf <- merge(MSK2017_maf, msk2017_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
MSK2017_maf <- split(MSK2017_maf, MSK2017_maf$`Gene Panel`)

MSK2018_maf <- merge(MSK2018_maf, msk2018_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
MSK2018_maf <- split(MSK2018_maf, MSK2018_maf$`Gene Panel`)

Genie_maf <- merge(Genie_maf, genie_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
Genie_maf <- split(Genie_maf, Genie_maf$`Sequence Assay ID`)

#CREATING PAN-DATASET CESA OBJECT FOR GENERAL MUTATION RATES
cesa_total <- CESAnalysis()
#consider using covered regions padding when variants are outside the intervals
cesa_total <- load_maf(cesa_total, maf = Broad_maf$WES)
cesa_total <- load_maf(cesa_total, maf = MSK2015_maf)
cesa_total <- load_maf(cesa_total, maf = OncoSG_maf)
cesa_total <- load_maf(cesa_total, maf = TCGA_maf)
cesa_total <- load_maf(cesa_total, maf = TracerX_maf)
cesa_total <- load_maf(cesa_total, maf = FMAD_maf, coverage = 'targeted', covered_regions = '../data/fmad_targets.bed', covered_regions_name = 'fmad_regions', covered_regions_padding = 100) #padding based on 23 variants having distance from interval between 10 and 100.

for(i in 12:length(Genie_maf)){
  cesa_total <- load_maf(cesa_total, maf = Genie_maf[i][[1]], coverage = 'targeted', covered_regions = genie_granges_list[names(Genie_maf)[i]][[1]], covered_regions_name = paste0(names(Genie_maf)[i], '_regions'), covered_regions_padding = 100)
}

cesa_total <- load_maf(cesa_total, maf = MSK2017_maf$IMPACT341, coverage = 'targeted', covered_regions = '../data/msk341_targets.bed', covered_regions_name = 'msk341_regions', covered_regions_padding = 100)
cesa_total <- load_maf(cesa_total, maf = MSK2017_maf$IMPACT410, coverage = 'targeted', covered_regions = '../data/msk410_targets.bed', covered_regions_name = 'msk410_regions', covered_regions_padding = 100)
cesa_total <- load_maf(cesa_total, maf = MSK2018_maf$IMPACT341, coverage = 'targeted', covered_regions = '../data/msk341_targets.bed', covered_regions_name = 'msk341_regions', covered_regions_padding = 100)
cesa_total <- load_maf(cesa_total, maf = MSK2018_maf$IMPACT410, coverage = 'targeted', covered_regions = '../data/msk410_targets.bed', covered_regions_name = 'msk410_regions', covered_regions_padding = 100)
cesa_total <- load_maf(cesa_total, maf = MSK2018_maf$IMPACT468, coverage = 'targeted', covered_regions = '../data/msk468_targets.bed', covered_regions_name = 'msk468_regions', covered_regions_padding = 100)
cesa_total <- load_maf(cesa_total, maf = TSP_maf, coverage = 'targeted', covered_regions = '../data/tsp_targets.bed', covered_regions_name = 'tsp_regions', covered_regions_padding = 100)
cesa_total <- load_maf(cesa_total, maf = Broad_maf$WGS, coverage = 'genome')


#CALCULATING MUTATION RATES
cesa_total <- trinuc_mutation_rates(cesa_total,
                                    signature_set = "COSMIC_v3.1",
                                    signatures_to_remove = signatures_to_remove
)
cesa_total <- gene_mutation_rates(cesa_total, covariates = "lung")

save_cesa(cesa_total, '../data/pan-dataset_samples_cesa.rds')
fwrite(cesa_total$gene_rates, '../data/pan-data_mutation_rates.txt')


#WES/WGS-ONLY CESA ANALYSIS BECAUSE ONLY THESE CAN HAVE SIGNATURE EXTRACTIONS BE PERFORMED ON THEM
#cesa_exome <- CESAnalysis()
#cesa_exome <- load_maf(cesa_exome, maf = Broad_maf$WGS, coverage = 'genome')
#cesa_exome <- load_maf(cesa_exome, maf = Broad_maf$WES)
#cesa_exome <- load_maf(cesa_exome, maf = MSK2015_maf)
#cesa_exome <- load_maf(cesa_exome, maf = OncoSG_maf)
#cesa_exome <- load_maf(cesa_exome, maf = TCGA_maf)
#cesa_exome <- load_maf(cesa_exome, maf = TracerX_maf)

#cesa_exome <- trinuc_mutation_rates(cesa_exome,
#                                    signature_set = "COSMIC_v3.1",
#                                    signatures_to_remove = signatures_to_remove
#)


#CALCULATING MUTATION RATES FOR COMPARISON PURPOSES WITH THE PAN-DATASET MUTATION RATES
#cesa_exome <- gene_mutation_rates(cesa_exome, covariates = "lung")
#save_cesa(cesa_exome, '../data/exome_samples_cesa.rds')
#fwrite(cesa_exome$gene_rates, '../data/exome_mutation_rates.txt')


#SUBSETTING TO SAMPLES WITH UNBLENDED SIGNATURE WEIGHTS (SEE MESSAGES WITH JEFF MANDELL) AND WITH GREATER THAN 50 SNVS
bio_weights <- cesa_exome$mutational_signatures$biological_weights
bio_weights_unblended <- bio_weights[bio_weights$group_avg_blended == F]
snv_counts <- cesa_exome$maf[variant_type == 'snv', .N, by = "Unique_Patient_Identifier"]
good_samples <- snv_counts[N > 50, Unique_Patient_Identifier]
good_sample_weights <- bio_weights_unblended[Unique_Patient_Identifier %in% good_samples,]
#SMOKING SAMPLES ARE ANY SAMPLES WITH >0 SIGNATURE WEIGHTS.
smoking_samples <- good_sample_weights[SBS4 > 0, Unique_Patient_Identifier]
nonsmoking_samples <- good_sample_weights[SBS4 == 0, Unique_Patient_Identifier]


#SMOKING GROUP
cesa_smoking <- CESAnalysis()
cesa_smoking <- load_maf(cesa_smoking, maf = Broad_maf$WES[Broad_maf$WES$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = MSK2015_maf[MSK2015_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = OncoSG_maf[OncoSG_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = TCGA_maf[TCGA_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = TracerX_maf[TracerX_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = Broad_maf$WGS[Broad_maf$WGS$Unique_Patient_Identifier %in% smoking_samples], coverage = 'genome')

cesa_smoking <- trinuc_mutation_rates(cesa_smoking,
                                      signature_set = "COSMIC_v3.1",
                                      signatures_to_remove = signatures_to_remove
)

cesa_smoking <- gene_mutation_rates(cesa_smoking, covariates = "lung")

save_cesa(cesa_smoking, '../data/smoking_samples_cesa.rds')
fwrite(cesa_smoking$gene_rates, '../data/smoking_mutation_rates.txt')


#NONSMOKING GROUP
cesa_nonsmoking <- CESAnalysis()
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = Broad_maf$WES[Broad_maf$WES$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = MSK2015_maf[MSK2015_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = OncoSG_maf[OncoSG_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = TCGA_maf[TCGA_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = TracerX_maf[TracerX_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = Broad_maf$WGS[Broad_maf$WGS$Unique_Patient_Identifier %in% nonsmoking_samples], coverage = 'genome')

cesa_nonsmoking <- trinuc_mutation_rates(cesa_nonsmoking,
                                      signature_set = "COSMIC_v3.1",
                                      signatures_to_remove = signatures_to_remove
)
cesa_nonsmoking <- gene_mutation_rates(cesa_nonsmoking, covariates = "lung")
save_cesa(cesa_nonsmoking, '../data/nonsmoking_samples_cesa.rds')
fwrite(cesa_nonsmoking$gene_rates, '../data/nonsmoking_mutation_rates.txt')