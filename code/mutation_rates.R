library(cancereffectsizeR)
library(ces.refset.hg19)
library(data.table)
library(dplyr)
library(ggplot2)
library(rtracklayer)
# All files not available in repo but still used here:
# Liftover,
# Gencode basic annotation,
# Genie genomic information

#READ IN ALL MAF FILES
maf_file <- read.csv("../output/merged_luad_maf.txt")
colnames(maf_file)[2] <- 'Tumor_Sample_Barcode'
colnames(maf_file)[7] <- 'Tumor_Allele'


#CHECKING COVERAGE INTERVALS TO SEE IF PADDING IS NECESSARY
temp_maf_list <- split(maf_file, maf_file$Source)
TSP_maf <- preload_maf(maf = temp_maf_list$TSP, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain", coverage_intervals_to_check = '../data/tsp_targets.bed')
TSP_maf[TSP_maf$dist_to_coverage_intervals %in% 100:500]


#PRELOADING ALL MAFS AT ONCE
#may need to preload each dataset individually
LUAD_maf <- cancereffectsizeR::preload_maf(maf = maf_file, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain", keep_extra_columns = c('Source'))


#REMOVING POSSIBLE PROBLEMATIC SAMPLES
LUAD_maf <- LUAD_maf[is.na(problem)]
LUAD_maf <- LUAD_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]
samples_to_remove <- c("luad_tsp_16929", "luad_tsp_16901", "luad_tsp_16875","luad_tsp_16915","LUAD-B01169","LUAD-D01382")
LUAD_maf <- LUAD_maf[!Unique_Patient_Identifier %in% samples_to_remove]


#SIGNATURES NOT RELEVANT TO LUAD
signatures_to_remove <- suggest_cosmic_signatures_to_remove(cancer_type = "LUAD", treatment_naive = F)

maf_list <- split(LUAD_maf, LUAD_maf$Source)

#READING IN INFORMATION NECESSARY FOR LOADING IN INDIVIDUAL MAF FILES (EXOME/GENOME AND COVERAGE INTERVALS FOR TGS)
broad_exome_or_genome <- fread('../data/luad_broad/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Platform')]
broad_exome_or_genome$Platform <- sapply(broad_exome_or_genome$Platform, function(x){if(str_detect(x, 'WGS')){return('WGS')} else return('WES')})
msk2017_panels_used <- fread('../data/lung_msk_2017/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Gene Panel')]
msk2018_panels_used <- fread('../data/nsclc_pd1_msk_2018/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Gene Panel')]
genie_panels_used <- fread('../data/genie_9/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Sequence Assay ID')]


#READING IN GENES INCLUDED IN EACH PANEL AND CREATING GRANGES OBJECT TO PASS INTO COVERED_REGIONS PARAMETER OF LOAD_MAF
gene_granges <- rtracklayer::import('/Users/Krishna1/Downloads/gencode.v38lift37.basic.annotation.gtf')
fmad_genes <- unique(fread('../gene_panels/fm-ad_genes.txt', sep = '\n', header = F))$V1
msk_341_genes <- unique(fread('../gene_panels/msk341.txt', sep = '\n', header = F)[V1 != 'nan'])$V1
msk_410_genes <- unique(fread('../gene_panels/msk410.txt', sep = '\n', header = F)[V1 != 'nan'])$V1
msk_468_genes <- unique(fread('../gene_panels/msk468.txt', sep = '\n', header = F)[V1 != 'nan'])$V1
tsp_genes <- unique(fread('../gene_panels/tsp.txt', sep = '\n', header = F))$V1

fmad_granges <- gene_granges[gene_granges$gene_name %in% fmad_genes, ]
fmad_granges <- fmad_granges[fmad_granges$type %in% c('CDS','stop_codon'),]
fmad_gr_clean = cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(), gr = fmad_granges)
  
tsp_granges <- gene_granges[gene_granges$gene_name %in% tsp_genes, ]
tsp_granges <- tsp_granges[tsp_granges$type %in% c('CDS','stop_codon'),]
tsp_gr_clean = cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(), gr = tsp_granges)

msk341_granges <- gene_granges[gene_granges$gene_name %in% msk_341_genes, ]
msk341_granges <- msk341_granges[msk341_granges$type %in% c('CDS','stop_codon'),]

msk410_granges <- gene_granges[gene_granges$gene_name %in% msk_410_genes, ]
msk410_granges <- msk410_granges[msk410_granges$type %in% c('CDS','stop_codon'),]

msk468_granges <- gene_granges[gene_granges$gene_name %in% msk_468_genes, ]
msk468_granges <- msk468_granges[msk468_granges$type %in% c('CDS','stop_codon'),]

genie_panel_genes <- fread('../data/genie_9/genomic_information.txt')[,c('Chromosome', 'Start_Position', 'End_Position', 'Hugo_Symbol', 'Feature_Type', 'SEQ_ASSAY_ID')]
genie_granges_list <- makeGRangesListFromDataFrame(genie_panel_genes, keep.extra.columns = T, ignore.strand = T, seqnames.field = 'Chromosome', start.field = 'Start_Position', end.field = 'End_Position', split.field = 'SEQ_ASSAY_ID')


#GENES INCLUDED IN EVERY PANEL (NOT TRUE FOR ALL THE PANELS IN GENIE THOUGH, SEE MESSAGES WITH JORGE)
genes_in_intersection <- c('BRAF', 'ERBB4', 'NF2', 'ETV6', 'BRCA1', 'TSHR', 'MPL', 'MSH2', 'PIK3R2', 'JUN', 'PDGFRA', 'MAP2K1', 'ABL1', 'CCND3', 'FGFR4', 'ATM', 'PTEN', 'SUFU', 'GNAS', 'NTRK1', 'KIT', 'BCL2', 'FLT1', 'STK11', 'CDKN1B', 'NRAS', 'DDR2', 'CCND1', 'JAK1', 'ERG', 'JAK2', 'MYC', 'AKT2', 'PIK3C2G', 'MSH6', 'CHEK1', 'NTRK3', 'ALK', 'NOTCH4', 'GSK3B', 'RAF1', 'EPHA5', 'TP53', 'NTRK2', 'FLT4', 'NOTCH3', 'PIK3C3', 'CDK8', 'PIK3R1', 'MEN1', 'EGFR', 'IGF1R', 'FGFR2', 'CTNNB1', 'IKBKE', 'CSF1R', 'BAP1', 'CDH1', 'CDKN2B', 'KRAS', 'AKT1', 'FGF3', 'HRAS', 'PTPN11', 'MAP2K2', 'BTK', 'MDM2', 'BCL6', 'CHEK2', 'NF1', 'FGF4', 'CDKN2C', 'CCNE1', 'KDR', 'RUNX1', 'CRKL', 'EPHB1', 'BARD1', 'BRCA2', 'MDM4', 'RB1', 'SRC', 'FBXW7', 'FGFR1', 'TSC2', 'FGFR3', 'ERBB3', 'CDKN2A', 'MAP3K1', 'PRKAR1A', 'REL', 'PIK3CA', 'MLH1', 'SMAD4', 'AXL', 'CDK6', 'MAP3K13', 'ATR', 'TSC1', 'VHL', 'RET', 'AURKB', 'APC', 'MAP2K4', 'CDK4', 'CCND2', 'MYCN', 'PIK3CG', 'SMAD2', 'MET', 'WT1', 'NOTCH1', 'JAK3', 'NOTCH2', 'TGFBR2', 'ARAF', 'EPHA3', 'CBL', 'IRS2', 'ERBB2', 'PDGFRB', 'SMO', 'GATA1', 'AKT3', 'FLT3', 'SYK', 'RICTOR')


#CREATING PAN-DATASET CESA OBJECT FOR GENERAL MUTATION RATES
# cesa_total <- CESAnalysis()
# #consider using covered regions padding when variants are outside the intervals
# cesa_total <- load_maf(cesa_total, maf = maf_list$Broad)
# cesa_total <- load_maf(cesa_total, maf = maf_list$MSK2015)
# cesa_total <- load_maf(cesa_total, maf = maf_list$OncoSG)
# cesa_total <- load_maf(cesa_total, maf = maf_list$TCGA)
# cesa_total <- load_maf(cesa_total, maf = maf_list$TracerX)
# cesa_total <- load_maf(cesa_total, maf = maf_list$FM-AD, coverage = 'targeted', covered_regions = '../data/fmad_targets.bed', covered_regions_name = 'fmad_regions', covered_regions_padding = 100) #padding based on 23 variants having distance from interval between 10 and 100.
# cesa_total <- load_maf(cesa_total, maf = maf_list$Genie, coverage = 'targeted')
# cesa_total <- load_maf(cesa_total, maf = maf_list$MSK2017, coverage = 'targeted')
# cesa_total <- load_maf(cesa_total, maf = maf_list$MSK2018, coverage = 'targeted')
# cesa_total <- load_maf(cesa_total, maf = maf_list$TSP, coverage = 'targeted', covered_regions = '../data/tsp_targets.bed', covered_regions_name = 'tsp_regions', covered_regions_padding = 100)
# 
# cesa_total <- trinuc_mutation_rates(cesa_total,
#                                     signature_set = "COSMIC_v3.2",
#                                     signatures_to_remove = signatures_to_remove
# )
# cesa_total <- gene_mutation_rates(cesa_total, covariates = "lung")
# 
# save_cesa(cesa_total, 'pan-dataset_samples_cesa.rds')
# fwrite(cesa_total$gene_rates[gene %in% genes_in_intersection], '../data/pan-data_mutation_rates.txt')


#WES/WGS-ONLY CESA ANALYSIS BECAUSE ONLY THESE CAN HAVE SIGNATURE EXTRACTIONS BE PERFORMED ON THEM
cesa_exome <- CESAnalysis()
cesa_exome <- load_maf(cesa_exome, maf = maf_list$Broad)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$MSK2015)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$OncoSG)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$TCGA)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$TracerX)

cesa_exome <- trinuc_mutation_rates(cesa_exome,
                              signature_set = "COSMIC_v3.2",
                              signatures_to_remove = signatures_to_remove
)


#CALCULATING MUTATION RATES FOR COMPARISON PURPOSES WITH THE PAN-DATASET MUTATION RATES
cesa_exome <- gene_mutation_rates(cesa_exome, covariates = "lung")
save_cesa(cesa_exome, 'exome_samples_cesa.rds')
fwrite(cesa_exome$gene_rates[gene %in% genes_in_intersection], '../data/exome_mutation_rates.txt')


#SUBSETTING TO SAMPLES WITH UNBLENDED SIGNATURE WEIGHTS (SEE MESSAGES WITH JEFF MANDELL) AND WITH GREATER THAN 50 SNVS
bio_weights <- cesa_exome$mutational_signatures$biological_weights
bio_weights_unblended <- bio_weights[bio_weights$group_avg_blended == F]
snv_counts <- cesa_exome$maf[variant_type == 'snv', .N, by = "Unique_Patient_Identifier"]
good_samples <- snv_counts[N > 50, Unique_Patient_Identifier]
good_sample_weights <- bio_weights_unblended[Unique_Patient_Identifier %in% good_samples,]
#SMOKING SAMPLES ARE ANY SAMPLES WITH >0 SIGNATURE WEIGHTS, A HEURISTIC.
smoking_samples <- good_sample_weights[SBS4 > 0, Unique_Patient_Identifier]
nonsmoking_samples <- good_sample_weights[SBS4 == 0, Unique_Patient_Identifier]


#SMOKING GROUP
cesa_smoking <- CESAnalysis()
maf_smoking <- data.table()
maf_smoking <- rbind(maf_smoking, maf_list$Broad[maf_list$Broad$Unique_Patient_Identifier %in% smoking_samples])
maf_smoking <- rbind(maf_smoking, maf_list$MSK2015[maf_list$MSK2015$Unique_Patient_Identifier %in% smoking_samples])
maf_smoking <- rbind(maf_smoking, maf_list$OncoSG[maf_list$OncoSG$Unique_Patient_Identifier %in% smoking_samples])
maf_smoking <- rbind(maf_smoking, maf_list$TCGA[maf_list$TCGA$Unique_Patient_Identifier %in% smoking_samples])
maf_smoking <- rbind(maf_smoking, maf_list$TracerX[maf_list$TracerX$Unique_Patient_Identifier %in% smoking_samples])

cesa_smoking <- load_maf(cesa_smoking, maf = maf_smoking)

cesa_smoking <- trinuc_mutation_rates(cesa_smoking,
                                    signature_set = "COSMIC_v3.2",
                                    signatures_to_remove = signatures_to_remove
)

cesa_smoking <- gene_mutation_rates(cesa_smoking, covariates = "lung")

save_cesa(cesa_smoking, 'smoking_samples_cesa.rds')


#NONSMOKING GROUP
cesa_nonsmoking <- CESAnalysis()
maf_nonsmoking <- data.table()
maf_nonsmoking <- rbind(maf_nonsmoking, maf_list$Broad[maf_list$Broad$Unique_Patient_Identifier %in% nonsmoking_samples])
maf_nonsmoking <- rbind(maf_nonsmoking, maf_list$MSK2015[maf_list$MSK2015$Unique_Patient_Identifier %in% nonsmoking_samples])
maf_nonsmoking <- rbind(maf_nonsmoking, maf_list$OncoSG[maf_list$OncoSG$Unique_Patient_Identifier %in% nonsmoking_samples])
maf_nonsmoking <- rbind(maf_nonsmoking, maf_list$TCGA[maf_list$TCGA$Unique_Patient_Identifier %in% nonsmoking_samples])
maf_nonsmoking <- rbind(maf_nonsmoking, maf_list$TracerX[maf_list$TracerX$Unique_Patient_Identifier %in% nonsmoking_samples])

cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = maf_nonsmoking)

cesa_nonsmoking <- trinuc_mutation_rates(cesa_nonsmoking,
                                      signature_set = "COSMIC_v3.2",
                                      signatures_to_remove = signatures_to_remove
)
cesa_nonsmoking <- gene_mutation_rates(cesa_nonsmoking, covariates = "lung")

save_cesa(cesa_nonsmoking, 'nonsmoking_samples_cesa.rds')

smoking_mus <- cesa_smoking$gene_rates
nonsmoking_mus <- cesa_nonsmoking$gene_rates
smoking_mus[, ratio := rate / nonsmoking_mus[,rate]]
summary(smoking_mus$ratio)

plot <- ggplot(smoking_mus, aes(x = 1:nrow(smoking_mus), y = ratio)) + geom_violin()