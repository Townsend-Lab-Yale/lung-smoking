library(cancereffectsizeR)
library(ces.refset.hg19)
library(data.table)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(stringr)
# All files not available in repo but still used here:
# Gencode basic annotation,
# Genie genomic information

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

Genie_maf <- cancereffectsizeR::preload_maf(maf = maf_list$Genie, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
Genie_maf <- Genie_maf[is.na(problem)]
Genie_maf <- Genie_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

MSK2015_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2015, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
MSK2015_maf <- MSK2015_maf[is.na(problem)]
MSK2015_maf <- MSK2015_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

MSK2017_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2017, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
MSK2017_maf <- MSK2017_maf[is.na(problem)]
MSK2017_maf <- MSK2017_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

MSK2018_maf <- cancereffectsizeR::preload_maf(maf = maf_list$MSK2018, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
MSK2018_maf <- MSK2018_maf[is.na(problem)]
MSK2018_maf <- MSK2018_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

OncoSG_maf <- cancereffectsizeR::preload_maf(maf = maf_list$OncoSG, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
OncoSG_maf <- OncoSG_maf[is.na(problem)]
OncoSG_maf <- OncoSG_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TCGA_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TCGA, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
TCGA_maf <- TCGA_maf[is.na(problem)]
TCGA_maf <- TCGA_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TracerX_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TracerX, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
TracerX_maf <- TracerX_maf[is.na(problem)]
TracerX_maf <- TracerX_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

TSP_maf <- cancereffectsizeR::preload_maf(maf = maf_list$TSP, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")
TSP_maf <- TSP_maf[is.na(problem)]
TSP_maf <- TSP_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]
TSP_maf <- TSP_maf[!Unique_Patient_Identifier %in% c("luad_tsp_16929", "luad_tsp_16901", "luad_tsp_16875","luad_tsp_16915")]

#maf_list <- sapply(maf_list, function(x){preload_maf(maf = x, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain")})


#SIGNATURES NOT RELEVANT TO LUAD
signatures_to_remove <- suggest_cosmic_signatures_to_remove(cancer_type = "LUAD", treatment_naive = F)


#READING IN INFORMATION NECESSARY FOR LOADING MAF FILES (EXOME/GENOME AND COVERAGE INTERVALS FOR TGS)
broad_exome_or_genome <- fread('../data/luad_broad/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Platform')]
broad_exome_or_genome$Platform <- sapply(broad_exome_or_genome$Platform, function(x){if(str_detect(x, 'WGS')){return('WGS')} else return('WES')})
msk2017_panels_used <- fread('../data/lung_msk_2017/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Gene Panel')]
msk2018_panels_used <- fread('../data/nsclc_pd1_msk_2018/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Gene Panel')]
genie_panels_used <- fread('../data/genie_9/data_clinical_sample.txt')[-(1:4),c('Sample Identifier', 'Sequence Assay ID')]


#READING IN GENES INCLUDED IN EACH PANEL AND CREATING GRANGES OBJECT TO PASS INTO COVERED_REGIONS PARAMETER OF LOAD_MAF
#once the granges are exported once, these functions don't need to be run anymore
gene_granges <- rtracklayer::import('/Users/Krishna1/Downloads/gencode.v38lift37.basic.annotation.gtf')

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


#GENES INCLUDED IN EVERY PANEL (NOT TRUE FOR ALL THE PANELS IN GENIE THOUGH, SEE MESSAGES WITH JORGE)
genes_in_intersection <- c('BRAF', 'ERBB4', 'NF2', 'ETV6', 'BRCA1', 'TSHR', 'MPL', 'MSH2', 'PIK3R2', 'JUN', 'PDGFRA', 'MAP2K1', 'ABL1', 'CCND3', 'FGFR4', 'ATM', 'PTEN', 'SUFU', 'GNAS', 'NTRK1', 'KIT', 'BCL2', 'FLT1', 'STK11', 'CDKN1B', 'NRAS', 'DDR2', 'CCND1', 'JAK1', 'ERG', 'JAK2', 'MYC', 'AKT2', 'PIK3C2G', 'MSH6', 'CHEK1', 'NTRK3', 'ALK', 'NOTCH4', 'GSK3B', 'RAF1', 'EPHA5', 'TP53', 'NTRK2', 'FLT4', 'NOTCH3', 'PIK3C3', 'CDK8', 'PIK3R1', 'MEN1', 'EGFR', 'IGF1R', 'FGFR2', 'CTNNB1', 'IKBKE', 'CSF1R', 'BAP1', 'CDH1', 'CDKN2B', 'KRAS', 'AKT1', 'FGF3', 'HRAS', 'PTPN11', 'MAP2K2', 'BTK', 'MDM2', 'BCL6', 'CHEK2', 'NF1', 'FGF4', 'CDKN2C', 'CCNE1', 'KDR', 'RUNX1', 'CRKL', 'EPHB1', 'BARD1', 'BRCA2', 'MDM4', 'RB1', 'SRC', 'FBXW7', 'FGFR1', 'TSC2', 'FGFR3', 'ERBB3', 'CDKN2A', 'MAP3K1', 'PRKAR1A', 'REL', 'PIK3CA', 'MLH1', 'SMAD4', 'AXL', 'CDK6', 'MAP3K13', 'ATR', 'TSC1', 'VHL', 'RET', 'AURKB', 'APC', 'MAP2K4', 'CDK4', 'CCND2', 'MYCN', 'PIK3CG', 'SMAD2', 'MET', 'WT1', 'NOTCH1', 'JAK3', 'NOTCH2', 'TGFBR2', 'ARAF', 'EPHA3', 'CBL', 'IRS2', 'ERBB2', 'PDGFRB', 'SMO', 'GATA1', 'AKT3', 'FLT3', 'SYK', 'RICTOR')


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
                                    signature_set = "COSMIC_v3.2",
                                    signatures_to_remove = signatures_to_remove
)
cesa_total <- gene_mutation_rates(cesa_total, covariates = "lung")

save_cesa(cesa_total, '../data/pan-dataset_samples_cesa.rds')
fwrite(cesa_total$gene_rates[gene %in% genes_in_intersection], '../data/pan-data_mutation_rates.txt')


#WES/WGS-ONLY CESA ANALYSIS BECAUSE ONLY THESE CAN HAVE SIGNATURE EXTRACTIONS BE PERFORMED ON THEM
cesa_exome <- CESAnalysis()
cesa_exome <- load_maf(cesa_exome, maf = Broad_maf$WGS, coverage = 'genome')
cesa_exome <- load_maf(cesa_exome, maf = Broad_maf$WES)
cesa_exome <- load_maf(cesa_exome, maf = MSK2015_maf)
cesa_exome <- load_maf(cesa_exome, maf = OncoSG_maf)
cesa_exome <- load_maf(cesa_exome, maf = TCGA_maf)
cesa_exome <- load_maf(cesa_exome, maf = TracerX_maf)

cesa_exome <- trinuc_mutation_rates(cesa_exome,
                              signature_set = "COSMIC_v3.2",
                              signatures_to_remove = signatures_to_remove
)


#CALCULATING MUTATION RATES FOR COMPARISON PURPOSES WITH THE PAN-DATASET MUTATION RATES
cesa_exome <- gene_mutation_rates(cesa_exome, covariates = "lung")
save_cesa(cesa_exome, '../data/exome_samples_cesa.rds')
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
cesa_smoking <- load_maf(cesa_smoking, maf = Broad_maf$WES[Broad_maf$WES$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = MSK2015_maf[MSK2015_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = OncoSG_maf[OncoSG_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = TCGA_maf[TCGA_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = TracerX_maf[TracerX_maf$Unique_Patient_Identifier %in% smoking_samples])
cesa_smoking <- load_maf(cesa_smoking, maf = Broad_maf$WGS[Broad_maf$WGS$Unique_Patient_Identifier %in% smoking_samples], coverage = 'genome')

cesa_smoking <- trinuc_mutation_rates(cesa_smoking,
                                    signature_set = "COSMIC_v3.2",
                                    signatures_to_remove = signatures_to_remove
)

cesa_smoking <- gene_mutation_rates(cesa_smoking, covariates = "lung")

save_cesa(cesa_smoking, '../data/smoking_samples_cesa.rds')
fwrite(cesa_smoking$gene_rates[gene %in% genes_in_intersection], '../data/smoking_mutation_rates.txt')


#NONSMOKING GROUP
cesa_nonsmoking <- CESAnalysis()
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = Broad_maf$WES[Broad_maf$WES$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = MSK2015_maf[MSK2015_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = OncoSG_maf[OncoSG_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = TCGA_maf[TCGA_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = TracerX_maf[TracerX_maf$Unique_Patient_Identifier %in% nonsmoking_samples])
cesa_nonsmoking <- load_maf(cesa_nonsmoking, maf = Broad_maf$WGS[Broad_maf$WGS$Unique_Patient_Identifier %in% nonsmoking_samples], coverage = 'genome')

cesa_nonsmoking <- trinuc_mutation_rates(cesa_nonsmoking,
                                      signature_set = "COSMIC_v3.2",
                                      signatures_to_remove = signatures_to_remove
)
cesa_nonsmoking <- gene_mutation_rates(cesa_nonsmoking, covariates = "lung")
save_cesa(cesa_nonsmoking, '../data/nonsmoking_samples_cesa.rds')
fwrite(cesa_nonsmoking$gene_rates[gene %in% genes_in_intersection], '../data/nonsmoking_mutation_rates.txt')

# smoking_mus <- cesa_smoking$gene_rates
# nonsmoking_mus <- cesa_nonsmoking$gene_rates
# smoking_mus[, ratio := rate / nonsmoking_mus[,rate]]
# summary(smoking_mus$ratio)
# 
# plot <- ggplot(smoking_mus, aes(x = 1:nrow(smoking_mus), y = ratio)) + geom_violin() + labs(x = '', y = 'Mutation Rate Ratio\n(Smoking signature mu / Nonsmoking signature mu)', title = 'Comparison of Gene Mutation Rates Among Smokers\nand Nonsmokers in Lung Adenocarcinoma') + theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank(), axis.title.y = element_text(size = 8),plot.title = element_text(size = 10, hjust = 0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
# 
# subsetted_smoking_mus <- cesa_smoking$gene_rates[gene %in% genes_in_intersection]
# subsetted_nonsmoking_mus <- cesa_nonsmoking$gene_rates[gene %in% genes_in_intersection]
# subsetted_smoking_mus[, ratio := rate / subsetted_nonsmoking_mus[,rate]]
# summary(subsetted_smoking_mus$ratio)
# 
# subsetted_ratios_plot <- ggplot(subsetted_smoking_mus, aes(x = 1:nrow(subsetted_smoking_mus), y = ratio)) + geom_violin() + labs(x = '', y = 'Mutation Rate Ratio\n(Smoking signature mu / Nonsmoking signature mu)', title = 'Comparison of Gene Mutation Rates Among Smokers\nand Nonsmokers in Lung Adenocarcinoma') + theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank(), axis.title.y = element_text(size = 8),plot.title = element_text(size = 10, hjust = 0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
# 
# pan_v_exome_ratios <- cesa_total$gene_rates$rate / cesa_exome$gene_rates$rate
# summary(pan_v_exome_ratios)
# 
# genie_panel_frequency <-as.data.table(table(genie_panels_used$`Sequence Assay ID`))
# names(genie_panel_frequency) <- c('SEQ_ASSAY_ID','Num_Patients')
# #genie_panel_numgenes <- as.data.table(table(genie_panel_genes[!duplicated(genie_panel_genes, by = 'Hugo_Symbol')]$SEQ_ASSAY_ID))
# #names(genie_panel_numgenes) <- c('SEQ_ASSAY_ID','Num_Genes')
# #genie_panel_info <- merge(genie_panel_frequency, genie_panel_numgenes, by = 'SEQ_ASSAY_ID')
# #genie_panel_info[order(-Num_Genes, -Num_Patients)]
# 
# rate_comps <- merge(cesa_smoking$gene_rates, cesa_nonsmoking$gene_rates, by = 'gene')
# names(rate_comps) <- c('gene', 'smoking_rate', 'nonsmoking_rate')
# rate_comps[, ratio := smoking_rate / nonsmoking_rate]
# rate_comps[gene == 'ATM']
