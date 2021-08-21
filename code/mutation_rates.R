library(cancereffectsizeR)
library(ces.refset.hg19)
library(data.table)
library(dplyr)
library(ggplot2)

maf_file <- read.csv("/Users/Krishna1/Desktop/Research/lung-smoking.nosync/output/merged_final.txt")
LUAD_maf <- cancereffectsizeR::preload_maf(maf = maf_file, refset = ces.refset.hg19, chain_file = "/Users/Krishna1/Desktop/Research/hg38ToHg19.over.chain", keep_extra_columns = c('Source'))
LUAD_maf <- LUAD_maf[is.na(problem)]

LUAD_maf <- LUAD_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

# possible_dups <- check_sample_overlap(LUAD_maf)
# possible_dups_filtered <- possible_dups[variants_shared > 2,]
# 
# source_A <- maf_file[maf_file$Tumor_Sample_Barcode %in% possible_dups_filtered$sample_A, c("Tumor_Sample_Barcode","Source")]
# source_A <- source_A %>% distinct(Tumor_Sample_Barcode, .keep_all = T)
# 
# source_B <- maf_file[maf_file$Tumor_Sample_Barcode %in% possible_dups_filtered$sample_B, c("Tumor_Sample_Barcode","Source")]
# source_B <- source_B %>% distinct(Tumor_Sample_Barcode, .keep_all = T)
# 
# possible_dups_filtered$Source_A <- ''
# 
# for (i in 1:dim(possible_dups_filtered)[1]) {
#   for (j in 1:dim(source_A)[1]) {
#     if (possible_dups_filtered[i, 'sample_A'] == source_A[j, 'Tumor_Sample_Barcode']) {
#       possible_dups_filtered[i, 'Source_A'] <- source_A[j, 'Source']
#     }
#   }
# }
# 
# possible_dups_filtered$Source_B <- ''
# 
# for (i in 1:dim(possible_dups_filtered)[1]) {
#   for (j in 1:dim(source_B)[1]) {
#     if (possible_dups_filtered[i, 'sample_B'] == source_B[j, 'Tumor_Sample_Barcode']) {
#       possible_dups_filtered[i, 'Source_B'] <- source_B[j, 'Source']
#     }
#   }
# }

samples_to_remove <- c("luad_tsp_16929", "luad_tsp_16901", "luad_tsp_16875","luad_tsp_16915","LUAD-B01169","LUAD-D01382")
LUAD_maf <- LUAD_maf[!Unique_Patient_Identifier %in% samples_to_remove]

maf_list <- split(LUAD_maf, LUAD_maf$Source)

cesa_exome <- CESAnalysis()
cesa_exome <- load_maf(cesa_exome, maf = maf_list$Broad)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$MSK2015)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$OncoSG)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$TCGA)
cesa_exome <- load_maf(cesa_exome, maf = maf_list$TracerX)

signatures_to_remove <- suggest_cosmic_signatures_to_remove(cancer_type = "LUAD", treatment_naive = F)

cesa_exome <- trinuc_mutation_rates(cesa_exome,
                              signature_set = "COSMIC_v3.2",
                              signatures_to_remove = signatures_to_remove
)

bio_weights <- cesa_exome$mutational_signatures$biological_weights
bio_weights_unblended <- bio_weights[bio_weights$group_avg_blended == F]
bio_weights_blended <- bio_weights[bio_weights$group_avg_blended == T]

snv_counts = cesa_exome$maf[variant_type == 'snv', .N, by = "Unique_Patient_Identifier"]
#snv_counts[N <50,]# quantile(N, seq(0, 1, .1))]
  
sorted_bio_weights_unblended <- bio_weights_unblended[order(SBS4),]
#sorted_bio_weights_unblended[SBS4 > 0, SBS4]

good_samples <- snv_counts[N > 50, Unique_Patient_Identifier]
good_sample_weights <- sorted_bio_weights_unblended[Unique_Patient_Identifier %in% good_samples,]
smoking_samples <- good_sample_weights[SBS4 > 0, Unique_Patient_Identifier]
nonsmoking_samples <- good_sample_weights[SBS4 == 0, Unique_Patient_Identifier]


#SMOKING GROUP
cesa_smoking <- CESAnalysis()
# maf_smoking <- data.table()
# sapply(maf_list, function(x){
#   print(x[x$Unique_Patient_Identifier %in% smoking_samples])
#   })
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

save_cesa(cesa_smoking, 'smoking_samples_cesa.rds')
save_cesa(cesa_nonsmoking, 'nonsmoking_samples_cesa.rds')

smoking_mus <- cesa_smoking$gene_rates
nonsmoking_mus <- cesa_nonsmoking$gene_rates
smoking_mus[, ratio := rate / nonsmoking_mus[,rate]]
summary(smoking_mus$ratio)

plot <- ggplot(smoking_mus, aes(x = 1:nrow(smoking_mus), y = ratio)) + geom_violin()

#cesa_nonsmoking$gene_rates[gene %in% c('TP53','KRAS','STK11','KEAP1','EGFR','BRAF'), ]

# cesa <- CESAnalysis()
# cesa <- load_maf(cesa, maf = maf_list)
# 
# signatures_to_remove <- suggest_cosmic_signatures_to_remove(cancer_type = "LUAD", treatment_naive = F)
# 
# cesa <- trinuc_mutation_rates(cesa,
#                               signature_set = "COSMIC_v3.2",
#                               signatures_to_remove = signatures_to_remove
# )
# 
# bio_weights <- cesa$mutational_signatures$biological_weights
# bio_weights_unblended <- bio_weights[bio_weights$group_avg_blended == F]
# sorted_bio_weights_unblended <- bio_weights_unblended[order(SBS4),]
# sorted_bio_weights_unblended[SBS4 > 0, SBS4]
# 
# plot <- ggplot(sorted_bio_weights_unblended, aes(x = 1:nrow(sorted_bio_weights_unblended), y = SBS4)) + geom_point()
# 
# cesa <- gene_mutation_rates(cesa, covariates = "lung")
# 
# #genes_in_intersection <- c('BRAF', 'ERBB4', 'NF2', 'ETV6', 'BRCA1', 'TSHR', 'MPL', 'MSH2', 'PIK3R2', 'JUN', 'PDGFRA', 'MAP2K1', 'ABL1', 'CCND3', 'FGFR4', 'ATM', 'PTEN', 'SUFU', 'GNAS', 'NTRK1', 'KIT', 'BCL2', 'FLT1', 'STK11', 'CDKN1B', 'NRAS', 'DDR2', 'CCND1', 'JAK1', 'ERG', 'JAK2', 'MYC', 'AKT2', 'PIK3C2G', 'MSH6', 'CHEK1', 'NTRK3', 'ALK', 'NOTCH4', 'GSK3B', 'RAF1', 'EPHA5', 'TP53', 'NTRK2', 'FLT4', 'NOTCH3', 'PIK3C3', 'CDK8', 'PIK3R1', 'MEN1', 'EGFR', 'IGF1R', 'FGFR2', 'CTNNB1', 'IKBKE', 'CSF1R', 'BAP1', 'CDH1', 'CDKN2B', 'KRAS', 'AKT1', 'FGF3', 'HRAS', 'PTPN11', 'MAP2K2', 'BTK', 'MDM2', 'BCL6', 'CHEK2', 'NF1', 'FGF4', 'CDKN2C', 'CCNE1', 'KDR', 'RUNX1', 'CRKL', 'EPHB1', 'BARD1', 'BRCA2', 'MDM4', 'RB1', 'SRC', 'FBXW7', 'FGFR1', 'TSC2', 'FGFR3', 'ERBB3', 'CDKN2A', 'MAP3K1', 'PRKAR1A', 'REL', 'PIK3CA', 'MLH1', 'SMAD4', 'AXL', 'CDK6', 'MAP3K13', 'ATR', 'TSC1', 'VHL', 'RET', 'AURKB', 'APC', 'MAP2K4', 'CDK4', 'CCND2', 'MYCN', 'PIK3CG', 'SMAD2', 'MET', 'WT1', 'NOTCH1', 'JAK3', 'NOTCH2', 'TGFBR2', 'ARAF', 'EPHA3', 'CBL', 'IRS2', 'ERBB2', 'PDGFRB', 'SMO', 'GATA1', 'AKT3', 'FLT3', 'SYK', 'RICTOR')
# 
# #fwrite(cesa$gene_rates[gene %in% genes_in_intersection], '../data/mutation_rates.txt')
# 
# cesa <- ces_variant(cesa = cesa)
# save_cesa(cesa, '../data/lung_smoking.rds')

