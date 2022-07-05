library(cancereffectsizeR)
library(data.table)
library(ggplot2)
library(ces.refset.hg19)

location_data = '~/Desktop/Research/lung-smoking/data/'

smoking_samples = fread(paste0(location_data, 'smoking_sample_ids.txt'), header = F)$V1
nonsmoking_samples = fread(paste0(location_data, 'nonsmoking_sample_ids.txt'), header = F)$V1
panel_smoking_samples = fread(paste0(location_data, 'panel_smoking_sample_ids.txt'), header = F)$V1
panel_nonsmoking_samples = fread(paste0(location_data, 'panel_nonsmoking_sample_ids.txt'), header = F)$V1

cesa = load_cesa(paste0(location_data, 'pan_data_cesa_for_cancer_epistasis.rds'))

cesa = clear_trinuc_rates_and_signatures(cesa)
cesa = clear_gene_rates(cesa)

signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = 'LUAD', treatment_naive = F)

#' Rates for smokers
cesa <- trinuc_mutation_rates(cesa, samples = c(smoking_samples, panel_smoking_samples),
                              signature_set = "COSMIC_v3.2",
                              signature_exclusions = signature_exclusions)
cesa <- gene_mutation_rates(cesa, samples = c(smoking_samples, panel_smoking_samples), covariates = "lung")

#' Rates for nonsmokers
cesa <- trinuc_mutation_rates(cesa, samples = c(nonsmoking_samples, panel_nonsmoking_samples),
                              signature_set = "COSMIC_v3.2",
                              signature_exclusions = signature_exclusions)
cesa <- gene_mutation_rates(cesa, samples = c(nonsmoking_samples, panel_nonsmoking_samples), covariates = "lung")




compute_gene_selection = function(cesa, genes){
  #' Setting maf_prevalence > 0, thought > 1 may be possible and/or preferable as well
  print('Compounding variants...')
  for_compound = cesa$variants[gene %in% genes & maf_prevalence >0]
  #' Since the only samples we are considering are those with smoker determinable
  #' status (only exome/genome and MSK TGS samples), we only add variants that have
  #' coverage by these panels. We can ignore the other panel's not overlapping
  #' because we won't be using their panels anyway.
  #' 
  #' For the genes that are covered by the MSK-IMPACT panels (20 genes), 
  #' only add the variants that are covered by the panels. 
  #' For the genes that we want that
  for_compound_1 = for_compound[unlist(lapply(for_compound$covered_in, function(x){any(c('msk341_regions','msk410_regions','msk468_regions') %in% x)}))]
  for_compound_2 = for_compound[gene %in% setdiff(unique(for_compound$gene), unique(for_compound_1$gene))]
  for_compound_2 = for_compound_2[unlist(lapply(for_compound_2$covered_in, function(x){'exome+' %in% x}))]
  for_compound = rbind(for_compound_1, for_compound_2)
  
  comp = define_compound_variants(cesa, variant_table = for_compound, by = 'gene', merge_distance = Inf)
  print('Done compounding variants\nComputing cancer effect size...')
  
  cesa = ces_variant(cesa = cesa, variants = comp,
                     samples = c(smoking_samples, panel_smoking_samples), run_name = 'smoking')
  cesa = ces_variant(cesa = cesa, variants = comp,
                     samples = c(nonsmoking_samples, panel_nonsmoking_samples), run_name = 'nonsmoking')
  cesa = ces_variant(cesa = cesa, variants = comp,
                     samples = c(smoking_samples, nonsmoking_samples, panel_smoking_samples, panel_nonsmoking_samples), run_name = 'combined')
  
  return(cesa)
}

#' This line outputs the proportion of variants retained for each gene
#table(for_compound$gene) / 
#  table(cesa$variants[gene %in% genes & maf_prevalence > 0, gene])

genes = as.vector(unlist(sorted_genes)) #sorted_genes from representing_selection_intensities.R
cesa = compute_gene_selection(cesa, genes)


two_group_lik_sum = cesa$selection$smoking$loglikelihood + cesa$selection$nonsmoking$loglikelihood
combined_lik = cesa$selection$combined$loglikelihood
chisquared = -2 * (combined_lik - two_group_lik_sum)
p = pchisq(chisquared, df = 1, lower.tail = F)

signif_genes = cesa$selection$combined$variant_name[p<.05]

si_tables = lapply(cesa$selection[1:2], function(x){x[variant_name %in% signif_genes, .(variant_name, selection_intensity, ci_low_95, ci_high_95, included_with_variant)]})
si_tables$smoking[, status := 'smoking']
si_tables$nonsmoking[, status := 'nonsmoking']

si_table = rbindlist(si_tables)
colnames(si_table) = c('gene', colnames(si_table)[2:4], 'prevalence', 'status')
si_table$gene = gsub('\\.1','',si_table$gene)
si_table$gene = factor(si_table$gene, levels = c('PBRM1','ATF7IP','ARID1A','ATM','STK11','TLR4','EGFR'))

breaks <- unique(as.numeric(round(quantile(si_table$prevalence))))
ggplot(data = si_table, aes(x = gene, y = selection_intensity, fill = status)) + 
  coord_flip() +
  scale_x_discrete(limits = rev) + 
  geom_point(aes(color = prevalence, shape = status), size = 2, position = position_dodge(0.5)) + 
  scale_shape_manual(values = c(16, 18)) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = "dark grey", position = position_dodge(0.5)) + 
  geom_vline(xintercept = 1:(nrow(si_table)/2) + 0.5, size = 0.2, color = 'grey') +
  geom_vline(xintercept = 2.5, color = 'red') +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = 'light grey', size = 0.2), panel.background = element_rect(fill = 'white')) + 
  ggtitle('Cancer Effect Sizes of TP53 and KRAS in Smoker and Nonsmoker Lung Adenocarcinoma') + 
  xlab('Gene') +
  ylab('Effect Size') + 
  scale_color_viridis_c(
    name = "variant prevalence", guide = "colorbar", trans = "log10",
    option = "plasma", breaks = breaks
  ) + 
  guides(color = guide_colourbar(ticks = T))




#' APC
apc_comp = define_compound_variants(cesa, variant_table = cesa$variants[gene == 'APC' & samples_covering > 6000 & maf_prevalence > 0], merge_distance = Inf)

cesa = ces_variant(cesa = cesa, variants = apc_comp,
                   samples = c(smoking_samples, panel_smoking_samples), run_name = 'APC_smoking')
cesa = ces_variant(cesa = cesa, variants = apc_comp,
                   samples = c(nonsmoking_samples, panel_nonsmoking_samples), run_name = 'APC_nonsmoking')
cesa = ces_variant(cesa = cesa, variants = apc_comp,
                   samples = c(smoking_samples, nonsmoking_samples, panel_smoking_samples, panel_nonsmoking_samples), run_name = 'APC_combined')

two_group_lik_sum = cesa$selection$APC_smoking$loglikelihood + cesa$selection$APC_nonsmoking$loglikelihood
combined_lik = cesa$selection$APC_combined$loglikelihood
chisquared = -2 * (combined_lik - two_group_lik_sum)
p = pchisq(chisquared, df = 1, lower.tail = F)
#p = 0.3372228
si_table = as.data.table(rbind(cesa$selection$APC_smoking, cesa$selection$APC_nonsmoking, cesa$selection$APC_combined)[,.(selection_intensity,ci_low_95,ci_high_95,included_with_variant)])
si_table = cbind(c('smoking','nonsmoking','all'), si_table)
si_table 



#' SETD2
setd2_comp = define_compound_variants(cesa, variant_table = cesa$variants[gene == 'SETD2' & samples_covering > 6000 & maf_prevalence > 0], merge_distance = Inf)

cesa = ces_variant(cesa = cesa, variants = setd2_comp,
                   samples = c(smoking_samples, panel_smoking_samples), run_name = 'SETD2_smoking')
cesa = ces_variant(cesa = cesa, variants = setd2_comp,
                   samples = c(nonsmoking_samples, panel_nonsmoking_samples), run_name = 'SETD2_nonsmoking')
cesa = ces_variant(cesa = cesa, variants = setd2_comp,
                   samples = c(smoking_samples, nonsmoking_samples, panel_smoking_samples, panel_nonsmoking_samples), run_name = 'SETD2_combined')

two_group_lik_sum = cesa$selection$SETD2_smoking$loglikelihood + cesa$selection$SETD2_nonsmoking$loglikelihood
combined_lik = cesa$selection$SETD2_combined$loglikelihood
chisquared = -2 * (combined_lik - two_group_lik_sum)
p = pchisq(chisquared, df = 1, lower.tail = F)
#p = 0.1391215
si_table = as.data.table(rbind(cesa$selection$SETD2_smoking, cesa$selection$SETD2_nonsmoking, cesa$selection$SETD2_combined)[,.(selection_intensity,ci_low_95,ci_high_95,included_with_variant)])
si_table = cbind(c('smoking','nonsmoking','all'), si_table)
si_table 

