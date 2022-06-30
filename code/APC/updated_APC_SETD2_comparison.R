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
  for_compound = rbindlist(lapply(genes, function(gene_x){cesa$variants[gene == gene_x & samples_covering >= cesa$variants[gene == gene_x, quantile(samples_covering)]['50%'] & maf_prevalence > 0]}))

  #' To preserve power, I am making the minimal samples_covering the first quartile because 
  #' ces_variant will only use samples with coverage at all sites, so if I include
  #' sites with less samples_covering, I will only be able to use that smaller
  #' set of samples.
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
table(for_compound$gene) / 
  table(cesa$variants[gene %in% genes & maf_prevalence > 0, gene])

#' Detects which genes have variants that are covered by different sets of panels
problem_genes = sapply(unique(for_compound$gene), function(x){uniqueN(lengths(for_compound[gene ==x, covered_in])) != 1})
problem_genes = names(problem_genes)[problem_genes == T]

genes = as.vector(unlist(sorted_genes)) #sorted_genes from representing_selection_intensities.R
cesa = compute_gene_selection(cesa, genes)






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



#' Any Gene


#' Commented lines still need to be vectorized
#two_group_lik_sum = cesa$selection$APC_smoking$loglikelihood + cesa$selection$APC_nonsmoking$loglikelihood
#combined_lik = cesa$selection$APC_combined$loglikelihood
chisquared = -2 * (combined_lik - two_group_lik_sum)
p = pchisq(chisquared, df = 1, lower.tail = F)

#si_table = as.data.table(rbind(cesa$selection$APC_smoking, cesa$selection$APC_nonsmoking, cesa$selection$APC_combined)[,.(selection_intensity,ci_low_95,ci_high_95,included_with_variant)])
si_table = cbind(c('smoking','nonsmoking','all'), si_table)
si_table 

