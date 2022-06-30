library(ggplot2)
library(data.table)
library(stringr)


if(!file.exists(paste0(location_data,'gene_prevalences.rds'))){
  pan_data_cesa = load_cesa('../../data/pan_data_cesa_for_cancer_epistasis.rds')
  #' The "plus" refers to the inclusion of clinically annotated panel data
  smoking_plus_cesa = load_cesa('../../data/smoking_samples_w_panel_cesa.rds')
  nonsmoking_plus_cesa = load_cesa('../../data/nonsmoking_samples_w_panel_cesa.rds') 
  
  
  gene_prevalences = list(pan_data_gene_prevalence = unique(pan_data_cesa$variants[,.(gene, maf_prevalence)][,pd_prevalence := sum(maf_prevalence),by=gene][,.(gene, pd_prevalence)]),
                          smoking_plus_gene_prevalence = unique(smoking_plus_cesa$variants[,.(gene, maf_prevalence)][,`s+_prevalence` := sum(maf_prevalence),by=gene][,.(gene, `s+_prevalence`)]),
                          nonsmoking_plus_gene_prevalence = unique(nonsmoking_plus_cesa$variants[,.(gene, maf_prevalence)][,`ns+_prevalence` := sum(maf_prevalence),by=gene][,.(gene, `ns+_prevalence`)]))
  
  saveRDS(gene_prevalences, paste0(location_data,'gene_prevalences.rds'))
} else {
  gene_prevalences = readRDS(paste0(location_data,'gene_prevalences.rds'))
}



genes_by_gamma = fread('../../output/genes_by_gamma.csv')
colnames(genes_by_gamma) = c('gene',colnames(genes_by_gamma)[2:ncol(genes_by_gamma)])
genes_by_gamma = merge.data.table(genes_by_gamma, gene_prevalences$pan_data_gene_prevalence,
                                  by = 'gene',all.x = T, all.y = F)
genes_by_gamma = merge.data.table(genes_by_gamma, gene_prevalences$smoking_plus_gene_prevalence,
                                  by = 'gene',all.x = T, all.y = F)
genes_by_gamma = merge.data.table(genes_by_gamma, gene_prevalences$nonsmoking_plus_gene_prevalence,
                                  by = 'gene',all.x = T, all.y = F)
genes_by_gamma[, ':=' (pd_from_WT_prevalence = pd_prevalence,
                       `s+_from_WT_prevalence` = `s+_prevalence`,
                       `ns+_from_WT_prevalence` = `ns+_prevalence`)]

pan_data_gammas = genes_by_gamma[,.(gene, pd_from_110_gamma, pd_from_110_low, pd_from_110_high, pd_from_WT_gamma, pd_from_WT_low, pd_from_WT_high, pd_prevalence, pd_from_WT_prevalence)]

#' PAN-DATA
significant_gammas = pan_data_gammas[pd_from_110_low > pd_from_WT_high | pd_from_WT_low > pd_from_110_high]
temp1 = significant_gammas[pd_from_110_gamma > pd_from_WT_gamma][order(-pd_from_110_gamma)]
temp2 = significant_gammas[pd_from_110_gamma < pd_from_WT_gamma][order(-pd_from_WT_gamma)]
significant_gammas = rbind(temp1, temp2)

pivoted_table = data.table::melt(significant_gammas, 
                                 id.vars = 'gene', measure.vars = patterns("low$","gamma$",'high$','prevalence$'),
                                 variable.name = 'dataset', variable.n, value.name = c('ci_low_95','selection_intensity','ci_high_95','prevalence'))
levels = c('from_110', 'from_WT')
names(levels) = 1:2
pivoted_table$dataset = as.vector(levels[pivoted_table$dataset])
pivoted_table$dataset = factor(pivoted_table$dataset, levels = as.vector(levels[c(2,1)]))
pivoted_table$gene = factor(pivoted_table$gene, levels = significant_gammas$gene)

breaks <- unique(as.numeric(round(quantile(pivoted_table$prevalence))))
ggplot(data = pivoted_table, aes(x = gene, y = selection_intensity, fill = dataset)) + 
  coord_flip() +
  scale_x_discrete(limits = rev) + 
  geom_point(aes(color = prevalence, shape = dataset), size = 2, position = position_dodge(0.5)) + 
  scale_shape_manual(values = c(16, 18)) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = "dark grey", position = position_dodge(0.5)) + 
  geom_vline(xintercept = 1:nrow(significant_gammas) + 0.5, size = 0.2, color = 'grey') +
  geom_vline(xintercept = 2.5, color = 'red') +
  geom_vline(xintercept = seq(1, nrow(significant_gammas), 5) -0.5, size = 0.3, color = 'dark grey') + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = 'light grey', size = 0.2), panel.background = element_rect(fill = 'white')) + 
  ggtitle('Epistatic Effects of TP53 and KRAS in Lung Adenocarcinoma') + 
  xlab('Gene') +
  ylab('Gamma (Selection Intensity)') + 
  scale_color_viridis_c(
    name = "variant prevalence", guide = "colorbar", trans = "log10",
    option = "plasma", breaks = breaks
  ) + 
  guides(color = guide_colourbar(ticks = T))


#' SMOKING/NONSMOKING
smoking_gammas = genes_by_gamma[,-c('pd_from_110_gamma', 'pd_from_110_low', 'pd_from_110_high', 'pd_from_WT_gamma', 'pd_from_WT_low', 'pd_from_WT_high', 'pd_prevalence', 'pd_from_WT_prevalence')]
smoking_gammas = smoking_gammas[gene != 'BRD3']

ci_lows = as.data.frame(smoking_gammas)[,str_which(colnames(smoking_gammas), 'low$')]
ci_highs = as.data.frame(smoking_gammas)[,str_which(colnames(smoking_gammas), 'high$')]

#' The line below tests significance by checking if the highest low of a CI overlaps
#' with the smallest high of a CI. However it lets the majority of variants through
#' so a more conservative approach is needed
#significant_gammas = smoking_gammas[apply(ci_lows, 1, max, na.rm = T) > apply(ci_highs, 1, min, na.rm = T)]

#' More conservative condition for significant gammas
#' If the CI of the max value for any gene overlaps with the CIs of any other values 
#' for that gene, remove that gene
#' Implemented by checking if max of ci_lows is greater than 2nd highest ci_high
significant_gammas = smoking_gammas[apply(ci_lows, 1, max, na.rm = T) > apply(ci_highs, 1, function(x){sort(x, decreasing = T)[2]})]

#' Excluding genes that don't have values for all 4 gammas
significant_gammas = significant_gammas[!(is.na(`ns+_from_WT_gamma`) | is.na(`s+_from_WT_gamma`))]

gamma_values = as.data.frame(significant_gammas)[,str_which(colnames(significant_gammas), 'gamma$')]
gamma_values = cbind(gene = significant_gammas$gene, gamma_values)
#identifies column index which has the max value for that row
gamma_values$max_val_col = apply(gamma_values[,-1], 1, which.max) + 1
gamma_values = as.data.table(gamma_values)

sorted_genes = list('smoking_from_110' = gamma_values[max_val_col == 2][order(-gamma_values[max_val_col == 2,2]), gene],
                     'nonsmoking_from_110' = gamma_values[max_val_col == 3][order(-gamma_values[max_val_col == 3,3]), gene],
                     'smoking_from_WT' = gamma_values[max_val_col == 4][order(-gamma_values[max_val_col == 4,4]), gene],
                     'nonsmoking_from_WT' = gamma_values[max_val_col == 5][order(-gamma_values[max_val_col == 5,5]), gene])
  
pivoted_table = data.table::melt(significant_gammas, 
                                 id.vars = 'gene', measure.vars = patterns("low$","gamma$",'high$','prevalence$'),
                                 variable.name = 'dataset', variable.n, value.name = c('ci_low_95','selection_intensity','ci_high_95','prevalence'))
levels = c('smoking_plus_from_110','nonsmoking_plus_from_110', 'smoking_plus_from_WT', 'nonsmoking_plus_from_WT')
names(levels) = 1:4
pivoted_table$dataset = as.vector(levels[pivoted_table$dataset])
pivoted_table$dataset = factor(pivoted_table$dataset, levels = as.vector(levels[c(3,1,4,2)]))
pivoted_table$gene = factor(pivoted_table$gene, levels = as.vector(unlist(sorted_genes)))

breaks <- unique(as.numeric(round(quantile(pivoted_table$prevalence, na.rm = T))))
ggplot(data = pivoted_table, aes(x = gene, y = selection_intensity, fill = dataset)) + 
  coord_flip() +
  scale_x_discrete(limits = rev) + 
  geom_point(aes(color = prevalence, shape = dataset), size = 3, position = position_dodge(1)) + 
  scale_shape_manual(values = 15:18) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = "darkgrey", position = position_dodge(1)) + 
  geom_vline(xintercept = 1:nrow(significant_gammas) + 0.5, size = 0.2, color = 'grey') +
  geom_vline(xintercept = nrow(significant_gammas) - length(sorted_genes[[1]]) + 0.5, color = 'red') +
  geom_vline(xintercept = nrow(significant_gammas) - (length(sorted_genes[[1]]) + length(sorted_genes[[2]])) + 0.5, color = 'red') +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = 'light grey', size = 0.2), panel.background = element_rect(fill = 'white')) + 
  ggtitle('Epistatic Effects of TP53 and KRAS in Lung Adenocarcinoma Between Smokers and Nonsmokers') + 
  xlab('Gene') +
  ylab('Gamma (Selection Intensity)') + 
  scale_color_viridis_c(
    name = "variant prevalence", guide = "colorbar", trans = "log10",
    option = "plasma", breaks = breaks
  ) + 
  guides(color = guide_colourbar(ticks = T))
