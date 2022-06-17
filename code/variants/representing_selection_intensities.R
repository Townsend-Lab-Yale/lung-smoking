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
  gene_prevalences = read_RDS(paste0(location_data,'gene_prevalences.rds'))
}



genes_by_gamma = fread('../../output/genes_by_gamma.csv')
colnames(genes_by_gamma) = c('gene',colnames(genes_by_gamma)[2:ncol(genes_by_gamma)])
genes_by_gamma = merge.data.table(genes_by_gamma, gene_prevalences$pan_data_gene_prevalence,
                                  by = 'gene',all.x = T, all.y = F)
genes_by_gamma = merge.data.table(genes_by_gamma, gene_prevalences$smoking_plus_gene_prevalence,
                                  by = 'gene',all.x = T, all.y = F)
genes_by_gamma = merge.data.table(genes_by_gamma, gene_prevalences$nonsmoking_plus_gene_prevalence,
                                  by = 'gene',all.x = T, all.y = F)

pivoted_table = data.table::melt(genes_by_gamma, 
                                 id.vars = 'gene', measure.vars = patterns("low$","gamma$",'high$','prevalence$'),
                                 variable.name = 'dataset', variable.n, value.name = c('ci_low_95','selection_intensity','ci_high_95','prevalence'))
pivoted_table[,dataset := fifelse(dataset == 1, 'pan_data', fifelse(dataset == 2, 'smoking_plus', 'nonsmoking_plus'))]

pivoted_table$dataset = factor(pivoted_table$dataset, levels = unique(pivoted_table$dataset))

breaks <- unique(as.numeric(round(quantile(pivoted_table$pd_prevalence))))
ggplot(data = pivoted_table, aes(x = reorder(gene, -selection_intensity), y = selection_intensity, fill = dataset)) + 
  geom_point(aes(color = prevalence, shape = dataset), size = 3, position = position_dodge(1)) + 
  scale_shape_manual(values = 16:18) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = "darkgrey", position = position_dodge(1)) + 
  geom_vline(xintercept = 1:nrow(genes_by_gamma) + 0.5, color = 'black') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = 'light grey'), panel.background = element_rect(fill = 'white')) + 
  ggtitle('Selection Intensity of Genes in Smoker and Nonsmoker Lung Adenocarcinoma') + 
  xlab('Gene') +
  ylab('Gamma (Selection Intensity)') + 
  scale_color_viridis_c(
    name = "variant prevalence", guide = "colorbar", trans = "log10",
    option = "plasma", breaks = breaks
  ) + 
  guides(color = guide_colourbar(ticks = FALSE))
