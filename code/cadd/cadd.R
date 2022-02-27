library(data.table)
library(cancereffectsizeR)
library(stringr)

#' Integrates pre-computed cadd scores into a cesa variant set
#' Returns a table with variant information and raw and scaled cadd scores

#' Note for further consideration: 
#'   since cadd will produce scores for variants that are not even represented in the dataset, 
#'   it can give us an indication of important variants we may still want to consider 
#'   even if they aren't represented in the dataset


integrate_cadd_scores = function(cesa, cadd_table){
  variants = c(unlist(cesa$variants[variant_type == 'aac', constituent_snvs]), unlist(cesa$variants[variant_type == 'snv', variant_id]))
  temp1 = str_split_fixed(variants, pattern = ':', n=2)
  temp2 = str_split_fixed(temp1[,2], pattern= '_', n=2)
  temp3 = str_split_fixed(temp2[,2], pattern= '>', n=2)
  
  variants = as.data.table(cbind(temp1[,1],temp2[,1],temp3))
  colnames(variants) = c('chr','pos','ref','tum')
  variants$pos = strtoi(variants$pos)
  variants = merge(variants, cadd_table, by.x = c('chr', 'pos', 'ref', 'tum'), by.y = c('#Chrom','Pos','Ref','Alt'), all.x=T)
  variants = variants[, variant := paste0(chr,':',pos,'_',ref,'>',tum)]
  
  cadd_variants = na.omit(variants)
  cadd_variants[, variant := paste0(chr,':',pos,'_',ref,'>',tum)]
  
  snv_variants = cesa$variants[variant_type == 'snv', .(gene, variant_name, variant_id, variant_type, aachange)]
  snv_variants = merge.data.table(snv_variants, cadd_variants, by.x = 'variant_id', by.y = 'variant')
  
  aac_variants = cesa$variants[variant_type == 'aac', .(constituent_snvs = unlist(constituent_snvs), gene, variant_id, variant_type, aachange), by = "variant_name"]
  aac_variants[,variant_id := constituent_snvs][,constituent_snvs := NULL]
  aac_variants = merge.data.table(aac_variants, cadd_variants, by.x = 'variant_id', by.y = 'variant')
  
  integrated_scores = rbind(snv_variants, aac_variants)
  return(integrated_scores)
}

# gene_phred_scores = split(cadd_variants_2$PHRED, cadd_variants_2$gene)
# boxplot(gene_phred_scores[gene_list])

# ggplot(cadd_variants_2[gene %in% gene_list,PHRED,by=gene], aes(PHRED, fill = gene)) + geom_density(alpha = 0.2)
# ggplot(cadd_variants_2[gene %in% gene_list,PHRED,by=gene], aes(PHRED, fill = gene)) + 
#  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')


# cadd_variants_2[gene %in% gene_list & PHRED >= 20, .(gene, variant_id)]

# lapply(gene_list, function(x) nrow(cesa$variants[gene == x]) - nrow(cadd_variants_2[gene == x]))