library(cancereffectsizeR)
library(data.table)

#' Produces a vcf from cesa's maf
#' Zip and upload resultant vcf to https://cadd.gs.washington.edu/score 

cesa = load_cesa('../data/pan-dataset_samples_cesa.rds')

genes_list = fread('../../data/genes_list.txt',header=F)$V1
variants = cesa$maf[genes %in% genes_list & variant_type %in% c('snv','del','ins'), .(Chromosome, Start_Position, variant_id, Reference_Allele, Tumor_Allele)]#, variant_type, top_gene, top_consequence)]
names(variants) = c("#CHROM","POS","ID","REF","ALT")
fwrite(variants, '../temp/CADD_variants.vcf',sep="\t")
