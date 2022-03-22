library(cancereffectsizeR)
library(data.table)

#' Produces a vcf from cesa's maf
#' Zip and upload resultant vcf to https://cadd.gs.washington.edu/score 

cesa = load_cesa('../data/pan-dataset_samples_cesa.rds')

genes_list = fread('../../data/genes_list.txt',header=F)$V1
vcf = cesa$maf[top_gene %in% genes_list & variant_type %in% c('snv','del','ins'), .(Chromosome, Start_Position, variant_id, Reference_Allele, Tumor_Allele, variant_type)]#, top_gene, top_consequence)]
names(vcf) = c("#CHROM","POS","ID","REF","ALT",'TYPE')
fwrite(vcf, '../../temp/CADD_variants.vcf',sep="\t")
