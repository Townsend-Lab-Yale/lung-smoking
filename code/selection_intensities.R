library(cancereffectsizeR)
library(data.table)

#PAN-DATA
cesa_total = load_cesa('../data/pan-dataset_samples_cesa.rds')
cesa_total = ces_variant(cesa_total)

gene_selection = split(cesa_total$selection$selection.1$selection_intensity, cesa_total$selection$selection.1$gene)
gene_selection = sapply(gene_selection, function(x){sum(x)})
fwrite(setDT(as.data.frame(gene_selection), keep.rownames = TRUE), '../output/pan_data_selections_no_epistasis.txt')
save_cesa(cesa_total, '../data/pan-dataset_samples_cesa.rds')

#SMOKING
cesa_smoking = load_cesa('../data/smoking_samples_cesa.rds')
cesa_smoking = ces_variant(cesa_smoking)

gene_selection = split(cesa_smoking$selection$selection.1$selection_intensity, cesa_smoking$selection$selection.1$gene)
gene_selection = sapply(gene_selection, function(x){sum(x)})
fwrite(setDT(as.data.frame(gene_selection), keep.rownames = TRUE), '../output/smoking_selections_no_epistasis.txt')
save_cesa(cesa_smoking, '../data/smoking_samples_cesa.rds')

#NONSMOKING
cesa_nonsmoking = load_cesa('../data/nonsmoking_samples_cesa.rds')
cesa_nonsmoking = ces_variant(cesa_nonsmoking)

gene_selection = split(cesa_nonsmoking$selection$selection.1$selection_intensity, cesa_nonsmoking$selection$selection.1$gene)
gene_selection = sapply(gene_selection, function(x){sum(x)})
fwrite(setDT(as.data.frame(gene_selection), keep.rownames = TRUE), '../output/nonsmoking_selections_no_epistasis.txt')
save_cesa(cesa_nonsmoking, '../data/nonsmoking_samples_cesa.rds')




ratioed_selections = as.data.frame(gene_selection / gene_selection['KRAS'])
ratioed_selections = setDT(ratioed_selections, keep.rownames = TRUE)[] 
fwrite(ratioed_selections, '../temp/ratioed_selections.txt')
