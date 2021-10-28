library(cancereffectsizeR)
library(data.table)

#PAN-DATA
cesa_total = load_cesa('../data/pan-dataset_samples_cesa.rds')
cesa_total = ces_variant(cesa_total)
save_cesa(cesa_total, '../data/pan-dataset_samples_cesa.rds')

pan_data_si = cesa_total$selection$selection.1[, .(variant_name, gene, selection_intensity, ci_low_95, ci_high_95)]
pan_data_si = pan_data_si[,.(selection_intensity = sum(selection_intensity), ci_low_95 = sum(ci_low_95), ci_high_95 =sum(ci_high_95)),by=.(gene)]
fwrite(pan_data_si, '../output/pan_data_selections_no_epistasis.txt')

#SMOKING
cesa_smoking = load_cesa('../data/smoking_samples_cesa.rds')
cesa_smoking = ces_variant(cesa_smoking)
save_cesa(cesa_smoking, '../data/smoking_samples_cesa.rds')

smoking_si = cesa_smoking$selection$selection.1[, .(variant_name, gene, selection_intensity, ci_low_95, ci_high_95)]
smoking_si = smoking_si[,.(selection_intensity = sum(selection_intensity), ci_low_95 = sum(ci_low_95), ci_high_95 =sum(ci_high_95)),by=.(gene)]
fwrite(smoking_si, '../output/smoking_selections_no_epistasis.txt')


#NONSMOKING
cesa_nonsmoking = load_cesa('../data/nonsmoking_samples_cesa.rds')
cesa_nonsmoking = ces_variant(cesa_nonsmoking)
save_cesa(cesa_nonsmoking, '../data/nonsmoking_samples_cesa.rds')

nonsmoking_si = cesa_nonsmoking$selection$selection.1[, .(variant_name, gene, selection_intensity, ci_low_95, ci_high_95)]
nonsmoking_si = nonsmoking_si[,.(selection_intensity = sum(selection_intensity), ci_low_95 = sum(ci_low_95), ci_high_95 =sum(ci_high_95)),by=.(gene)]
fwrite(nonsmoking_si, '../output/nonsmoking_selections_no_epistasis.txt')





ratioed_selections = as.data.frame(gene_selection / gene_selection['KRAS'])
ratioed_selections = setDT(ratioed_selections, keep.rownames = TRUE)[] 
fwrite(ratioed_selections, '../temp/ratioed_selections.txt')



