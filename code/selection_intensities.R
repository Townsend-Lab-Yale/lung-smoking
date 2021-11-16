library(cancereffectsizeR)
library(data.table)
library(stringr)

#PAN-DATA
cesa_total = load_cesa('../data/pan-dataset_samples_cesa.rds')
cesa_total = ces_variant(cesa_total)
save_cesa(cesa_total, '../data/pan-dataset_samples_cesa.rds')

# pan_data_si = merge(cesa_total$variants[,.(variant_id,gene)],cesa_total$selection$selection.1[,.(variant_id,selection_intensity,ci_low_95,ci_high_95)],by='variant_id')
# pan_data_si = pan_data_si[,.(selection_intensity = sum(selection_intensity), ci_low_95 = sum(ci_low_95), ci_high_95 =sum(ci_high_95)),by=.(gene)]
# comp_sis = fread('../output/pan_data_selections_no_epistasis.txt')
# comp_df_2 = merge(pan_data_si, comp_sis, by = 'gene')
# comp_df_2[selection_intensity.y > ci_high_95.x | selection_intensity.y < ci_low_95.x]
fwrite(pan_data_si, '../output/pan_data_selections_no_epistasis.txt')

#SMOKING
cesa_smoking = load_cesa('../data/smoking_samples_cesa.rds')
cesa_smoking = ces_variant(cesa_smoking)
save_cesa(cesa_smoking, '../data/smoking_samples_cesa.rds')

# smoking_si = merge(cesa_smoking$variants[,.(variant_id,gene)],cesa_smoking$selection$selection.1[,.(variant_id,selection_intensity,ci_low_95,ci_high_95)],by='variant_id')
# smoking_si = smoking_si[,.(selection_intensity = sum(selection_intensity), ci_low_95 = sum(ci_low_95), ci_high_95 =sum(ci_high_95)),by=.(gene)]
# comp_smoking_sis = fread('../output/smoking_selections_no_epistasis.txt')
# comp_df_smoking_sis = merge(smoking_si, comp_smoking_sis, by = 'gene')
# comp_df_smoking_sis[selection_intensity.y > ci_high_95.x | selection_intensity.y < ci_low_95.x]
fwrite(smoking_si, '../output/smoking_selections_no_epistasis.txt')


#NONSMOKING
cesa_nonsmoking = load_cesa('../data/nonsmoking_samples_cesa.rds')
cesa_nonsmoking = ces_variant(cesa_nonsmoking)
save_cesa(cesa_nonsmoking, '../data/nonsmoking_samples_cesa.rds')

# nonsmoking_si = merge(cesa_nonsmoking$variants[,.(variant_id,gene)],cesa_nonsmoking$selection$selection.1[,.(variant_id,selection_intensity,ci_low_95,ci_high_95)],by='variant_id')
# nonsmoking_si = nonsmoking_si[,.(selection_intensity = sum(selection_intensity), ci_low_95 = sum(ci_low_95), ci_high_95 =sum(ci_high_95)),by=.(gene)]
# comp_nonsmoking_sis = fread('../output/nonsmoking_selections_no_epistasis.txt')
# comp_df_nonsmoking_sis = merge(nonsmoking_si, comp_nonsmoking_sis, by = 'gene')
# comp_df_nonsmoking_sis[selection_intensity.y > ci_high_95.x | selection_intensity.y < ci_low_95.x]
fwrite(nonsmoking_si, '../output/nonsmoking_selections_no_epistasis.txt')





ratioed_selections = as.data.frame(gene_selection / gene_selection['KRAS'])
ratioed_selections = setDT(ratioed_selections, keep.rownames = TRUE)[] 
fwrite(ratioed_selections, '../temp/ratioed_selections.txt')



