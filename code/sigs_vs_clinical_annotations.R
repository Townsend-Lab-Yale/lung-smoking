smoking_sig_vs_clinic = unique(maf_clinical, by='Sample ID')[`Sample ID` %in% smoking_samples, .(`Sample ID`,Source,Smoker)]
nonsmoking_sig_vs_clinic = unique(maf_clinical, by='Sample ID')[`Sample ID` %in% nonsmoking_samples, .(`Sample ID`,Source,Smoker)]

smoking_matches = nrow(smoking_sig_vs_clinic[Smoker == T])/length(smoking_samples)
nonsmoking_matches = nrow(nonsmoking_sig_vs_clinic[Smoker == F])/length(nonsmoking_samples)

smoking_disagreements_by_source = round(table(smoking_sig_vs_clinic[Smoker == F,Source])/table(smoking_sig_vs_clinic$Source),3)
nonsmoking_disagreements_by_source = round(table(nonsmoking_sig_vs_clinic[Smoker == T,Source])/table(nonsmoking_sig_vs_clinic$Source),3)
