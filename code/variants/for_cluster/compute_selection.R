#' This will calculate selection for variants using the entire dataset
#' Note that cancereffectsizeR and ces.refset.hg19 and/or hg38must be installed
#'  prior to running the code (login --> salloc --> module load R; R --> 
#'  installation tutorial from cancereffectsizeR website)

library(cancereffectsizeR)

cesa = load_cesa("../data/pan_data_cesa_for_cancer_epistasis.rds")
cesa = clear_gene_rates(cesa)
cesa = clear_trinuc_rates_and_signatures(cesa)

signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = "LUAD", treatment_naive = F)

cesa = trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                              signature_exclusions = signature_exclusions,
                              cores=30
)

cesa = gene_mutation_rates(cesa, covariates = "lung")

cesa = ces_variant(cesa, 
                   variants = select_variants(cesa, min_freq = 2),
                   cores = 30,
                   run_name = "pan_data")

save_cesa(cesa, "../data/pan_data_cesa_for_cancer_epistasis.rds")