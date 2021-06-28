if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packages = c("cBioPortalData", "AnVIL")

## Check if packages are installed, install them if not, and then load them
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x)
      library(x, character.only = TRUE)
    }
  }
)

luad_broad <- cBioDataPack('luad_broad', ask = FALSE)
lung_msk_2017 <- cBioDataPack('lung_msk_2017', ask = FALSE)
luad_tsp <- cBioDataPack('luad_tsp', ask = FALSE)
luad_oncosg_2020 <- cBioDataPack('luad_oncosg_2020', ask = FALSE)
luad_mskcc_2015 <- cBioDataPack('luad_mskcc_2015', ask = FALSE)
nsclc_pd1_msk_2018 <- cBioDataPack('nsclc_pd1_msk_2018', ask = FALSE)
nsclc_tracerx_2017 <- cBioDataPack('nsclc_tracerx_2017', ask = FALSE)

#to get TCGA data: put this in terminal: curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/0458c57f-316c-4a7c-9294-ccd11c97c2f9'

# study_ids <- c(
#   'luad_broad',
#   'lung_msk_2017',
#   'luad_tsp',
#   'luad_oncosg_2020',
#   'luad_mskcc_2015',
#   'nsclc_pd1_msk_2018',
#   'nsclc_tracerx_2017')
#
# luad_studies = rep(NA, length(study_ids))
#
# for(i in 1:length(study_ids)) {
#   luad_studies[i] <- cBioDataPack(study_ids[i], ask = FALSE)
# }
