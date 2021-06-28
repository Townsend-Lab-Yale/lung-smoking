#!/bin/bash

# Rscript cBioPortal_API.R
## TCGA data is not available in cBioPortal
curl --output-dir '../data/' -o 'luad_tcga.gz' 'https://api.gdc.cancer.gov/data/0458c57f-316c-4a7c-9294-ccd11c97c2f9'
gunzip -c '../data/luad_tcga.gz' > '../data/tsp.luad.maf.txt'
rm '../data/luad_tcga.gz'
