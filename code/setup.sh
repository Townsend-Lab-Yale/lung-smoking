#!/bin/bash

cd ../data

## cBioPortal data
URL="https://cbioportal-datahub.s3.amazonaws.com/"
EXTENSION=".tar.gz"

declare -a DATASETS=("luad_broad"
                     "lung_msk_2017"
                     "luad_tsp"
                     "luad_oncosg_2020"
                     "luad_mskcc_2015"
                     "nsclc_pd1_msk_2018"
                     "nsclc_tracerx_2017")

for DATASET in ${DATASETS[@]}; do
   curl "${URL}${DATASET}${EXTENSION}" | gunzip -dc | tar -xf -
done

## TCGA data is not available in cBioPortal
# curl --output-dir '../data/' -o 'luad_tcga.gz' 'https://api.gdc.cancer.gov/data/0458c57f-316c-4a7c-9294-ccd11c97c2f9'
# gunzip -c '../data/luad_tcga.gz' > '../data/tsp.luad.maf.txt'
# rm '../data/luad_tcga.gz'
