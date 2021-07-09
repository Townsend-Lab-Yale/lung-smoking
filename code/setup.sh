#!/bin/bash

## * Install python requirements
python3.9 -m venv .venv
.venv/bin/pip install -r requirements.txt


## * Obtain gene coordinates
source .venv/bin/activate
python obtain_gene_coordinates.py
deactivate


## * Download data
cd ../data


## ** cBioPortal data
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


## ** TCGA data
mkdir 'luad_tcga'
cd 'luad_tcga'
curl 'https://api.gdc.cancer.gov/data/0458c57f-316c-4a7c-9294-ccd11c97c2f9' | gunzip -c > 'data_mutations_extended.txt'
## clinical data file is available from portal.gdc.cancer.gov
## https://portal.gdc.cancer.gov/files/0458c57f-316c-4a7c-9294-ccd11c97c2f9
## but it does not have a direct link to download (you have to "add to
## cart", go to the cart, and then download Clinical TSV). So we
## upload it to our website to use a direct URL to download:
curl 'https://misc.cidma.us/data/clinical.TCGA.LUAD.tar.gz' | gunzip -dc | tar -xf -
cd ..


## ** FM-AD TCGA
mkdir 'luad_FM-AD'
cd 'luad_FM-AD'
## clinical data file is available from portal.gdc.cancer.gov
## https://portal.gdc.cancer.gov/files/e23dae1d-4795-4f94-872a-76da22135a6a
## but it does not have a direct link to download (you have to "add to
## cart", go to the cart, and then download Clinical TSV). So we
## upload it to our website to use a direct URL to download:
curl 'https://misc.cidma.us/data/clinical.FM-AD_SNV.Bronchus_And_Lung.tar.gz' | gunzip -c | tar -xf -
## The maf file is protected and one has to gain access to obtain it
## https://gdc.cancer.gov/access-data/obtaining-access-controlled-data
## for now we make it available for ourselves in an easier way
curl 'https://misc.cidma.us/data/FM-AD_SNV.Bronchus_And_Lung.protected.maf.gz' | gunzip -c > 'data_mutations_extended.txt'
cd ..
