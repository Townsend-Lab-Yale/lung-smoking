#!/bin/bash

echo "Downloading data..."
echo ""
cd ../data ## assuming we are running the code from /code

echo "Downloading gencode basic annotation..."
wget 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.basic.annotation.gtf.gz'
gunzip -d gencode.v38lift37.basic.annotation.gtf.gz
echo "...done."
echo ""
echo ""

echo "Downloading hg38 to hg19 chain file..."
wget 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O hg38ToHg19.over.chain.gz
gunzip -d hg38ToHg19.over.chain.gz
echo "...done."
echo ""
echo ""



echo "Downloading GENIE data..."
echo ""
## GENIE does not allow to download previous versions of data, so we
## provide it from our website:
mkdir genie_9
cd genie_9
wget https://misc.cidma.us/data/genie_9/data_clinical_patient.txt
wget https://misc.cidma.us/data/genie_9/data_clinical_sample.txt
wget https://misc.cidma.us/data/genie_9/data_mutations_extended_lifted.txt
wget https://misc.cidma.us/data/genie_9/data_mutations_extended.txt
wget https://misc.cidma.us/data/genie_9/genomic_information.txt
echo "...done downloading GENIE data."
echo ""
echo ""
cd ../


echo "Downloading cBioPortal data..."
echo ""
URL="https://cbioportal-datahub.s3.amazonaws.com/"
EXTENSION=".tar.gz"
declare -a DATASETS=("luad_broad"
                     "lung_msk_2017"
                     "luad_tsp"
                     "luad_oncosg_2020"
                     "luad_mskcc_2015"
                     "nsclc_pd1_msk_2018"
                     "nsclc_tracerx_2017"
                     "lung_nci_2022"
                     "luad_cptac_2020")
for DATASET in ${DATASETS[@]}; do
    echo "Obtaining data set $DATASET..."
    curl "${URL}${DATASET}${EXTENSION}" | gunzip -dc | tar -xf -
    echo "done."
    echo ""
done
echo "Renaming data/luad_broad/data_mutations.txt to "`
     `"data/luad_mskcc_2015/data_mutations_extended.txt for consistency."
mv luad_broad/data_mutations.txt luad_broad/data_mutations_extended.txt
echo "Renaming data/luad_mskcc_2015/data_mutations.txt to "`
     `"data/luad_mskcc_2015/data_mutations_extended.txt for consistency."
mv luad_mskcc_2015/data_mutations.txt luad_mskcc_2015/data_mutations_extended.txt
echo "Renaming data/luad_oncosg_2020/data_mutations.txt to "`
     `"data/luad_oncosg_2020/data_mutations_extended.txt for consistency."
mv luad_oncosg_2020/data_mutations.txt luad_oncosg_2020/data_mutations_extended.txt
echo "Renaming data/lung_msk_2017/data_mutations.txt to "`
     `"data/lung_msk_2017/data_mutations_extended.txt for consistency."
mv lung_msk_2017/data_mutations.txt lung_msk_2017/data_mutations_extended.txt
echo "Renaming data/nsclc_pd1_msk_2018/data_mutations.txt to "`
     `"data/nsclc_pd1_msk_2018/data_mutations_extended.txt for consistency."
mv nsclc_pd1_msk_2018/data_mutations.txt nsclc_pd1_msk_2018/data_mutations_extended.txt
echo "Renaming data/lung_nci_2022/data_mutations.txt to "`
     `"data/lung_nci_2022/data_mutations_extended.txt for consistency."
mv lung_nci_2022/data_mutations.txt lung_nci_2022/data_mutations_extended.txt
echo "Renaming data/nsclc_tracerx_2017/data_mutations.txt to "`
     `"data/nsclc_tracerx_2017/data_mutations_extended.txt for consistency."
mv nsclc_tracerx_2017/data_mutations.txt nsclc_tracerx_2017/data_mutations_extended.txt
echo "Renaming data/luad_tsp/data_mutations.txt to "`
     `"data/luad_tsp/data_mutations_extended.txt for consistency."
mv luad_tsp/data_mutations.txt luad_tsp/data_mutations_extended.txt
echo "Renaming data/luad_cptac_2020/data_mutations.txt to "`
     `"data/luad_cptac_2020/data_mutations_extended.txt for consistency."
mv luad_cptac_2020/data_mutations.txt luad_cptac_2020/data_mutations_extended.txt
echo ""
echo "...done downloading cBioPortal data."
echo ""
echo ""


echo "Downloading Yale data..."
echo ""
mkdir 'yale_luad'
cd 'luad_tcga'
wget https://misc.cidma.us/data/yale_luad/data_clinical_sample.txt
wget https://misc.cidma.us/data/yale_luad/data_mutations_extended.txt
echo "...done downloading Yale data."
echo ""
echo ""
cd ..


echo "Downloading TCGA data..."
echo ""
mkdir 'luad_tcga'
cd 'luad_tcga'
echo "Obtaining main mutations data..."
curl 'https://api.gdc.cancer.gov/data/0458c57f-316c-4a7c-9294-ccd11c97c2f9' | gunzip -c > 'data_mutations_extended.txt'
echo "done."
echo ""
## clinical data file is available from portal.gdc.cancer.gov
## https://portal.gdc.cancer.gov/files/0458c57f-316c-4a7c-9294-ccd11c97c2f9
## but it does not have a direct link to download (you have to "add to
## cart", go to the cart, and then download Clinical TSV). So we
## upload it to our website to use a direct URL to download:
echo "Obtainining clinical data file..."
curl 'https://misc.cidma.us/data/clinical.TCGA.LUAD.tar.gz' | gunzip -dc | tar -xf -
echo "done."
echo ""
echo "...done downloading TCGA data."
echo ""
echo ""
cd ..

echo "Downloading FM-AD data..."
echo ""
mkdir 'luad_fm-ad'
cd 'luad_fm-ad'
echo "FM-AD data is protected. You have to gain access to obtain it."
echo "If you do not have access proceed here: https://gdc.cancer.gov/access-data/obtaining-access-controlled-data"
read -p "If you already have gained access press Y, if not press N (the rest of the code will still work without FM-AD data): " answer

# Check the user's response
if [[ "$answer" == [Yy] ]]; then
    curl 'https://misc.cidma.us/data/FM-AD_SNV.Bronchus_And_Lung.protected.maf.gz' | gunzip -c > 'data_mutations_extended.txt'
    ## clinical data file is available from portal.gdc.cancer.gov
    ## https://portal.gdc.cancer.gov/files/e23dae1d-4795-4f94-872a-76da22135a6a
    ## but it does not have a direct link to download (you have to "add to
    ## cart", go to the cart, and then download Clinical TSV). So we
    ## upload it to our website to use a direct URL to download:
    echo "Obtainining clinical data file..."
    curl 'https://misc.cidma.us/data/clinical.FM-AD_SNV.Bronchus_And_Lung.tar.gz' | gunzip -dc | tar -xf -
    echo "done."
    echo ""
    echo "...done downloading FM-AD data."
else
    # Do not run the command if the user inputs 'N' or 'n'
    echo "Skipping FM-AD data download."
fi
echo ""
echo ""


echo "...done downloading all data."
echo ""
echo ""
