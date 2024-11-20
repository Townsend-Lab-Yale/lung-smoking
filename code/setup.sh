#!/bin/bash

echo "Checking for python version..."
python_version=`python --version`
if [ "$python_version" != "Python 3.9.5" ]
then
    echo -e "Warning:\r\n We are working with python version 3.9.5, "`
    `"it has to be installed before running this setup.\r\n\r\n"`
    `"A good way to have multiple python versions on your system is "`
    `"by installing pyenv (https://github.com/pyenv/pyenv).\r\n\r\n"`
    `"If pyenv is installed, run:\r\n\r\n"`
    `"env PYTHON_CONFIGURE_OPTS=\"--enable-shared\" pyenv install 3.9.5\r\n"
    read -p "If you still want to continue with a different python version, "`
    `"press enter to continue or Cntrl-C to abort setup."
fi
echo $python_version
echo ""


echo "Creating a python virtual environment..."
python -m venv .venv
echo "done."
echo ""


echo "Installing all python packages..."
echo ""
.venv/bin/pip install --upgrade pip
.venv/bin/pip install -r requirements.txt
echo "done installing python packages."
echo ""
echo ""


echo "Downloading data..."
echo ""


echo "Downloading GENIE data..."
echo ""
# echo -e "GENIE requires a username and password, "`
#      `"because you have to accept their terms to "`
#      `"download the data set first. If you do not "`
#      `"have done so yet, go to:\r\n"`
#      `"https://www.synapse.org/#!Synapse:syn24179660\r\n"`
#      `"to register."
# source .venv/bin/activate
# python get_synapse_data.py
# deactivate
# echo "...done downloading GENIE data."
# echo ""
# echo ""


cd ../data # keep this line if using synapse (above), it is required
           # for the downloading other files to the proper directory
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
## TODO: Update to latest genie data


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
# ## clinical data file is available from portal.gdc.cancer.gov
# ## https://portal.gdc.cancer.gov/files/e23dae1d-4795-4f94-872a-76da22135a6a
# ## but it does not have a direct link to download (you have to "add to
# ## cart", go to the cart, and then download Clinical TSV). So we
# ## upload it to our website to use a direct URL to download.
# ## The maf file is protected and one has to gain access to obtain it
# ## https://gdc.cancer.gov/access-data/obtaining-access-controlled-data
# ## for now we make it available for ourselves in an easier way
echo "Obtaining main mutations data..."
curl 'https://misc.cidma.us/data/FM-AD_SNV.Bronchus_And_Lung.protected.maf.gz' | gunzip -c > 'data_mutations_extended.txt'
echo "done."
echo ""
echo "Obtainining clinical data file..."
curl 'https://misc.cidma.us/data/clinical.FM-AD_SNV.Bronchus_And_Lung.tar.gz' | gunzip -dc | tar -xf -
echo "done."
echo ""
echo "...done downloading FM-AD data."
echo ""
echo ""
cd ../


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

echo "...done downloading all data."
echo ""
echo ""
cd ../code


echo "Merging MAF and clinical files..."
echo ""
source .venv/bin/activate
python prepare_maf_data.py
python merge_MAF_clinical.py
python store_panel_info.py
deactivate
echo "...done."
echo ""
echo ""


echo "Installing R packages..."
cd variants
mkdir .Rlibs
./ces_installation.R
echo ""
echo "...done installing R packages."
echo ""


echo "Setup ready."
