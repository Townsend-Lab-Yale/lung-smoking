import pandas as pd
import os

'''Checks all studies to find which ones included multiple samples from the same patient.'''

files = ["luad_oncosg_2020/","luad_broad/", "luad_mskcc_2015/", "lung_msk_2017/", "nsclc_msk_2018/","luad_tsp/","nsclc_tracerx/"]
for i in files:
    #REASONING: If there are multiple samples matching to the same patient, then the patient ID will be repeated for each sample that comes from them.
    df = pd.read_csv(os.path.join("data/clinical_data/",i, "data_clinical_sample.txt"), sep="\t", comment="#")
    if any(df.duplicated(subset = ['PATIENT_ID'])):
        print(i)

tcga_df = pd.read_csv('data/clinical_data/luad_tcga/clinical.tsv', sep = "\t")
if any(tcga_df[::2].duplicated(subset = ['case_submitter_id'])):
    print('tcga')

fmad_df = pd.read_csv('data/luad_FM-AD/clinical.tsv', sep = '\t')
fmad_df = fmad_df[fmad_df['primary_diagnosis'] == 'Adenocarcinoma, NOS']
fmad_df = fmad_df[fmad_df['classification_of_tumor'] == 'primary']
if any(fmad_df.duplicated(subset = ['case_submitter_id'])):
    print('fmad')