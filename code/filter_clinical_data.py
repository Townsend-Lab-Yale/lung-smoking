import os
import pandas as pd
import numpy as np

from locations import merged_clinical_file_name
from locations import location_data

random_seeds = [777, 778]
"""Random seed to feed the random generators, to be able to replicate
results."""



files = ["luad_oncosg_2020","luad_broad", "luad_mskcc_2015", "lung_msk_2017", "nsclc_pd1_msk_2018",
         "nsclc_tracerx_2017","genie_9","lung_nci_2022","luad_cptac_2020"]
#files that can be auto-merged
nm_files = ["tsp.luad.maf.txt"]
#files that don't need merging

clinical_files = {i:None for i in files}

for i in files:
    df1 = pd.read_csv(os.path.join(location_data,i, "data_clinical_patient.txt"), sep="\t", comment="#")
    df2 = pd.read_csv(os.path.join(location_data,i, "data_clinical_sample.txt"), sep="\t", comment="#")
    clinical_files[i] = pd.merge(df1, df2, on="PATIENT_ID")
#merging on basis of patient id

s_files = ["luad_tcga"]
#files that need special merging procedure
df1 = pd.read_csv(os.path.join(location_data,'luad_tcga/clinical.tsv'), sep="\t")
df2 = pd.read_csv(os.path.join(location_data,'luad_tcga/exposure.tsv'), sep="\t")
clinical_files['luad_tcga'] = pd.merge(df1, df2, on="case_id")
#merging on basis of patient ID

clinical_files['yale_luad'] = pd.read_csv(os.path.join(location_data,'yale_luad/data_clinical_sample.txt'), sep="\t") \
                                .rename(columns={'patient_name':'Patient ID','tumor_name':'Sample ID'})

stage_dict = {'I':'1','IA':'1a','IB':'1b','II':'2','IIA':'2a','IIB':'2b','III':'3','IIIA':'3a','IIIB':'3b','IV':'4'}

'''
COLUMN SUMMARIES:
Sample ID: ID of tumor sample
Smoker: True if ever smoked, False if never smoked
Stage: Tumor stage
Survivals: self-explanatory
Treatment: True if received targeted therapy treatment, False otherwise
'''

"""filtering requires individual steps because each file is too different"""

broad_df = clinical_files.get('luad_broad')
#converting smoking to 1/0
broad_df["SMOKER"] = broad_df["SMOKER"].apply(lambda x: True if x in ("Heavy Smoker", "Light Smoker") else False if x == "Never Smoker" else np.NaN)
#converting stages to one format amongst all datasets
broad_df["STAGE"] = broad_df["STAGE"].map(stage_dict).fillna(broad_df["STAGE"])
broad_df = broad_df[["PATIENT_ID", "SMOKER", "STAGE", "PFS_MONTHS"]]
broad_df.columns = ["Sample ID","Smoker","Stage","Progression Free Survival (months)"]
#"All tumors were chemotherapy-naive, primary resection specimens except for one case with whole genome sequence data (LU-A08-43) that was a post-chemotherapy metastatic tumor from a never-smoker" (Broad, 2012).
broad_df.loc[:, 'Treatment'] = False
#removing non-primary samples
broad_df = broad_df.drop(broad_df.index[broad_df['Sample ID'] == 'LU-A08-43'])


tcga_df = clinical_files.get('luad_tcga')
#replacing '-- with -1 should a temporary fix, reason for it is that NaN and None can't be used for >/< comparisons
tcga_df = tcga_df.replace({"'--":-1, "not reported":''})
tcga_df["tumor_stage"] = tcga_df["tumor_stage"].apply(lambda x: x[6:])
#convert smoking to 1/0. multiple columns for smoking data (pack years, years, cigarettes per day), but pack years provided the least NA ('--) values
#unfortunately, TCGA did not specify when 0 packs were smoked, leaving us unsure when '-- refers to NA or to never smoked. We will assume it means never smoked.
tcga_df["pack_years_smoked"] = tcga_df["pack_years_smoked"].apply(lambda x: True if float(x) > 0 else False)
#converting days to death to months to death for comparison with PFS_months and OS_months
tcga_df["months_to_death"] = tcga_df["days_to_death"].apply(lambda x: float(x)/30 if float(x) >= 0 else np.NaN)
#converting stages to one format amongst all datasets
tcga_df["tumor_stage"] = tcga_df["tumor_stage"].str.upper().map(stage_dict).fillna(tcga_df["tumor_stage"])
#removing duplicate rows for radiation therapy
tcga_df = tcga_df[tcga_df['treatment_type'] == 'Pharmaceutical Therapy, NOS']
#adding treatment column
tcga_df['Treatment'] = tcga_df['prior_treatment'].apply(lambda x: True if x == 'Yes' else False if x == 'No' else np.NaN)
tcga_df = tcga_df[["case_id","pack_years_smoked","tumor_stage", "months_to_death", "Treatment"]]
tcga_df.columns = ["Sample ID","Smoker","Stage","Overall Survival (months)","Treatment"]
#metastatis_at_diagnosis/metastasis_at_diagnosis_site columns are only NA values


oncosg_df = clinical_files.get('luad_oncosg_2020')
#converting smoking to 1/0
oncosg_df["SMOKING_STATUS"] = oncosg_df["SMOKING_STATUS"].apply(lambda x: True if x == "Yes" else False if x == "No" else np.NaN)
#converting stages to one format amongst all datasets
oncosg_df["STAGE"] = oncosg_df["STAGE"].map(stage_dict).fillna(oncosg_df["STAGE"])
#adding treatment column
oncosg_df['Treatment'] = oncosg_df.apply(lambda x: True if (x.TKI_TREATMENT == "Yes" or x.CHEMOTHERAPY == "Yes") else False if (x.TKI_TREATMENT == "No" and x.CHEMOTHERAPY == "No") else np.NaN, axis=1)
oncosg_df['Treatment'] = oncosg_df['TKI_TREATMENT'].apply(lambda x: True if x == 'Yes' else False if x == 'No' else np.NaN)
oncosg_df = oncosg_df[["PATIENT_ID", "SMOKING_STATUS", "STAGE", "OS_MONTHS", "Treatment"]]
oncosg_df.columns = ["Sample ID","Smoker","Stage","Overall Survival (months)", "Treatment"]
#no indication in paper of whether tumors were primary or metastatic


msk2015_df = clinical_files.get('luad_mskcc_2015')
msk2015_df["SMOKING_HISTORY"] = msk2015_df["SMOKING_HISTORY"].apply(lambda x: True if x in ("Current", "Former") else False if x == "Never" else np.NaN)
# According to supplementary methods, all patients had stage 4 cancer
msk2015_df['Stage'] = "4"
msk2015_df['Treatment'] = True
# msk2015_df['Treatment'] = False
# msk2015_df['Treatment'][msk2015_df['Sample ID'] == 'DM123062'] = True
non_luad_sample_ids_msk2015 = msk2015_df[msk2015_df['HISTOLOGY'] != 'Adenocarcinoma']['SAMPLE_ID']
msk2015_df = msk2015_df[msk2015_df['HISTOLOGY'] == 'Adenocarcinoma']
# keeping only primary samples from lung tissue, according to supplementary table 4
# samples LO3793, SC6470, VA1330, VA7859, WA7899 are from lung tissue but are metastatic
msk2015_df = msk2015_df[msk2015_df["SAMPLE_ID"].isin(["AU5884","CA9903","GR0134","JB112852","LO5004","MA7027","NI9507","RO3338","SB010944","SC0899","SR070761","TU0428"])]
msk2015_df = msk2015_df[["SAMPLE_ID", "Stage", "SMOKING_HISTORY", "PFS_MONTHS", "Treatment"]]
#using sample_id instead of patient_id because that is what matches with patient_id in the maf file
msk2015_df.columns = ["Sample ID","Stage","Smoker","Progression Free Survival (months)","Treatment"]


msk2017_df = clinical_files.get('lung_msk_2017')
# Excluding "Former light" from consideration
light_smoker_sample_ids_2017 = msk2017_df[msk2017_df['SMOKING_HISTORY'] == 'Former light']['SAMPLE_ID']
msk2017_df["SMOKING_HISTORY"] = msk2017_df["SMOKING_HISTORY"].apply(lambda x: True if x in ("Current heavy", "Former heavy") else False if x == "Never" else np.NaN)
#converting stages to one format amongst all datasets
msk2017_df["STAGE_AT_DIAGNOSIS"] = msk2017_df["STAGE_AT_DIAGNOSIS"].map(stage_dict | {'IA L,IV R':'1a'}).fillna(msk2017_df["STAGE_AT_DIAGNOSIS"])
#adding treatment column
msk2017_df['Treatment'] = msk2017_df['SAMPLE_PRE_LUNG_THERAPY'].apply(lambda x: False if x == 'YES' else True if x == 'NO' else np.NaN)
#removing non-primary samples
metastatic_sample_ids_2017 = msk2017_df[msk2017_df['SAMPLE_TYPE'] != 'Primary']['SAMPLE_ID']
msk2017_df = msk2017_df[msk2017_df['SAMPLE_TYPE'] == 'Primary']
#removing multiple samples from same patient in msk 2017 data, keeping most recent and then most coverage.
msk2017_df = msk2017_df.sort_values(by = ['SAMPLE_TYPE','LINES_OF_TX_PRIOR_IMPACT','TUMOR_PURITY'], ascending=[False, True, False])
multi_sample_ids_2017 = msk2017_df[msk2017_df.duplicated(subset = ['PATIENT_ID'], keep = 'first')]['SAMPLE_ID']
msk2017_df = msk2017_df.drop_duplicates(subset = ['PATIENT_ID'], keep = 'first')
msk2017_df = msk2017_df[["PATIENT_ID","SAMPLE_ID", "SMOKING_HISTORY", "STAGE_AT_DIAGNOSIS", "VITAL_STATUS","Treatment"]]
#using sample_id instead of patient_id because that is what matches with patient_id in the maf file
#msk2017 only has vital status rather than PFS or OS.
msk2017_df.columns = ["Patient ID","Sample ID","Smoker","Stage","Vital Status","Treatment"]

msk2018_df = clinical_files.get('nsclc_pd1_msk_2018')
msk2018_df['SMOKER'] = np.where(msk2018_df['SAMPLE_ID'].isin(light_smoker_sample_ids_2017),np.NaN,msk2018_df['SMOKER'])
msk2018_df['SMOKER'] = msk2018_df['SMOKER'].apply(lambda x: True if x == 'Ever' else False if x == 'Never' else np.NaN)
msk2018_df['Stage'] = np.NaN
#all patients treated with ICI
msk2018_df['Treatment'] = True
non_luad_sample_ids_msk2018 = msk2018_df[msk2018_df['CANCER_TYPE_DETAILED'] != 'Lung Adenocarcinoma']['SAMPLE_ID']
msk2018_df = msk2018_df[msk2018_df['CANCER_TYPE_DETAILED'] == 'Lung Adenocarcinoma']
msk2018_df = msk2018_df[['PATIENT_ID','SAMPLE_ID','SMOKER','Stage','PFS_MONTHS','Treatment']]
msk2018_df.columns = ["Patient ID","Sample ID","Smoker","Stage","Progression Free Survival (months)","Treatment"]

tracer_df = clinical_files.get('nsclc_tracerx_2017')
tracer_df['SMOKING_HISTORY'] = tracer_df['SMOKING_HISTORY'].apply(lambda x: True if x in ("Current Smoker", "Ex-Smoker", "Recent Ex-Smoker") else False if x == "Never Smoked" else np.NaN)
tracer_df['TUMOR_STAGE'] = tracer_df['TUMOR_STAGE'].map(stage_dict).fillna(tracer_df['TUMOR_STAGE'])
tracer_df['Treatment'] = tracer_df['SAMPLE_COLLECTION_TIMEPOINT'].apply(lambda x: True if x == 'Post-treatment' else False if x == 'Pre-treatment' else np.NaN)
non_luad_sample_ids_tracer = tracer_df[tracer_df['HISTOLOGY'] != 'Invasive adenocarcinoma']
tracer_df = tracer_df[tracer_df['HISTOLOGY'] == 'Invasive adenocarcinoma']
#removing non-primary samples
metastatic_sample_ids_tracer = tracer_df[tracer_df['SAMPLE_TYPE'] != 'Primary']['SAMPLE_ID']
tracer_df = tracer_df[tracer_df['SAMPLE_TYPE'] == 'Primary']
#removing cfDNA samples
tracer_df = tracer_df[~(tracer_df['SAMPLE_ID'].str.contains('DNA'))]
#randomly sampling from multiple samples per patient (multi-region sampling)
tracer_df_sampled = pd.DataFrame()
for id in pd.unique(tracer_df['PATIENT_ID']):
    tracer_df_sampled = pd.concat([tracer_df_sampled, tracer_df[tracer_df['PATIENT_ID'] == id].sample(random_state=random_seeds[0])])
tracer_df_sampled = tracer_df_sampled[['SAMPLE_ID','SMOKING_HISTORY','TUMOR_STAGE','RFS_MONTHS','Treatment']]
#not sure if regression free survival = progression free survival
tracer_df_sampled.columns = ["Sample ID","Smoker","Stage","Progression Free Survival (months)","Treatment"]
keep_tracer_samples = tracer_df_sampled['Sample ID']


genie_df = clinical_files.get('genie_9')
genie_df['Smoker'] = np.NaN
genie_df['Stage'] = np.NaN
genie_df['Treatment'] = np.NaN
#survival column is also np.NaN (will automatically fill in when dfs are concatenated)
#removing non-LUAD samples and collecting indices to remove from maf file
non_luad_sample_ids_genie = genie_df[genie_df['CANCER_TYPE_DETAILED'] != 'Lung Adenocarcinoma']['SAMPLE_ID']
genie_df = genie_df[genie_df['CANCER_TYPE_DETAILED'] == 'Lung Adenocarcinoma']
#removing non-primary samples
metastatic_sample_ids_genie = genie_df[genie_df['SAMPLE_TYPE'] != 'Primary']['SAMPLE_ID']
genie_df = genie_df[genie_df['SAMPLE_TYPE'] == 'Primary']
#TEMPORARY SOLUTION to remove multiple samples for one patient
genie_df = genie_df.sample(frac=1, random_state=random_seeds[1])
genie_df = genie_df.sort_values(by = ['AGE_AT_SEQ_REPORT'])
multi_sample_ids_genie = genie_df[genie_df.duplicated(subset = ['PATIENT_ID'], keep = 'first')]['SAMPLE_ID']
genie_df = genie_df.drop_duplicates(subset = ['PATIENT_ID'], keep = 'first')
#this is only for preserving MSK patient IDs and will likely interfere with other patient IDs, but only Sample ID matters in the end
genie_df['Patient ID'] = genie_df['PATIENT_ID'].str.slice(10)
genie_df = genie_df[['Patient ID','SAMPLE_ID','Smoker','Stage','Treatment']]
genie_df.columns = ['Patient ID','Sample ID','Smoker','Stage','Treatment']


nci_df = clinical_files.get('lung_nci_2022')
non_luad_sample_ids_nci = nci_df[nci_df['HISTOLOGY'] != 'Adenocarcinomas']['SAMPLE_ID']
nci_df = nci_df[nci_df["HISTOLOGY"] == "Adenocarcinomas"]
nci_df.loc[:, 'Smoker'] = False
# all samples were of primary tumor
nci_df.loc[:, 'Stage'] = nci_df["TUMOR_STAGE"].map(stage_dict).fillna(nci_df["TUMOR_STAGE"])
nci_df.loc[:, 'Overall Survival (Months)'] = nci_df["OS_MONTHS"].fillna("Living")
nci_df.loc[:, 'Treatment'] = False
nci_df = nci_df[["PATIENT_ID",'SAMPLE_ID','Smoker','Stage','Overall Survival (Months)','Treatment']]
nci_df.columns = ['Patient ID','Sample ID','Smoker','Stage','Overall Survival (Months)','Treatment']


cptac_df = clinical_files.get('luad_cptac_2020')
cptac_df['SMOKING_STATUS'] = cptac_df['SMOKING_STATUS'].apply(lambda x: True if x == 'Smoker' else False if x=='Non-Smoker' else np.NaN)
cptac_df['STAGE'] = cptac_df['STAGE'].apply(lambda x: str(x).lower())
cptac_df['Treatment'] = False
cptac_df = cptac_df[['PATIENT_ID','SAMPLE_ID','SMOKING_STATUS','STAGE','Treatment']]
cptac_df.columns = ['Patient ID','Sample ID','Smoker','Stage','Treatment']


fmad_df = pd.read_csv(os.path.join(location_data,'luad_fm-ad/clinical.tsv'), sep = '\t')
non_luad_sample_ids_fmad = fmad_df[fmad_df['primary_diagnosis'] != 'Adenocarcinoma, NOS']['case_id']
fmad_df = fmad_df[fmad_df['primary_diagnosis'] == 'Adenocarcinoma, NOS']
fmad_non_primary = fmad_df[fmad_df['classification_of_tumor'] != 'primary']['case_id']
fmad_df = fmad_df[fmad_df['classification_of_tumor'] == 'primary']
fmad_df = fmad_df[['case_id']]
fmad_df.columns = ['Sample ID']

yale_df = clinical_files.get('yale_luad')
non_yale_sample_ids = yale_df['Sample ID'][~yale_df['Sample ID'].str.startswith('MDA')]
yale_df = yale_df[yale_df['Sample ID'].str.startswith('MDA')]
yale_df['Smoker'] = np.NaN
yale_df['Stage'] = np.NaN
yale_df['Treatment'] = np.NaN
yale_df = yale_df[['Sample ID','Smoker','Stage','Treatment']]

#removing repeated patients between msk 2017 and msk 2018, keeping msk 2018
merged_temp = pd.merge(msk2017_df, msk2018_df, on = 'Patient ID', how = 'inner')
dup_patient_ids_1718 = merged_temp['Patient ID']
dup_sample_ids_1718 = msk2017_df[msk2017_df['Patient ID'].isin(dup_patient_ids_1718)]['Sample ID']
msk2017_df = msk2017_df[~msk2017_df['Patient ID'].isin(dup_patient_ids_1718)]

#removing repeated patients between genie and msk 2018, keeping msk 2018
merged_temp = pd.merge(genie_df, msk2018_df, on = 'Patient ID', how = 'inner')
dup_patient_ids_gen18 = merged_temp['Patient ID']
dup_sample_ids_gen18 = genie_df[genie_df['Patient ID'].isin(dup_patient_ids_gen18)]['Sample ID']
genie_df = genie_df[~genie_df['Patient ID'].isin(dup_patient_ids_gen18)]

#removing repeated patients between genie and msk 2017, keeping msk 2017
merged_temp = pd.merge(genie_df, msk2017_df, on = 'Patient ID', how = 'inner')
dup_patient_ids_gen17 = merged_temp['Patient ID']
dup_sample_ids_gen17 = genie_df[genie_df['Patient ID'].isin(dup_patient_ids_gen17)]['Sample ID']
genie_df = genie_df[~genie_df['Patient ID'].isin(dup_patient_ids_gen17)]

genie_df = genie_df.drop(columns='Patient ID')
msk2017_df = msk2017_df.drop(columns='Patient ID')
msk2018_df = msk2018_df.drop(columns='Patient ID')


"""TSP not included because desired information is not in dataset and is not currently accessible"""

all_clinical_df = pd.concat([broad_df, tcga_df, oncosg_df, msk2015_df, msk2017_df, msk2018_df, 
                             tracer_df_sampled, genie_df, fmad_df, nci_df, cptac_df, yale_df])

all_clinical_df.to_csv(merged_clinical_file_name)
