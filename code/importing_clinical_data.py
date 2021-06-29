import os
import pandas as pd
import numpy as np

if '__file__' not in globals():
    __file__ = '.'

location_data = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "data/"
                 "clinical_data"))

files = ["luad_oncosg_2020","luad_broad", "luad_mskcc_2015", "lung_msk_2017", "nsclc_msk_2018", "nsclc_tracerx"]
#files that can be auto-merged
nm_files = ["tsp.luad.maf.txt"]
#files that don't need merging


for i in files:
    df1 = pd.read_csv(os.path.join(location_data,i, "data_clinical_patient.txt"), sep="\t", comment="#")
    df2 = pd.read_csv(os.path.join(location_data,i, "data_clinical_sample.txt"), sep="\t", comment="#")
    merged_df = pd.merge(df1, df2, on="PATIENT_ID")
    merged_df.to_csv("output/unmerged_clinical/" + i + "_clinical.txt")
#merging on basis of patient id

s_files = ["luad_tcga"]
#files that need special merging procedure
df1 = pd.read_csv("data/clinical_data/luad_tcga/clinical.tsv", sep="\t")
df2 = pd.read_csv("data/clinical_data/luad_tcga/exposure.tsv", sep="\t")
merged_df = pd.merge(df1, df2, on="case_id")
merged_df.to_csv("output/unmerged_clinical/luad_tcga_clinical.txt")
#merging on basis of patient ID

stage_dict = {'I':'1','IA':'1a','IB':'1b','II':'2','IIA':'2a','IIB':'2b','III':'3','IIIA':'3a','IIIB':'3b','IV':'4'}


'''
COLUMN SUMMARIES:
Sample ID: ID of tumor sample
Smoker: True if ever smoked, False if never smoked
Stage: Tumor stage
Survivals: self-explanatory
Treatment: True if received targeted therapy treatment, False otherwise
Metastatic: True if tumor is metastatic, False otherwise
is_LUAD: True if tumor histology indicates lung adenocarcinoma, False otherwise
'''

"""filtering requires individual steps because each file is too different"""

broad_df = pd.read_csv("output/unmerged_clinical/luad_broad_clinical.txt")
#converting smoking to 1/0
broad_df["SMOKER"] = broad_df["SMOKER"].apply(lambda x: True if x in ("Heavy Smoker", "Light Smoker") else False if x == "Never Smoker" else np.NaN)
#converting stages to one format amongst all datasets
broad_df["STAGE"] = broad_df["STAGE"].map(stage_dict).fillna(broad_df["STAGE"])
broad_df = broad_df[["PATIENT_ID", "SMOKER", "STAGE", "PFS_MONTHS"]]
broad_df.columns = ["Sample ID","Smoker","Stage","Progression Free Survival (months)"]
#"All tumors were chemotherapy-naive, primary resection specimens except for one case with whole genome sequence data (LU-A08-43) that was a post-chemotherapy metastatic tumor from a never-smoker" (Broad, 2012).
broad_df['Treatment'] = False
broad_df['Metastatic'] = False
broad_df.loc[broad_df.index[broad_df['Sample ID'] == 'LU-A08-43'],['Treatment', 'Metastatic']] = [True,True]
broad_df['is_LUAD'] = True
broad_df.to_csv("output/unmerged_clinical/luad_broad_clinical.txt")



tcga_df = pd.read_csv("output/unmerged_clinical/luad_tcga_clinical.txt")
#replacing '-- with -1 should a temporary fix, reason for it is that NaN and None can't be used for >/< comparisons
tcga_df = tcga_df.replace({"'--":-1, "not reported":''})
tcga_df["tumor_stage"] = tcga_df["tumor_stage"].apply(lambda x: x[6:])
#convert smoking to 1/0. multiple columns for smoking data (pack years, years, cigarettes per day), but pack years provided the least NA ('--) values
#unfortunately, TCGA did not specify when 0 packs were smoked, leaving us unsure when '-- refers to NA or to never smoked. We will assume it means never smoked.
tcga_df["pack_years_smoked"] = tcga_df["pack_years_smoked"].apply(lambda x: True if float(x) > 0 else 0)
#converting days to death to months to death for comparison with PFS_months and OS_months
tcga_df["months_to_death"] = tcga_df["days_to_death"].apply(lambda x: float(x)/30 if float(x) >= 0 else np.NaN)
#converting stages to one format amongst all datasets
tcga_df["tumor_stage"] = tcga_df["tumor_stage"].str.upper().map(stage_dict).fillna(tcga_df["tumor_stage"])
#removing duplicate rows for radiation therapy
tcga_df = tcga_df[tcga_df['treatment_type'] == 'Pharmaceutical Therapy, NOS']
#adding treatment column
tcga_df['Treatment'] = tcga_df['treatment_or_therapy'].apply(lambda x: True if x == 'yes' else False if x == 'no' else np.NaN)
'''
#tcga_df['treatment_or_therapy'] = tcga_df['treatment_or_therapy'].apply(lambda x: 1 if one is yes else 0 if both no else np.NaN)
tcga_df['Treatment'] = -1.0
for row in tcga_df.itertuples():
    if row.Index %2 == 0 and row.case_id == tcga_df.loc[row.Index + 1, 'case_id']:
        if row.treatment_or_therapy == 'yes' or tcga_df.loc[row.Index + 1, 'treatment_or_therapy'] == 'yes':
            tcga_df.at[row.Index,'Treatment'] = 1
        elif row.treatment_or_therapy == 'no' and tcga_df.loc[row.Index + 1, 'treatment_or_therapy'] == 'no':
            tcga_df.at[row.Index,'Treatment'] = 0
        else:
            tcga_df.at[row.Index,'Treatment'] = np.NaN
#remove duplicate rows
tcga_df = tcga_df[tcga_df['Treatment'] != -1]
'''
tcga_df = tcga_df[["case_id","pack_years_smoked","tumor_stage", "months_to_death", "Treatment"]]
tcga_df.columns = ["Sample ID","Smoker","Stage","Overall Survival (months)","Treatment"]
#metastatis_at_diagnosis/metastasis_at_diagnosis_site columns are only NA values
tcga_df["Metastatic"] = np.NaN
tcga_df['is_LUAD'] = True
tcga_df.to_csv("output/unmerged_clinical/luad_tcga_clinical.txt")



oncosg_df = pd.read_csv("output/unmerged_clinical/luad_oncosg_2020_clinical.txt")
#converting smoking to 1/0
oncosg_df["SMOKING_STATUS"] = oncosg_df["SMOKING_STATUS"].apply(lambda x: True if x == "Yes" else False if x == "No" else np.NaN)
#converting stages to one format amongst all datasets
oncosg_df["STAGE"] = oncosg_df["STAGE"].map(stage_dict).fillna(oncosg_df["STAGE"])
#adding treatment column
oncosg_df['Treatment'] = oncosg_df['TKI_TREATMENT'].apply(lambda x: True if x == 'Yes' else False if x == 'No' else np.NaN)
'''
oncosg_df['Treatment'] = np.NaN
for row in oncosg_df.itertuples():
    if row.TKI_TREATMENT == 'Yes' or row.CHEMOTHERAPY == 'Yes':
        oncosg_df.at[row.Index,'Treatment'] = 1
    elif row.TKI_TREATMENT == 'No' and row.CHEMOTHERAPY == 'No':
        oncosg_df.at[row.Index,'Treatment'] = 0
'''
'''
conditions = [
    (oncosg_df['TKI_TREATMENT'] == 'Yes' or oncosg_df['CHEMOTHERAPY'] == 'Yes'), 
    (oncosg_df['TKI_TREATMENT'] == 'No' and oncosg_df['CHEMOTHERAPY'] == 'No'),
]
values = [1, 0]
oncosg_df['Treatment'] = np.select(conditions, values, default = np.NaN)
'''
oncosg_df = oncosg_df[["PATIENT_ID", "SMOKING_STATUS", "STAGE", "OS_MONTHS", "Treatment"]]
oncosg_df.columns = ["Sample ID","Smoker","Stage","Overall Survival (months)", "Treatment"]
#no indication in paper of whether tumors were primary or metastatic
oncosg_df['Metastatic'] = np.NaN
oncosg_df['is_LUAD'] = True
oncosg_df.to_csv("output/unmerged_clinical/luad_oncosg_2020_clinical.txt")



msk2015_df = pd.read_csv("output/unmerged_clinical/luad_mskcc_2015_clinical.txt")
#Removing unspecified NSCLC tumor subtypes and LUSC tumors.
#msk2015_df = msk2015_df[msk2015_df['HISTOLOGY'] == "Adenocarcinoma"]
msk2015_df["SMOKING_HISTORY"] = msk2015_df["SMOKING_HISTORY"].apply(lambda x: True if x in ("Current", "Former") else False if x == "Never" else np.NaN)
#msk2015 does not have stage
msk2015_df['Stage'] = np.NaN
#All patients treated with pembrolizumab, no indication of primary/metastasis in paper
msk2015_df['Treatment'] = True
msk2015_df['Metastatic'] = np.NaN
msk2015_df['is_LUAD'] = msk2015_df['HISTOLOGY'].apply(lambda x: True if x == 'Adenocarcinoma' else False)
msk2015_df = msk2015_df[["SAMPLE_ID", "Stage", "SMOKING_HISTORY", "PFS_MONTHS", "Treatment", "Metastatic", "is_LUAD"]]
#using sample_id instead of patient_id because that is what matches with patient_id in the maf file
msk2015_df.columns = ["Sample ID","Stage","Smoker","Progression Free Survival (months)","Treatment", "Metastatic", "is_LUAD"]
msk2015_df.to_csv("output/unmerged_clinical/luad_mskcc_2015_clinical.txt")


msk2017_df = pd.read_csv("output/unmerged_clinical/lung_msk_2017_clinical.txt")
msk2017_df["SMOKING_HISTORY"] = msk2017_df["SMOKING_HISTORY"].apply(lambda x: True if x in ("Current heavy", "Former heavy", "Former light") else False if x == "Never" else np.NaN)
#converting stages to one format amongst all datasets
msk2017_df["STAGE_AT_DIAGNOSIS"] = msk2017_df["STAGE_AT_DIAGNOSIS"].map(stage_dict | {'IA L,IV R':'1a'}).fillna(msk2017_df["STAGE_AT_DIAGNOSIS"])
#adding treatment column
msk2017_df['Treatment'] = msk2017_df['TARGET_THERAPY'].apply(lambda x: True if x == 'YES' else False if x == 'NO' else np.NaN)
'''
msk2017_df['Treatment'] = np.NaN
for row in msk2017_df.itertuples():
    if row.TARGET_THERAPY == 'YES' or row.IMMUNE_TREATMENT == 'YES' or row.TREATMENT_CHEMOTHERAPY == 'YES':
        msk2017_df.at[row.Index,'Treatment'] = 1
    elif row.TARGET_THERAPY == 'NO' and row.IMMUNE_TREATMENT == 'NO' and row.TREATMENT_CHEMOTHERAPY == 'NO':
        msk2017_df.at[row.Index,'Treatment'] = 0
'''
'''
conditions = [
    (msk2017_df['TARGET_THERAPY'] == 'YES' or msk2017_df['IMMUNE_TREATMENT'] == 'YES' or msk2017_df['TREATMENT_CHEMOTHERAPY'] == 'YES'), 
    (msk2017_df['TARGET_THERAPY'] == 'NO' and msk2017_df['IMMUNE_TREATMENT'] == 'NO' and msk2017_df['TREATMENT_CHEMOTHERAPY'] == 'NO'),
]
values = [1, 0]
msk2017_df['Treatment'] = np.select(conditions, values, default = np.NaN)
'''
#converting primary/metastasis column to 1/0
msk2017_df["SAMPLE_TYPE"] = msk2017_df["SAMPLE_TYPE"].apply(lambda x: True if x == 'Metastasis' else 0 if False == 'Primary' else np.NaN)
msk2017_df = msk2017_df[["PATIENT_ID","SAMPLE_ID", "SMOKING_HISTORY", "STAGE_AT_DIAGNOSIS", "VITAL_STATUS","Treatment","SAMPLE_TYPE"]]
#using sample_id instead of patient_id because that is what matches with patient_id in the maf file
#msk2017 only has vital status rather than PFS or OS.
msk2017_df.columns = ["Patient ID","Sample ID","Smoker","Stage","Vital Status","Treatment","Metastatic"]
msk2017_df['is_LUAD'] = True

msk2018_df = pd.read_csv('output/unmerged_clinical/nsclc_msk_2018_clinical.txt')
msk2018_df['SMOKER'] = msk2018_df['SMOKER'].apply(lambda x: True if x == 'Ever' else False if x == 'Never' else np.NaN)
msk2018_df['Stage'] = np.NaN
#all patients treated with ICI
msk2018_df['Treatment'] = True
msk2018_df['Metastatic'] = np.NaN
msk2018_df['is_LUAD'] = msk2018_df['CANCER_TYPE_DETAILED'].apply(lambda x: True if x == 'Lung Adenocarcinoma' else False)
msk2018_df = msk2018_df[['PATIENT_ID','SAMPLE_ID','SMOKER','Stage','PFS_MONTHS','Treatment','Metastatic','is_LUAD']]
msk2018_df.columns = ["Patient ID","Sample ID","Smoker","Stage","Progression Free Survival (months)","Treatment","Metastatic","is_LUAD"]

tracer_df = pd.read_csv('output/unmerged_clinical/nsclc_tracerx_clinical.txt')
tracer_df['SMOKING_HISTORY'] = tracer_df['SMOKING_HISTORY'].apply(lambda x: True if x in ("Current Smoker", "Ex-Smoker", "Recent Ex-Smoker") else False if x == "Never Smoked" else np.NaN)
tracer_df['TUMOR_STAGE'] = tracer_df['TUMOR_STAGE'].map(stage_dict).fillna(tracer_df['TUMOR_STAGE'])
tracer_df['Treatment'] = tracer_df['SAMPLE_COLLECTION_TIMEPOINT'].apply(lambda x: True if x == 'Post-treatment' else False if x == 'Pre-treatment' else np.NaN)
#none of the tumors are metastatic but this just checks to make sure
tracer_df['Metastatic'] = tracer_df['SAMPLE_TYPE'].apply(lambda x: True if 'met' in x.lower() else False if x in ('Normal','Primary','Recurrence') else np.NaN)
tracer_df['is_LUAD'] = tracer_df['HISTOLOGY'].apply(lambda x: True if x == 'Invasive adenocarcinoma' else False if x in ('Adenosquamous carcinoma','Carcinosarcoma','Large cell carcinoma','Large Cell Neuroendocrine','Squamous cell carcinoma') else np.NaN)
tracer_df = tracer_df[['SAMPLE_ID','SMOKING_HISTORY','TUMOR_STAGE','RFS_MONTHS','Treatment','Metastatic','is_LUAD']]
#not sure if regression free survival = progression free survival
tracer_df.columns = ["Sample ID","Smoker","Stage","Progression Free Survival (months)","Treatment","Metastatic","is_LUAD"]
tracer_df.to_csv('output/unmerged_clinical/nsclc_tracerx_clinical.txt')

#removing repeated patients between msk 2017 and msk 2018
merged_temp = pd.merge(msk2017_df, msk2018_df, on = 'Patient ID', how = 'inner')
dup_patient_ids = merged_temp['Patient ID']
dup_sample_ids = msk2017_df[msk2017_df['Patient ID'].isin(dup_patient_ids)]['Sample ID']
msk2017_df = msk2017_df[~msk2017_df['Patient ID'].isin(dup_patient_ids)]

msk2017_df = msk2017_df.drop(columns='Patient ID')
msk2018_df = msk2018_df.drop(columns='Patient ID')

msk2017_df.to_csv("output/unmerged_clinical/lung_msk_2017_clinical.txt")
msk2018_df.to_csv('output/unmerged_clinical/nsclc_msk_2018_clinical.txt')


"""TSP not included because desired information is not in dataset and is not currently accessible"""

all_clinical_df = pd.concat([broad_df, tcga_df, oncosg_df, msk2015_df, msk2017_df, msk2018_df, tracer_df])

all_clinical_df.to_csv("output/luad_all_clinical.txt")