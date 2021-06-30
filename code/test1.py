import pandas as pd
from importing_clinical_data import dup_patient_ids

df1 = pd.read_csv('data/clinical_data/lung_msk_2017/data_clinical_patient.txt', sep = '\t', comment='#')
df2 = pd.read_csv('data/clinical_data/lung_msk_2017/data_clinical_sample.txt', sep = '\t', comment='#')
df = pd.merge(df1, df2, on="PATIENT_ID")
#merged_df = merged_df[merged_df['HISTOLOGY'] == 'Invasive adenocarcinoma']

df = df[~df['PATIENT_ID'].isin(dup_patient_ids)]
df = df[df.duplicated(subset=['PATIENT_ID'], keep=False)]
test_df = df.sort_values(by = ['LINES_OF_TX_PRIOR_IMPACT','SAMPLE_COVERAGE'], ascending=False)
check_df = test_df[test_df.duplicated(subset=['PATIENT_ID'], keep='first')]['SAMPLE_ID']
df = df[['PATIENT_ID','SAMPLE_ID','SAMPLE_TYPE','METASTATIC_SITE','PRIMARY_SITE','SO_COMMENTS','SAMPLE_COVERAGE','TUMOR_PURITY','SOMATIC_STATUS','GENE_PANEL','LINES_OF_TX_PRIOR_IMPACT','DRIVER_MUTATIONS','SAMPLE_PRE_LUNG_THERAPY']]
'''
for i in range(df.shape[0])[::2]:
    comparison_df = df.iloc[i].compare(df.iloc[i+1], align_axis = 0, keep_shape = True)
    comparison_df.to_csv('what.txt')
'''
df.to_csv('test_msk.txt')
check_df.to_csv('check_test.txt')