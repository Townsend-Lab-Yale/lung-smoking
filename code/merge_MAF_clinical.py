import pandas as pd

from locations import merged_maf_file_name
from locations import merged_clinical_file_name
from locations import merged_maf_clinical_file_name

maf = pd.read_csv(merged_maf_file_name)
clinical = pd.read_csv(merged_clinical_file_name)

#outer merge to keep MAF data without clinical data
maf_clinical = pd.merge(maf, clinical, on="Sample ID", how = 'outer')
#removing artifacts of merge
maf_clinical = maf_clinical.drop(columns=['Unnamed: 0_x','Unnamed: 0_y'])

maf_clinical.to_csv(merged_maf_clinical_file_name)

'''Below are checks to ensure that the columns contain expected values.'''
'''
print(maf_clinical['Chromosome'].unique())
#not sure how to check these three columns
#print (maf_clinical[pd.to_numeric(maf_clinical['Start_Position'], errors='coerce').isnull()])
#print (maf_clinical[pd.to_numeric(maf_clinical['Progression Free Survival (months)'], errors='coerce').isnull()])
#print (maf_clinical[pd.to_numeric(maf_clinical['Overall Survival (months)'], errors='coerce').isnull()])
print(maf_clinical['Variant_Classification'].unique())
print(maf_clinical['Smoker'].unique())
print(maf_clinical['Stage'].unique())
print(maf_clinical['Treatment'].unique())
print(maf_clinical['Metastatic'].unique())
print(maf_clinical['Vital Status'].unique())
'''
