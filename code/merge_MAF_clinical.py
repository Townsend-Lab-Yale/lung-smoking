import pandas as pd

maf = pd.read_csv("output/merged_luad_maf.txt")
clinical = pd.read_csv("output/luad_all_clinical.txt")

#outer merge to keep MAF data without clinical data
maf_clinical = pd.merge(maf, clinical, on="Sample ID", how = 'outer')
#removing artifacts of merge
maf_clinical = maf_clinical.drop(columns=['Unnamed: 0_x','Unnamed: 0_y'])
#removes non-LUAD NSCLC data
maf_clinical = maf_clinical[maf_clinical['is_LUAD'] == True]
#removes mutations with unknown chromosome locations
maf_clinical = maf_clinical[~(maf_clinical['Chromosome'].isin(['GL000230.1', 'hs37d5', 'GL000211.1','MT', 'GL000192.1', 'GL000214.1', 'GL000241.1', 'GL000220.1', 'GL000212.1','GL000205.1', 'GL000195.1', 'GL000218.1', 'GL000216.1', 'GL000226.1','GL000224.1', 'GL000231.1', 'GL000221.1', 'GL000234.1', 'GL000219.1','GL000191.1', 'GL000229.1', 'GL000238.1']))]

maf_clinical.to_csv("output/luad_maf_clinical.txt")

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