import pandas as pd

from locations import merged_maf_file_name
from locations import merged_clinical_file_name
from locations import merged_maf_clinical_file_name

maf = pd.read_csv(merged_maf_file_name, index_col=0)
clinical = pd.read_csv(merged_clinical_file_name, index_col=0)

#outer merge to keep MAF data without clinical data
maf_clinical = pd.merge(maf, clinical, on="Sample ID", how = 'outer')

maf_clinical.to_csv(merged_maf_clinical_file_name)
