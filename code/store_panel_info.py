import pandas as pd

from locations import merged_maf_file_name
from locations import all_panel_genes_file_name
from locations import all_panel_samples_file_name

maf = pd.read_csv(merged_maf_file_name)

'''
Store samples per panel
'''

samples_per_panel = maf[~maf[['Sample ID', 'Panel']].duplicated()][['Sample ID', 'Panel']]
samples_per_panel = samples_per_panel[samples_per_panel['Panel'].notnull()]
samples_per_panel.to_csv(all_panel_samples_file_name, index = False)