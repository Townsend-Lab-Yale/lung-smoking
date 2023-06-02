import pandas as pd

from locations import merged_maf_file_name
from locations import all_panel_samples_file_name

'''Creates 2-column table with sample ID and the panel used for TGS of that sample'''

maf = pd.read_csv(merged_maf_file_name)

samples_per_panel = maf[~maf[['Sample ID', 'Panel']].duplicated()][['Sample ID', 'Panel']]
samples_per_panel = samples_per_panel[samples_per_panel['Panel'].notnull()]
samples_per_panel.to_csv(all_panel_samples_file_name, index = False)