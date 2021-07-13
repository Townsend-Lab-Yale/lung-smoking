import lift_coords
## The package lift_coords requires to have the LiftOver command-line
## program. Obtain it from https://genome-store.ucsc.edu/
import pandas as pd
import os


if '__file__' not in globals():
    __file__ = '.'

location_data = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "data/"))
"""Location of directory that contains data for the model."""


g38hg38 = ['luad_tcga','luad_fm-ad']
g37hg38 = ['luad_oncosg_2020','luad_broad','luad_mskcc_2015','luad_tsp','lung_msk_2017','nsclc_pd1_msk_2018']#,'nsclc_tracerx_2017']
#tracerX doesn't work for some reason
hg38 = ['genie_9']

for db in g37hg38:
    df = pd.read_csv(os.path.join(location_data, db, "data_mutations_extended.txt"),
                     sep="\t",
                     comment="#")
    df2, failed = lift_coords.lift_over(df, 'grch37', 'grch38')#, keep_orig=True)
    #df3, failed = lift_coords.lift_over(df2, 'hg19', 'hg38')
    print(failed.shape)
    print(df2.head())
    #print(df3)
