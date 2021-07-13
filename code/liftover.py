import lift_coords
import pandas as pd
import os

g38hg38 = ['luad_tcga','luad_fm-ad']
g37hg38 = ['luad_oncosg_2020','luad_broad','luad_mskcc_2015','luad_tsp','lung_msk_2017','nsclc_pd1_msk_2018']#,'nsclc_tracerx_2017']
#tracerX doesn't work for some reason
hg38 = ['genie_9']

for i in g37hg38:
    df = pd.read_csv('data/' + i + '/data_mutations_extended.txt', sep="\t", comment="#")
    df2, failed = lift_coords.lift_over(df, 'grch37', 'grch38')#, keep_orig=True)
    #df3, failed = lift_coords.lift_over(df2, 'hg19', 'hg38')
    print(failed.shape)
    print(df2.head())
    #print(df3)