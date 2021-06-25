from os import defpath
import pandas as pd

df = pd.read_csv('output/luad_maf_clinical.txt')
print(df.shape)
print(df[~(df['Chromosome'].isin(['GL000230.1', 'hs37d5', 'GL000211.1','MT', 'GL000192.1', 'GL000214.1', 'GL000241.1', 'GL000220.1', 'GL000212.1','GL000205.1', 'GL000195.1', 'GL000218.1', 'GL000216.1', 'GL000226.1','GL000224.1', 'GL000231.1', 'GL000221.1', 'GL000234.1', 'GL000219.1','GL000191.1', 'GL000229.1', 'GL000238.1']))])
no_chr = df[df['Chromosome'].isin(['GL000230.1', 'hs37d5', 'GL000211.1','MT', 'GL000192.1', 'GL000214.1', 'GL000241.1', 'GL000220.1', 'GL000212.1','GL000205.1', 'GL000195.1', 'GL000218.1', 'GL000216.1', 'GL000226.1','GL000224.1', 'GL000231.1', 'GL000221.1', 'GL000234.1', 'GL000219.1','GL000191.1', 'GL000229.1', 'GL000238.1'])]['Sample ID']
no_chr.to_csv('no_chr.txt')
