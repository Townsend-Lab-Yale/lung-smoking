import pandas as pd
import numpy as np
from cancer_epistasis import compute_gammas

pan_data_mutation_rates_file = pd.read_csv('../data/pan-data_mutation_rates.txt')
pan_data_mutation_rates = dict(zip(pan_data_mutation_rates_file['gene'], pan_data_mutation_rates_file['rate']))

smoking_mutation_rates_file = pd.read_csv('../data/smoking_mutation_rates.txt')
smoking_mutation_rates = dict(zip(smoking_mutation_rates_file['gene'], smoking_mutation_rates_file['rate']))

nonsmoking_mutation_rates_file = pd.read_csv('../data/nonsmoking_mutation_rates.txt')
nonsmoking_mutation_rates = dict(zip(nonsmoking_mutation_rates_file['gene'], nonsmoking_mutation_rates_file['rate']))

pan_data_fluxes_mles = np.load('../output/pan_data_fluxes_mles.npy', allow_pickle = True).item()
smoking_fluxes_mles = np.load('../output/smoking_fluxes_mles.npy', allow_pickle = True).item()
nonsmoking_fluxes_mles = np.load('../output/nonsmoking_fluxes_mles.npy', allow_pickle = True).item()

'''
pan_data_gammas = compute_gammas(pan_data_fluxes_mles, pan_data_mutation_rates)
pan_data_gammas = pd.DataFrame.from_dict(pan_data_gammas, orient = 'index', columns = ['gamma'])
pan_data_gammas['gene'] = pan_data_gammas.index

smoking_gammas = compute_gammas(smoking_fluxes_mles, smoking_mutation_rates)
smoking_gammas = pd.DataFrame.from_dict(smoking_gammas, orient = 'index', columns = ['gamma'])
smoking_gammas['gene'] = smoking_gammas.index

nonsmoking_gammas = compute_gammas(nonsmoking_fluxes_mles, nonsmoking_mutation_rates)
nonsmoking_gammas = pd.DataFrame.from_dict(nonsmoking_gammas, orient = 'index', columns = ['gamma'])
nonsmoking_gammas['gene'] = nonsmoking_gammas.index

comparison_genes_by_gamma = pd.merge(pan_data_gammas, pd.merge(smoking_gammas, nonsmoking_gammas, how = 'outer',on = 'gene'), how = 'outer', on = 'gene')
comparison_genes_by_gamma = comparison_genes_by_gamma.rename(columns = {'gamma':'pan_data_gamma', 'flux_x':'smoking_gamma', 'flux_y': 'nonsmoking_gamma'}).set_index('gene')

print(comparison_genes_by_gamma)
'''

pan_data_gammas = dict()
missing_mus = set()
for i in pan_data_fluxes_mles.keys():
    if i not in pan_data_mutation_rates.keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): pan_data_mutation_rates[genes3_[0]],
            (0, 1, 0): pan_data_mutation_rates[genes3_[1]],
            (0, 0, 1): pan_data_mutation_rates[genes3_[2]]}
    pan_data_gammas[i] = compute_gammas(pan_data_fluxes_mles[i], mus3_)[((1, 1, 0), (1, 1, 1))]

smoking_gammas = dict()
for i in smoking_fluxes_mles.keys():
    if i not in smoking_mutation_rates.keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): smoking_mutation_rates[genes3_[0]],
            (0, 1, 0): smoking_mutation_rates[genes3_[1]],
            (0, 0, 1): smoking_mutation_rates[genes3_[2]]}
    smoking_gammas[i] = compute_gammas(smoking_fluxes_mles[i], mus3_)[((1, 1, 0), (1, 1, 1))]

nonsmoking_gammas = dict()
for i in nonsmoking_fluxes_mles.keys():
    if i not in nonsmoking_mutation_rates.keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): nonsmoking_mutation_rates[genes3_[0]],
            (0, 1, 0): nonsmoking_mutation_rates[genes3_[1]],
            (0, 0, 1): nonsmoking_mutation_rates[genes3_[2]]}
    nonsmoking_gammas[i] = compute_gammas(nonsmoking_fluxes_mles[i], mus3_)[((1, 1, 0), (1, 1, 1))]

with open('../data/missing_mus.txt','w') as missing_mus_file:
    for i in missing_mus:
        missing_mus_file.write(i + '\n')

pan_data_gammas = pd.DataFrame.from_dict(pan_data_gammas, orient = 'index', columns = ['pan_data_gamma'])
pan_data_gammas = pan_data_gammas.reset_index().rename(columns={'index':'gene'})

smoking_gammas = pd.DataFrame.from_dict(smoking_gammas, orient = 'index', columns = ['smoking_gamma'])
smoking_gammas = smoking_gammas.reset_index().rename(columns={'index':'gene'})

nonsmoking_gammas = pd.DataFrame.from_dict(nonsmoking_gammas, orient = 'index', columns = ['nonsmoking_gamma'])
nonsmoking_gammas = nonsmoking_gammas.reset_index().rename(columns={'index':'gene'})

comparison_genes_by_gamma = pd.merge(pan_data_gammas, pd.merge(smoking_gammas, nonsmoking_gammas, how = 'outer',on = 'gene'), how = 'outer', on = 'gene')

comparison_genes_by_gamma.to_csv('../data/genes_by_gamma.csv')


