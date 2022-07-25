## Get relevant gammas to make plotting in R easier. In python this is
## done with import_results.py

import pandas as pd
import numpy as np
from cancer_epistasis import compute_gammas
from cancer_epistasis import compute_CI_gamma
from locations import full_mutation_rate_file_names
from locations import full_flux_mle_file_names
from locations import full_flux_ci_file_names

mutation_rate_dict = {
    key: pd.read_csv(full_mutation_rate_file_names[key]).set_index('gene').to_dict()['rate']
    for key in full_mutation_rate_file_names.keys()
}

fluxes_mles_dict = {
    key: np.load(full_flux_mle_file_names[key], allow_pickle=True).item()
    for key in full_flux_mle_file_names.keys()
}

fluxes_cis_dict = {
    key: np.load(full_flux_ci_file_names[key], allow_pickle=True).item()
    for key in full_flux_ci_file_names.keys()
}

'''
pan_data_fluxes_mles = np.load(location_output + 'pan_data_fluxes_mles.npy', allow_pickle = True).item()
smoking_fluxes_mles = np.load(location_output + 'smoking_fluxes_mles.npy', allow_pickle = True).item()
nonsmoking_fluxes_mles = np.load(location_output + 'nonsmoking_fluxes_mles.npy', allow_pickle = True).item()
'''

missing_mus = set()

pan_data_gammas_from_WT = dict()
pan_data_gamma_cis_from_WT = dict()
pan_data_gammas_from_110 = dict()
pan_data_gamma_cis_from_110 = dict()
for i in fluxes_mles_dict['pan_data'].keys():
    if i not in mutation_rate_dict['pan_data'].keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): mutation_rate_dict['pan_data'][genes3_[0]],
            (0, 1, 0): mutation_rate_dict['pan_data'][genes3_[1]],
            (0, 0, 1): mutation_rate_dict['pan_data'][genes3_[2]]}
    pan_data_gammas_from_WT[i] = compute_gammas(fluxes_mles_dict['pan_data'][i], mus3_)[((0, 0, 0), (0, 0, 1))]
    pan_data_gamma_cis_from_WT[i] = compute_CI_gamma(fluxes_cis_dict['pan_data'][i], mus3_)[((0, 0, 0), (0, 0, 1))]
    pan_data_gammas_from_110[i] = compute_gammas(fluxes_mles_dict['pan_data'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
    pan_data_gamma_cis_from_110[i] = compute_CI_gamma(fluxes_cis_dict['pan_data'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
#np.save(os.path.join(location_output, 'pan_data' + '_gammas.npy'), pan_data_gammas)

smoking_gammas_from_WT = dict()
smoking_gamma_cis_from_WT = dict()
smoking_gammas_from_110 = dict()
smoking_gamma_cis_from_110 = dict()
for i in fluxes_mles_dict['smoking_plus'].keys():
    if i not in mutation_rate_dict['smoking'].keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): mutation_rate_dict['smoking'][genes3_[0]],
            (0, 1, 0): mutation_rate_dict['smoking'][genes3_[1]],
            (0, 0, 1): mutation_rate_dict['smoking'][genes3_[2]]}
    smoking_gammas_from_WT[i] = compute_gammas(fluxes_mles_dict['smoking_plus'][i], mus3_)[((0, 0, 0), (0, 0, 1))]
    smoking_gamma_cis_from_WT[i] = compute_CI_gamma(fluxes_cis_dict['smoking_plus'][i], mus3_)[((0, 0, 0), (0, 0, 1))]
    smoking_gammas_from_110[i] = compute_gammas(fluxes_mles_dict['smoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
    smoking_gamma_cis_from_110[i] = compute_CI_gamma(fluxes_cis_dict['smoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
#np.save(os.path.join(location_output, str(samples_used + '_fluxes_cis.npy')), lambdas_cis)

nonsmoking_gammas_from_WT = dict()
nonsmoking_gamma_cis_from_WT = dict()
nonsmoking_gammas_from_110 = dict()
nonsmoking_gamma_cis_from_110 = dict()
for i in fluxes_mles_dict['nonsmoking_plus'].keys():
    if i in ['TTF1', 'GOPC','RAD17','PDGFRA','CCND1','RBP1','LMTK2','CCNE1']: continue
    if i not in mutation_rate_dict['nonsmoking'].keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): mutation_rate_dict['nonsmoking'][genes3_[0]],
            (0, 1, 0): mutation_rate_dict['nonsmoking'][genes3_[1]],
            (0, 0, 1): mutation_rate_dict['nonsmoking'][genes3_[2]]}
    nonsmoking_gammas_from_WT[i] = compute_gammas(fluxes_mles_dict['nonsmoking_plus'][i], mus3_)[((0, 0, 0), (0, 0, 1))]
    nonsmoking_gamma_cis_from_WT[i] = compute_CI_gamma(fluxes_cis_dict['nonsmoking_plus'][i], mus3_)[((0, 0, 0), (0, 0, 1))]
    nonsmoking_gammas_from_110[i] = compute_gammas(fluxes_mles_dict['nonsmoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
    nonsmoking_gamma_cis_from_110[i] = compute_CI_gamma(fluxes_cis_dict['nonsmoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]

with open('../data/missing_mus.txt','w') as missing_mus_file:
    for i in missing_mus:
        missing_mus_file.write(i + '\n')

pan_data_gammas_from_110 = pd.DataFrame.from_dict(pan_data_gammas_from_110, orient = 'index', columns = ['pd_from_110_gamma'])
smoking_gammas_from_110 = pd.DataFrame.from_dict(smoking_gammas_from_110, orient = 'index', columns = ['s+_from_110_gamma'])
nonsmoking_gammas_from_110 = pd.DataFrame.from_dict(nonsmoking_gammas_from_110, orient = 'index', columns = ['ns+_from_110_gamma'])

pan_data_gamma_cis_from_110 = pd.DataFrame(pan_data_gamma_cis_from_110).transpose().rename(columns={0:'pd_from_110_low',1:'pd_from_110_high'})
smoking_gamma_cis_from_110 = pd.DataFrame(smoking_gamma_cis_from_110).transpose().rename(columns={0:'s+_from_110_low',1:'s+_from_110_high'})
nonsmoking_gamma_cis_from_110 = pd.DataFrame(nonsmoking_gamma_cis_from_110).transpose().rename(columns={0:'ns+_from_110_low',1:'ns+_from_110_high'})

pan_data_gammas_from_WT = pd.DataFrame.from_dict(pan_data_gammas_from_WT, orient = 'index', columns = ['pd_from_WT_gamma'])
smoking_gammas_from_WT = pd.DataFrame.from_dict(smoking_gammas_from_WT, orient = 'index', columns = ['s+_from_WT_gamma'])
nonsmoking_gammas_from_WT = pd.DataFrame.from_dict(nonsmoking_gammas_from_WT, orient = 'index', columns = ['ns+_from_WT_gamma'])

pan_data_gamma_cis_from_WT = pd.DataFrame(pan_data_gamma_cis_from_WT).transpose().rename(columns={0:'pd_from_WT_low',1:'pd_from_WT_high'})
smoking_gamma_cis_from_WT = pd.DataFrame(smoking_gamma_cis_from_WT).transpose().rename(columns={0:'s+_from_WT_low',1:'s+_from_WT_high'})
nonsmoking_gamma_cis_from_WT = pd.DataFrame(nonsmoking_gamma_cis_from_WT).transpose().rename(columns={0:'ns+_from_WT_low',1:'ns+_from_WT_high'})


comparison_genes_by_gamma_from_110 = pan_data_gammas_from_110.join(smoking_gammas_from_110.join(nonsmoking_gammas_from_110, how = 'outer'), how = 'outer')
comparison_gamma_cis_from_110 = pan_data_gamma_cis_from_110.join(smoking_gamma_cis_from_110.join(nonsmoking_gamma_cis_from_110, how = 'outer'), how = 'outer')
comparison_genes_by_gamma_from_110 = comparison_genes_by_gamma_from_110.join(comparison_gamma_cis_from_110, how='outer')

comparison_genes_by_gamma_from_WT = pan_data_gammas_from_WT.join(smoking_gammas_from_WT.join(nonsmoking_gammas_from_WT, how = 'outer'), how = 'outer')
comparison_gamma_cis_from_WT = pan_data_gamma_cis_from_WT.join(smoking_gamma_cis_from_WT.join(nonsmoking_gamma_cis_from_WT, how = 'outer'), how = 'outer')
comparison_genes_by_gamma_from_WT = comparison_genes_by_gamma_from_WT.join(comparison_gamma_cis_from_WT, how='outer')

all_comparison_genes_by_gamma = comparison_genes_by_gamma_from_110.join(comparison_genes_by_gamma_from_WT, how = 'outer')

all_comparison_genes_by_gamma.to_csv('../output/genes_by_gamma.csv')
