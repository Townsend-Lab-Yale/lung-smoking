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

pan_data_gammas = dict()
pan_data_gamma_cis = dict()
for i in fluxes_mles_dict['pan_data'].keys():
    if i not in mutation_rate_dict['pan_data'].keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): mutation_rate_dict['pan_data'][genes3_[0]],
            (0, 1, 0): mutation_rate_dict['pan_data'][genes3_[1]],
            (0, 0, 1): mutation_rate_dict['pan_data'][genes3_[2]]}
    pan_data_gammas[i] = compute_gammas(fluxes_mles_dict['pan_data'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
    pan_data_gamma_cis[i] = compute_CI_gamma(fluxes_cis_dict['pan_data'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
#np.save(os.path.join(location_output, 'pan_data' + '_gammas.npy'), pan_data_gammas)

smoking_gammas = dict()
smoking_gamma_cis = dict()
for i in fluxes_mles_dict['smoking_plus'].keys():
    if i not in mutation_rate_dict['smoking'].keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): mutation_rate_dict['smoking'][genes3_[0]],
            (0, 1, 0): mutation_rate_dict['smoking'][genes3_[1]],
            (0, 0, 1): mutation_rate_dict['smoking'][genes3_[2]]}
    smoking_gammas[i] = compute_gammas(fluxes_mles_dict['smoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
    smoking_gamma_cis[i] = compute_CI_gamma(fluxes_cis_dict['smoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
#np.save(os.path.join(location_output, str(samples_used + '_fluxes_cis.npy')), lambdas_cis)


nonsmoking_gammas = dict()
nonsmoking_gamma_cis = dict()
for i in fluxes_mles_dict['nonsmoking_plus'].keys():
    if i not in mutation_rate_dict['nonsmoking'].keys():
        missing_mus.add(i)
        continue
    genes3_ = ['TP53', 'KRAS', i]
    mus3_ = {(1, 0, 0): mutation_rate_dict['nonsmoking'][genes3_[0]],
            (0, 1, 0): mutation_rate_dict['nonsmoking'][genes3_[1]],
            (0, 0, 1): mutation_rate_dict['nonsmoking'][genes3_[2]]}
    nonsmoking_gammas[i] = compute_gammas(fluxes_mles_dict['nonsmoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]
    nonsmoking_gamma_cis[i] = compute_CI_gamma(fluxes_cis_dict['nonsmoking_plus'][i], mus3_)[((1, 1, 0), (1, 1, 1))]

with open('../data/missing_mus.txt','w') as missing_mus_file:
    for i in missing_mus:
        missing_mus_file.write(i + '\n')

pan_data_gammas = pd.DataFrame.from_dict(pan_data_gammas, orient = 'index', columns = ['pd_gamma'])
smoking_gammas = pd.DataFrame.from_dict(smoking_gammas, orient = 'index', columns = ['s_gamma'])
nonsmoking_gammas = pd.DataFrame.from_dict(nonsmoking_gammas, orient = 'index', columns = ['ns_gamma'])

pan_data_gamma_cis = pd.DataFrame(pan_data_gamma_cis).transpose().rename(columns={0:'pd_low',1:'pd_high'})
smoking_gamma_cis = pd.DataFrame(smoking_gamma_cis).transpose().rename(columns={0:'s_low',1:'s_high'})
nonsmoking_gamma_cis = pd.DataFrame(nonsmoking_gamma_cis).transpose().rename(columns={0:'ns_low',1:'ns_high'})

comparison_genes_by_gamma = pan_data_gammas.join(smoking_gammas.join(nonsmoking_gammas, how = 'outer'), how = 'outer')
comparison_gamma_cis = pan_data_gamma_cis.join(smoking_gamma_cis.join(nonsmoking_gamma_cis, how = 'outer'), how = 'outer')
comparison_genes_by_gamma = comparison_genes_by_gamma.join(comparison_gamma_cis, how='outer')

comparison_genes_by_gamma.to_csv('../output/genes_by_gamma.csv')