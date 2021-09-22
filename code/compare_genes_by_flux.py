import numpy as np
from locations import fluxes_mles_file_name
from locations import fluxes_cis_file_name
import pandas as pd

pan_data_fluxes_mles = np.load('../output/pan_data_fluxes_mles.npy', allow_pickle = True).item()
ranked_genes_by_fluxes = [i[0] for i in sorted(pan_data_fluxes_mles.items(), key = lambda item: item[1][((1, 1, 0), (1, 1, 1))])]
ranked_genes_by_fluxes = dict(zip(ranked_genes_by_fluxes, [pan_data_fluxes_mles[gene][((1, 1, 0), (1, 1, 1))] for gene in ranked_genes_by_fluxes]))

#print(ranked_genes_by_fluxes)

smoking_fluxes_mles = np.load('../output/smoking_fluxes_mles.npy', allow_pickle = True).item()
smoking_genes_by_fluxes = [i[0] for i in sorted(smoking_fluxes_mles.items(), key = lambda item: item[1][((1, 1, 0), (1, 1, 1))])]
smoking_genes_by_fluxes = dict(zip(smoking_genes_by_fluxes, [smoking_fluxes_mles[gene][((1, 1, 0), (1, 1, 1))] for gene in smoking_genes_by_fluxes]))

#print(smoking_genes_by_fluxes)

nonsmoking_fluxes_mles = np.load('../output/nonsmoking_fluxes_mles.npy', allow_pickle = True).item()
nonsmoking_genes_by_fluxes = [i[0] for i in sorted(nonsmoking_fluxes_mles.items(), key = lambda item: item[1][((1, 1, 0), (1, 1, 1))])]
nonsmoking_genes_by_fluxes = dict(zip(nonsmoking_genes_by_fluxes, [nonsmoking_fluxes_mles[gene][((1, 1, 0), (1, 1, 1))] for gene in nonsmoking_genes_by_fluxes]))

#print(nonsmoking_genes_by_fluxes)

ranked_genes_by_fluxes = pd.DataFrame.from_dict(ranked_genes_by_fluxes, orient = 'index', columns = ['pan_data_flux'])
ranked_genes_by_fluxes['gene'] = ranked_genes_by_fluxes.index
smoking_genes_by_fluxes = pd.DataFrame.from_dict(smoking_genes_by_fluxes, orient = 'index', columns = ['smoking_flux'])
smoking_genes_by_fluxes['gene'] = smoking_genes_by_fluxes.index
nonsmoking_genes_by_fluxes = pd.DataFrame.from_dict(nonsmoking_genes_by_fluxes, orient = 'index', columns = ['nonsmoking_flux'])
nonsmoking_genes_by_fluxes['gene'] = nonsmoking_genes_by_fluxes.index

#print(pd.merge(smoking_genes_by_fluxes, nonsmoking_genes_by_fluxes, how = 'outer', on = 'gene'))
comparison_genes_by_fluxes = pd.merge(ranked_genes_by_fluxes, pd.merge(smoking_genes_by_fluxes, nonsmoking_genes_by_fluxes, how = 'outer',on = 'gene'), how = 'outer', on = 'gene')
comparison_genes_by_fluxes = comparison_genes_by_fluxes.set_index('gene')

comparison_genes_by_fluxes.to_csv('../data/genes_by_fluxes.csv')

#print(comparison_genes_by_fluxes)