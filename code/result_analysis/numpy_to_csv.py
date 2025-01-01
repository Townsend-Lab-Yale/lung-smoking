import os
import sys
import pandas as pd

sys.path.insert(0, '../')
from theory import build_S_with_tuples
from load_results import load_results

def convert_all_numpy_to_csv(extension=None, mu_method="variant"):
    if extension is None: 
        subset_extension = "subset"
        all_extension = "all"
    else: 
        subset_extension = os.path.join(extension,"subset")
        all_extension = os.path.join(extension,"all")

    fluxes_mles = load_results('fluxes','mles',subset_extension)
    fluxes_cis = load_results('fluxes','cis',subset_extension)

    selection_mles = load_results('selections','mles',subset_extension)
    selection_cis = load_results('selections','cis',subset_extension)

    if mu_method == 'cesR': 
        location_output = "cesR_results"
        mu_dict = load_results('mutations','cesR')
    elif mu_method == 'variant': 
        location_output = "variant_results"
        mu_dict = load_results('mutations','variant')
    os.mkdir(location_output)

    all_samples = load_results('samples', extension=subset_extension)

    # Write mutation rates to CSV
    mu_df = pd.json_normalize({mu_method:mu_dict},sep='|').T.reset_index()
    mu_df.columns = ['index','mu']
    mu_df[['method','key','gene']] = mu_df['index'].str.split('|', expand=True)
    mu_df = mu_df[['method','key','gene','mu']]
    mu_df.to_csv(os.path.join(location_output, 'mutation_rates.csv'), index=False)

    # Write fluxes to CSV
    flux_df = pd.json_normalize(fluxes_mles,sep='|').T.reset_index()
    flux_df.columns = ['index','flux_mle']

    tmp = pd.json_normalize(fluxes_cis,sep='|').T.reset_index()
    tmp.columns = ['index','flux_cis']
    tmp['flux_ci_low'] = tmp['flux_cis'].apply(lambda x: x[0])
    tmp['flux_ci_high'] = tmp['flux_cis'].apply(lambda x: x[1])
    tmp = tmp.drop('flux_cis',axis=1)

    flux_df = pd.merge(flux_df,tmp,on='index')
    flux_df[['key','gene_set','mutation']] = flux_df['index'].str.split('|', expand=True)
    flux_df = flux_df[['key','gene_set','mutation','flux_mle','flux_ci_low','flux_ci_high']]
    flux_df[['first_gene','second_gene','third_gene']] = flux_df['gene_set'].str.strip(',()').str.replace("'","").str.split(', ',expand=True)
    flux_df = flux_df.drop('gene_set', axis=1)

    M1_flux_df = flux_df[flux_df['second_gene'].isnull()]
    M1_flux_df = M1_flux_df.drop(['second_gene','third_gene'],axis=1).rename(columns={'first_gene':'gene'})
    M1_flux_df = M1_flux_df[['key','gene','flux_mle','flux_ci_low','flux_ci_high']]
    M1_flux_df.to_csv(os.path.join(location_output, 'M1_gene_fluxes.csv'), index=False)

    M2_flux_df = flux_df[(flux_df['second_gene'].notnull()) & (flux_df['third_gene'].isnull())]
    M2_flux_df = M2_flux_df.drop('third_gene',axis=1)
    M2_flux_df = M2_flux_df[['key','first_gene','second_gene','mutation','flux_mle','flux_ci_low','flux_ci_high']]
    M2_flux_df.to_csv(os.path.join(location_output, 'M2_gene_fluxes.csv'), index=False)

    M3_flux_df = flux_df[flux_df['third_gene'].notnull()]
    M3_flux_df = M3_flux_df[['key','first_gene','second_gene','third_gene','mutation','flux_mle','flux_ci_low','flux_ci_high']]
    M3_flux_df.to_csv(os.path.join(location_output, 'M3_gene_fluxes.csv'), index=False)

    # Write gammas to CSV
    gamma_df = pd.json_normalize(selection_mles,sep='|').T.reset_index()
    gamma_df.columns = ['index','gamma_mle']

    tmp = pd.json_normalize(selection_cis,sep='|').T.reset_index()
    tmp.columns = ['index','gamma_cis']
    tmp['gamma_ci_low'] = tmp['gamma_cis'].apply(lambda x: x[0])
    tmp['gamma_ci_high'] = tmp['gamma_cis'].apply(lambda x: x[1])
    tmp = tmp.drop('gamma_cis',axis=1)

    gamma_df = pd.merge(gamma_df,tmp,on='index')

    gamma_df[['key','gene_set','mutation']] = gamma_df['index'].str.split('|', expand=True)
    gamma_df = gamma_df[['key','gene_set','mutation','gamma_mle','gamma_ci_low','gamma_ci_high']]

    gamma_df[['first_gene','second_gene','third_gene']] = gamma_df['gene_set'].str.strip(',()').str.replace("'","").str.split(', ',expand=True)
    gamma_df = gamma_df.drop('gene_set', axis=1)
    gamma_df['method'] = mu_method

    M1_gamma_df = gamma_df[gamma_df['second_gene'].isnull()]
    M1_gamma_df = M1_gamma_df.drop(['second_gene','third_gene'],axis=1).rename(columns={'first_gene':'gene'})
    M1_gamma_df = M1_gamma_df[['method','key','gene','gamma_mle','gamma_ci_low','gamma_ci_high']]
    M1_gamma_df.to_csv(os.path.join(location_output, 'M1_gene_gammas.csv'), index=False)

    M2_gamma_df = gamma_df[(gamma_df['second_gene'].notnull()) & (gamma_df['third_gene'].isnull())]
    M2_gamma_df = M2_gamma_df.drop('third_gene',axis=1)
    M2_gamma_df = M2_gamma_df[['method','key','first_gene','second_gene','mutation','gamma_mle','gamma_ci_low','gamma_ci_high']]
    M2_gamma_df.to_csv(os.path.join(location_output, 'M2_gene_gammas.csv'), index=False)

    M3_gamma_df = gamma_df[gamma_df['third_gene'].notnull()]
    M3_gamma_df = M3_gamma_df[['method','key','first_gene','second_gene','third_gene','mutation','gamma_mle','gamma_ci_low','gamma_ci_high']]
    M3_gamma_df.to_csv(os.path.join(location_output, 'M3_gene_gammas.csv'), index=False)

    # Write samples per combination to CSV
    sample_df = pd.json_normalize(all_samples,sep='|').T.reset_index().rename(columns={0:'counts'})
    sample_df[['key','gene_set']] = sample_df['index'].str.split('|', expand=True)
    sample_df = sample_df.drop('index',axis=1)
    sample_df['gene_set'] = sample_df['gene_set'].str.strip(',()')

    sample_dict = dict()
    for i in range(1,4): 
        sample_dict[i] = sample_df[sample_df['counts'].apply(len)==2**i]
        sample_dict[i][[str(x) for x in build_S_with_tuples(i)]] = sample_dict[i]['counts'].apply(pd.Series)
        
        if i == 1: 
            sample_dict[i]['gene'] = sample_dict[i]['gene_set'].str.replace("'","")
            sample_dict[i] = sample_dict[i][['key','gene'] + [str(x) for x in build_S_with_tuples(i)]]
        elif i == 2: 
            sample_dict[i][['first_gene','second_gene']] = sample_dict[i]['gene_set'].str.replace("'","").str.split(', ',expand=True)
            sample_dict[i] = sample_dict[i][['key','first_gene','second_gene'] + [str(x) for x in build_S_with_tuples(i)]]
        elif i == 3:
            sample_dict[i][['first_gene','second_gene','third_gene']] = sample_dict[i]['gene_set'].str.replace("'","").str.split(', ',expand=True)
            sample_dict[i] = sample_dict[i][['key','first_gene','second_gene','third_gene'] + [str(x) for x in build_S_with_tuples(i)]]

        sample_dict[i].to_csv(os.path.join(location_output, f'M{i}_samples_per_combination.csv'),index=False)

    # Write M=1 for 1200 genes results to CSV
    fluxes_mles = load_results('fluxes','mles',all_extension)
    fluxes_cis = load_results('fluxes','cis',all_extension)

    selection_mles = load_results('selections','mles',all_extension)
    selection_cis = load_results('selections','cis',all_extension)

    all_samples = load_results('samples',extension=all_extension)

    flux_df = pd.json_normalize(fluxes_mles,sep='|').T.reset_index()
    flux_df.columns = ['index','flux_mle']

    tmp = pd.json_normalize(fluxes_cis,sep='|').T.reset_index()
    tmp.columns = ['index','flux_cis']
    tmp['flux_ci_low'] = tmp['flux_cis'].apply(lambda x: x[0])
    tmp['flux_ci_high'] = tmp['flux_cis'].apply(lambda x: x[1])
    tmp = tmp.drop('flux_cis',axis=1)

    flux_df = pd.merge(flux_df,tmp,on='index')
    flux_df[['key','gene','mutation']] = flux_df['index'].str.split('|', expand=True)
    flux_df['gene'] = flux_df['gene'].str.strip(",()'")
    flux_df = flux_df[['key','gene','flux_mle','flux_ci_low','flux_ci_high']]
    flux_df.to_csv(os.path.join(location_output, 'M1_all_gene_fluxes.csv'), index=False)

    gamma_df = pd.json_normalize(selection_mles,sep='|').T.reset_index()
    gamma_df.columns = ['index','gamma_mle']

    tmp = pd.json_normalize(selection_cis,sep='|').T.reset_index()
    tmp.columns = ['index','gamma_cis']
    tmp['gamma_ci_low'] = tmp['gamma_cis'].apply(lambda x: x[0])
    tmp['gamma_ci_high'] = tmp['gamma_cis'].apply(lambda x: x[1])
    tmp = tmp.drop('gamma_cis',axis=1)

    gamma_df = pd.merge(gamma_df,tmp,on='index')
    gamma_df[['key','gene','mutation']] = gamma_df['index'].str.split('|', expand=True)
    gamma_df['gene'] = gamma_df['gene'].str.strip(",()'")
    gamma_df['method'] = mu_method
    gamma_df = gamma_df[['method','key','gene','gamma_mle','gamma_ci_low','gamma_ci_high']]
    gamma_df.to_csv(os.path.join(location_output, 'M1_all_gene_gammas.csv'), index=False)

    sample_df = pd.json_normalize(all_samples,sep='|').T.reset_index().rename(columns={0:'counts'})
    sample_df[['key','gene_set']] = sample_df['index'].str.split('|', expand=True)
    sample_df = sample_df.drop('index',axis=1)
    sample_df['gene_set'] = sample_df['gene_set'].str.strip(",()'")
    sample_df[[str(x) for x in build_S_with_tuples(1)]] = sample_df['counts'].apply(pd.Series)
    sample_df['gene'] = sample_df['gene_set'].str.replace("'","")
    sample_df = sample_df[['key','gene'] + [str(x) for x in build_S_with_tuples(1)]]

    sample_df.to_csv(os.path.join(location_output, 'M1_all_samples_per_combination.csv'),index=False)