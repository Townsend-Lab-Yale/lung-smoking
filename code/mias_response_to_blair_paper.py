## Mia's Blair response project

import pandas as pd
import numpy as np

import cancer_epistasis
import importlib
importlib.reload(cancer_epistasis)

from cancer_epistasis import estimate_lambdas

from cancer_epistasis import p_value_same_lambda_xy
from cancer_epistasis import p_value_same_gamma_xy
from cancer_epistasis import p_value_gamma_xy_equal_1

## Import mutation rates

mutation_rates = pd.read_csv("../data/mias_blair_response/gene-level_neutral_mutation_rates.csv",
                             index_col=["gene", "data_type"])

oncovariants_mutation_rates = np.load(
    "../data/mias_blair_response/mean_oncogene_variant_mutation_rates.npy",
    allow_pickle=True).item()

for key in oncovariants_mutation_rates.keys():
    oncovariants_mutation_rates[key].index = [
        x[0] for x in oncovariants_mutation_rates[key].index.str.split("_ENSP")]

dataset_no_epi = pd.read_csv("../data/mias_blair_response/target_effect_samples.csv",
                             index_col=["variant_name", "group"])

oncovariants = ['KRAS_G12D', 'KRAS_G12C', 'BRAF_V600E', 'EGFR_L858R']

oncogenes = ["BRAF", "EGFR", "KRAS"]


## No epistasis

print("Computing estimates under H0 of same selection (without "
      "epistasis) on smokers vs non-smokers and respective p-values...")

for oncovariant in oncovariants:
    p_value = p_value_same_gamma_xy(
        [dataset_no_epi["included_without_variant"][oncovariant, 'Smoker'],
         dataset_no_epi["included_with_variant"][oncovariant, 'Smoker']],
        [dataset_no_epi["included_without_variant"][oncovariant, 'Never-smoker'],
         dataset_no_epi["included_with_variant"][oncovariant, 'Never-smoker']],
        ((0,), (1,)),
        oncovariants_mutation_rates['smoking'][oncovariant],
        oncovariants_mutation_rates['nonsmoking'][oncovariant],
        verbose=True)
    print(f"{oncovariant} p-value: {p_value}")
    print("")

print("...done.")
print("")
print("")
print("")

## Pairwise epistasis

db_epi = {} # for now, will be eventually just a DataFrame

for oncogene in oncogenes:
    db_epi[oncogene] = pd.read_csv(
        f"../data/mias_blair_response/ForLRT_oncogene/{oncogene}_samples_forLRT.csv",
        index_col=["variant_A", "variant_B", "data_type"])

for oncovariant in oncovariants:
    db_epi[oncovariant] = pd.read_csv(
        f"../data/mias_blair_response/ForLRT_oncovariant/{oncovariant.split('_')[1]}_samples_forLRT.csv",
        index_col=["variant_A", "variant_B", "data_type"])

db_epi["TSGinOncogene"] = pd.read_csv(
        f"../data/mias_blair_response/sampleSize_forLRT_TSGinOncogene.csv",
        index_col=["variant_A", "variant_B", "data_type"])

db_epi["TSGinOncovariant"] = pd.read_csv(
        f"../data/mias_blair_response/sampleSize_forLRT_TSGinOncovariant.csv",
        index_col=["variant_A", "variant_B", "data_type"])


## Merging into a single DataFrame
db_epi = pd.concat([db_epi[oncogene] for oncogene in oncogenes] +
                   [db_epi[oncovariant] for oncovariant in oncovariants] +
                   [db_epi["TSGinOncogene"]] + [db_epi["TSGinOncovariant"]])

## Ordering samples in the wanted order for estimate_lambdas:
## [(0, 0),
##  (0, 1),
##  (1, 0),
##  (1, 1)]
db_epi["samples"] = db_epi.apply(
    lambda row: np.array([row['n00'], row['nB0'], row['nA0'], row['nAB']]),
    axis=1)

db_epi = db_epi.sort_index()

## Computing substitution rates estimates

print("Computing substitution rates estimates under H1 "
      "(regular pair-wise epistasis models)...")
db_epi["lambdas_h1"] = db_epi["samples"].apply(
    lambda samples: np.nan if samples[-1] == 0
    else estimate_lambdas(samples, draws=1)['lambdas'])
print("...done.")
print("")
print("")
print("")


unique_Bs = list(db_epi.index.get_level_values('variant_B').unique())

p_values_under_h0_no_diff_smoking = {}
p_values_under_h0_no_selection = {}
lambdas_under_h0_no_diff_smoking = {}
lambdas_under_h0_no_selection = {}

print("Computing substitution rates estimates under H0 of same "
      "selection on smokers vs non-smokers and respective p-values...")
print("... for oncogenes...")
for variant in oncogenes:
    for tsg in unique_Bs:
        if ((variant, tsg, "Smo") not in db_epi.index or
            (variant, tsg, "nonSmo") not in db_epi.index or
            db_epi.loc[(variant, tsg, "Smo"), "nAB"].values[0] == 0 or
            db_epi.loc[(variant, tsg, "nonSmo"), "nAB"].values[0] == 0):
            print(f"Not estimating p-values for {(variant, tsg)} model")
        else:
            print(f"Estimating p-values for {(variant, tsg)} model")
            (p_values_under_h0_no_diff_smoking[(variant, tsg, "selectionA_with_B")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionA_with_B", "lambdas1")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionA_with_B", "lambdas2")]) = (
                 p_value_same_gamma_xy(
                     db_epi.loc[(variant, tsg, "Smo"), "samples"].values[0],
                     db_epi.loc[(variant, tsg, "nonSmo"), "samples"].values[0],
                     ((0, 1), (1, 1)),
                     mutation_rates.loc[variant.split("_")[0], "Smo"]["rate_grp_1"], # remove .split
                     mutation_rates.loc[variant.split("_")[0], "nonSmo"]["rate_grp_1"], # remove .split
                     lambdas1_h1=db_epi.loc[(variant, tsg, "Smo"), "lambdas_h1"].values[0],
                     lambdas2_h1=db_epi.loc[(variant, tsg, "nonSmo"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
            print("")
            (p_values_under_h0_no_diff_smoking[(variant, tsg, "selectionB_with_A")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionB_with_A", "lambdas1")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionB_with_A", "lambdas2")]) = (
                 p_value_same_gamma_xy(
                     db_epi.loc[(variant, tsg, "Smo"), "samples"].values[0],
                     db_epi.loc[(variant, tsg, "nonSmo"), "samples"].values[0],
                     ((1, 0), (1, 1)),
                     mutation_rates.loc[tsg, "Smo"]["rate_grp_1"],
                     mutation_rates.loc[tsg, "nonSmo"]["rate_grp_1"],
                     lambdas1_h1=db_epi.loc[(variant, tsg, "Smo"), "lambdas_h1"].values[0],
                     lambdas2_h1=db_epi.loc[(variant, tsg, "nonSmo"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
print("... oncongenes done...")

print("")
print("")
print("... for oncovariantes...")
for variant in oncovariants:
    for tsg in unique_Bs:
        if ((variant, tsg, "Smo") not in db_epi.index or
            (variant, tsg, "nonSmo") not in db_epi.index or
            db_epi.loc[(variant, tsg, "Smo"), "nAB"].values[0] == 0 or
            db_epi.loc[(variant, tsg, "nonSmo"), "nAB"].values[0] == 0):
            print(f"Not estimating p-values for {(variant, tsg)} model")
        else:
            print(f"Estimating p-values for {(variant, tsg)} model")
            (p_values_under_h0_no_diff_smoking[(variant, tsg, "selectionA_with_B")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionA_with_B", "lambdas1")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionA_with_B", "lambdas2")]) = (
                 p_value_same_gamma_xy(
                     db_epi.loc[(variant, tsg, "Smo"), "samples"].values[0],
                     db_epi.loc[(variant, tsg, "nonSmo"), "samples"].values[0],
                     ((0, 1), (1, 1)),
                     oncovariants_mutation_rates["smoking"][variant],
                     oncovariants_mutation_rates["nonsmoking"][variant],
                     lambdas1_h1=db_epi.loc[(variant, tsg, "Smo"), "lambdas_h1"].values[0],
                     lambdas2_h1=db_epi.loc[(variant, tsg, "nonSmo"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
            print("")
            (p_values_under_h0_no_diff_smoking[(variant, tsg, "selectionB_with_A")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionB_with_A", "lambdas1")],
             lambdas_under_h0_no_diff_smoking[(variant, tsg, "selectionB_with_A", "lambdas2")]) = (
                 p_value_same_gamma_xy(
                     db_epi.loc[(variant, tsg, "Smo"), "samples"].values[0],
                     db_epi.loc[(variant, tsg, "nonSmo"), "samples"].values[0],
                     ((1, 0), (1, 1)),
                     mutation_rates.loc[tsg, "Smo"]["rate_grp_1"],
                     mutation_rates.loc[tsg, "nonSmo"]["rate_grp_1"],
                     lambdas1_h1=db_epi.loc[(variant, tsg, "Smo"), "lambdas_h1"].values[0],
                     lambdas2_h1=db_epi.loc[(variant, tsg, "nonSmo"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
        print("")
        print("")
print("... oncovariants done...")




print("Computing substitution rates estimates under H0 "
      "of no selection (gamma = 1) and p-values...")
print("... for oncogenes...")
for variant in oncogenes:
    for tsg in unique_Bs:
        if (db_epi.loc[(variant, tsg, "All"), "nAB"].values[0] == 0):
            print(f"Not estimating p-values for {(variant, tsg)} model")
        else:
            print(f"Estimating p-values for {(variant, tsg)} model")
            (p_values_under_h0_no_selection[(variant, tsg, "selectionA_with_B")],
             lambdas_under_h0_no_selection[(variant, tsg, "selectionA_with_B", "lambdas")]) = (
                 p_value_gamma_xy_equal_1(
                     db_epi.loc[(variant, tsg, "All"), "samples"].values[0],
                     ((0, 1), (1, 1)),
                     mutation_rates.loc[variant, "All"]["rate_grp_1"],
                     lambdas_h1=db_epi.loc[(variant, tsg, "All"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
            print("")
            (p_values_under_h0_no_selection[(variant, tsg, "selectionB_with_A")],
             lambdas_under_h0_no_selection[(variant, tsg, "selectionB_with_A", "lambdas")]) =(
                 p_value_gamma_xy_equal_1(
                     db_epi.loc[(variant, tsg, "All"), "samples"].values[0],
                     ((1, 0), (1, 1)),
                     mutation_rates.loc[tsg, "All"]["rate_grp_1"],
                     lambdas_h1=db_epi.loc[(variant, tsg, "All"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
        print("")
        print("")
        print("")


print("... for oncovariants...")
for variant in oncovariants:
    for tsg in unique_Bs:
        if (db_epi.loc[(variant, tsg, "All"), "nAB"].values[0] == 0):
            print(f"Not estimating p-values for {(variant, tsg)} model")
        else:
            print(f"Estimating p-values for {(variant, tsg)} model")
            print(f"Selection of {variant} with {tsg} present")
            (p_values_under_h0_no_selection[(variant, tsg, "selectionA_with_B")],
             lambdas_under_h0_no_selection[(variant, tsg, "selectionA_with_B", "lambdas")]) = (
                 p_value_gamma_xy_equal_1(
                     db_epi.loc[(variant, tsg, "All"), "samples"].values[0],
                     ((0, 1), (1, 1)),
                     oncovariants_mutation_rates['allSamples'][variant],
                     lambdas_h1=db_epi.loc[(variant, tsg, "All"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
            print("")
            print(f"Selection of {tsg} with {variant} present")
            (p_values_under_h0_no_selection[(variant, tsg, "selectionB_with_A")],
             lambdas_under_h0_no_selection[(variant, tsg, "selectionB_with_A", "lambdas")]) =(
                 p_value_gamma_xy_equal_1(
                     db_epi.loc[(variant, tsg, "All"), "samples"].values[0],
                     ((1, 0), (1, 1)),
                     mutation_rates.loc[tsg, "All"]["rate_grp_1"],
                     lambdas_h1=db_epi.loc[(variant, tsg, "All"), "lambdas_h1"].values[0],
                     verbose=True,
                     return_lambdas_estimates=True))
        print("")
        print("")


## Significant results


print("Significant results for difference between smokers and non-smokers:")
for key, value in p_values_under_h0_no_diff_smoking.items():
    if value < 0.05:
        if key[-1] == "selectionA_with_B":
            print(f"Selection of mutation in {key[0]} when {key[1]} is already mutated")
            print(f"lambda smoking: {db_epi.loc[(key[0], key[1], 'Smo'), 'lambdas_h1'].values[0][-2]}")
            print(f"lambda non-smoking: {db_epi.loc[(key[0], key[1], 'nonSmo'), 'lambdas_h1'].values[0][-2]}")
        else:
            print(f"Selection of mutation in {key[1]} when {key[0]} is already mutated")
            print(f"lambda smoking: {db_epi.loc[(key[0], key[1], 'Smo'), 'lambdas_h1'].values[0][-1]}")
            print(f"lambda non-smoking: {db_epi.loc[(key[0], key[1], 'nonSmo'), 'lambdas_h1'].values[0][-1]}")
        print(f"p-value: {value}")
        print("")


print("Significant results for difference in selection after a mutation in another gene (epistasis):")
for key, value in p_values_under_h0_no_selection.items():
    if value < 0.05:
        if key[-1] == "selectionA_with_B":
            print(f"Selection of mutation in {key[0]} when {key[1]} is already mutated")
            print(f"lambda under H1: {db_epi.loc[(key[0], key[1], 'All'), 'lambdas_h1'].values[0][-2]}")
        else:
            print(f"Selection of mutation in {key[1]} when {key[0]} is already mutated")
            print(f"lambda under H1: {db_epi.loc[(key[0], key[1], 'All'), 'lambdas_h1'].values[0][-1]}")
        print(f"p-value: {value}")
        print("")


## Reorder dictionaries of p_values by p_value

p_values_under_h0_no_diff_smoking = sorted(
    p_values_under_h0_no_diff_smoking.items(), key=lambda item: item[1])


p_values_under_h0_no_selection = sorted(
    p_values_under_h0_no_selection.items(), key=lambda item: item[1])
