import os
import multiprocessing as mp

from math import comb
from math import ceil

from itertools import combinations
from itertools import chain

import pandas as pd
import numpy as np

from count_combinations import updated_compute_samples

from cancer_epistasis import estimate_lambdas
from cancer_epistasis import asymp_CI_lambdas
from cancer_epistasis import convert_lambdas_to_dict

from cancer_epistasis import compute_gammas
from cancer_epistasis import compute_CI_gamma

from theory import build_S_with_tuples
from theory import build_S_as_array
from theory import order_pos_lambdas

from locations import gene_list_file
from locations import location_output
from locations import results_keys
from locations import samples_per_combination_files

from load_results import load_results

from filter_data import key_filtered_dbs
from filter_data import filter_samples_for_genes

from pymc3.exceptions import SamplingError


# gene_list = list(pd.read_csv(gene_list_file, header=None)[0])
# gene_list = [gene.upper() for gene in gene_list]

gene_list = ["TP53","KRAS","EGFR","BRAF","CTNNB1",
             "KEAP1","STK11","ATM","PIK3CA","RBM10",
             "SMARCA4","SMAD4","ALK","ARID1A","APC",
             "MET","RB1","SETD2","BRCA2","MGA",
             "GNAS"] #"NRAS","ROS1","RET","U2AF1","NTRK",

# hm_df = pd.read_csv("../data/gene_sets/hallmark_pathway_df.csv")

# #HALLMARK_TNFA_SIGNALING_VIA_NFKB = list(hm_df[hm_df["pathway"]=="HALLMARK_TNFA_SIGNALING_VIA_NFKB"]["genes"])[0].split("|")
# HALLMARK_HYPOXIA = list(hm_df[hm_df["pathway"]=="HALLMARK_HYPOXIA"]["genes"])[0].split("|")
# HALLMARK_WNT_BETA_CATENIN_SIGNALING = list(hm_df[hm_df["pathway"]=="HALLMARK_WNT_BETA_CATENIN_SIGNALING"]["genes"])[0].split("|")
# HALLMARK_DNA_REPAIR = list(hm_df[hm_df["pathway"]=="HALLMARK_DNA_REPAIR"]["genes"])[0].split("|")
# HALLMARK_G2M_CHECKPOINT = list(hm_df[hm_df["pathway"]=="HALLMARK_G2M_CHECKPOINT"]["genes"])[0].split("|")

# HALLMARK_PI3K_AKT_MTOR_SIGNALING = list(hm_df[hm_df["pathway"]=="HALLMARK_PI3K_AKT_MTOR_SIGNALING"]["genes"])[0].split("|")
# HALLMARK_MYC_TARGETS_V2 = list(hm_df[hm_df["pathway"]=="HALLMARK_MYC_TARGETS_V2"]["genes"])[0].split("|")
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = list(hm_df[hm_df["pathway"]=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]["genes"])[0].split("|")
# HALLMARK_INFLAMMATORY_RESPONSE = list(hm_df[hm_df["pathway"]=="HALLMARK_INFLAMMATORY_RESPONSE"]["genes"])[0].split("|")
# HALLMARK_METABOLISM = list(hm_df[hm_df["pathway"]=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"]["genes"])[0].split("|") + list(hm_df[hm_df["pathway"]=="HALLMARK_GLYCOLYSIS"]["genes"])[0].split("|")

# HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = list(hm_df[hm_df["pathway"]=="HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]["genes"])[0].split("|")
# HALLMARK_P53_PATHWAY = list(hm_df[hm_df["pathway"]=="HALLMARK_P53_PATHWAY"]["genes"])[0].split("|")
# HALLMARK_ANGIOGENESIS = list(hm_df[hm_df["pathway"]=="HALLMARK_ANGIOGENESIS"]["genes"])[0].split("|")
# HALLMARK_KRAS = list(hm_df[hm_df["pathway"]=="HALLMARK_KRAS_SIGNALING_UP"]["genes"])[0].split("|") + list(hm_df[hm_df["pathway"]=="HALLMARK_KRAS_SIGNALING_DN"]["genes"])[0].split("|")

# gene_list = {"HALLMARK_HYPOXIA":HALLMARK_HYPOXIA,
#              "HALLMARK_WNT_BETA_CATENIN_SIGNALING":HALLMARK_WNT_BETA_CATENIN_SIGNALING,
#              "HALLMARK_DNA_REPAIR":HALLMARK_DNA_REPAIR,
#              "HALLMARK_G2M_CHECKPOINT":HALLMARK_G2M_CHECKPOINT,
#              "HALLMARK_PI3K_AKT_MTOR_SIGNALING":HALLMARK_PI3K_AKT_MTOR_SIGNALING,
#              "HALLMARK_MYC_TARGETS_V2":HALLMARK_MYC_TARGETS_V2,
#              "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION":HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
#              "HALLMARK_INFLAMMATORY_RESPONSE":HALLMARK_INFLAMMATORY_RESPONSE,
#              "HALLMARK_METABOLISM":HALLMARK_METABOLISM,
#              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY":HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
#              "HALLMARK_P53_PATHWAY":HALLMARK_P53_PATHWAY,
#              "HALLMARK_ANGIOGENESIS":HALLMARK_ANGIOGENESIS,
#              "HALLMARK_KRAS":HALLMARK_KRAS}

n_cores = 50

def filter_and_compute_samples(combo, key, pathways=False, print_info=False):
    if pathways:
        all_genes = list(chain(*combo.values()))
        db = filter_samples_for_genes(all_genes, key_filtered_dbs[key], print_info=print_info)
        for pathway, genes in combo.items():
            db[f'{pathway}']=db[genes].sum(axis='columns').apply(lambda x: 1 if x > 1 else x)
        counts = updated_compute_samples(db,
                                         mutations=list(combo.keys()),
                                         print_info=True)
        combo = tuple(combo.keys())
    else:
        db = filter_samples_for_genes(combo, key_filtered_dbs[key], print_info=print_info)
        counts = updated_compute_samples(db,
                                         mutations=list(combo),
                                         print_info=print_info)


    return {combo: counts}


def compute_samples_for_all_combinations(genes=None,
                                         key=None,
                                         num_per_combo=3,
                                         mus_keys=None,
                                         pathways=False,
                                         chunksize=None,
                                         print_info=False,
                                         save_results=True):

    """Compute number of patients with each combination of
    const:`genes`.

    :type genes: list or NoneType
    :param genes: List of genes from which gene-set
        combinations will be made.


    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking
        - smoking_plus (includes panel data)
        - nonsmoking_plus (includes panel data)

    :type num_per_combo: int or list
    :param num_per_combo: How many genes to include in
        each set of genes. It can also be list with integers.

    :type save_results: bool
    :param save_results: If True (default) save results as a csv file.

    :rtype: dict
    :return: Dictionary with keys being the gene-sets, and values
        being arrays with number of samples per mutation combination,
        ordered as in :func:`build_S_with_tuples`.

    """

    if key is None:
        key = "pan_data"
    if genes is None:
        raise IOError("Please input a list of genes.")

    if pathways:
        all_genes = list(chain(*genes.values()))
    else:
        all_genes = genes

    if isinstance(num_per_combo, list):
        counts = {}
        for M in num_per_combo:
            print(f"Computing all samples with M={M}...")
            counts.update(compute_samples_for_all_combinations(
                genes=genes,
                key=key,
                num_per_combo=M,
                mus_keys=mus_keys,
                pathways=pathways,
                print_info=print_info,
                save_results=False)) # results will be saved at the
                                     # end but in a single dictionary

            print("")
            print("")


    else:

        if not isinstance(num_per_combo, int):
            raise ValueError("The number of genes in each combination must be an integer.")
        print("Computing samples for all combinations of "
              f"{num_per_combo} from {len(genes)} genes/pathways for {key}...")
        print("\n")
        unrepresented_genes = [gene for gene in all_genes if gene not in key_filtered_dbs[key].columns]
        if len(unrepresented_genes) > 0:
            print(f"No sample availability information for the following genes: "
                  f"{str(unrepresented_genes)}."
                  "\nThis may be because there are not mutations in those genes in the data, "
                  "or because the genes are not correctly represented in sequencing data")
        no_mutation_rate_genes = [gene for gene in all_genes if gene not in mus_keys]
        if len(no_mutation_rate_genes) > 0:
            print(f"No mutation rate available for the following genes: "
                  f"{str(no_mutation_rate_genes)}."
                  "\nThis may be a problem with the cancereffectsizeR mutation rate estimation")
        genes_to_remove = set(unrepresented_genes + no_mutation_rate_genes)
        print("\n")

        if pathways:
            pathway_prop_missing = {
                pathway: (len(set.intersection(set(pathway_genes), genes_to_remove))
                          / len(pathway_genes)) # proportion of missing genes in each pathway
                for pathway, pathway_genes in genes.items()}
            print("Proportion of genes in each "
                  f"pathway for which we can't calculate fluxes: {pathway_prop_missing}")
            pathways_to_remove = [pathway for pathway, prop_missing in pathway_prop_missing.items()
                                  if prop_missing > 0.05]
            if len(pathways_to_remove) > 0:
                print("\n")
                print("Dropping the following pathways "
                      "because more than 5% of genes have incalculable fluxes: "
                      f"{pathways_to_remove}")

            for pathway in pathways_to_remove:
                genes.pop(pathway)
            genes = {pathway: list(filter(lambda gene: gene not in genes_to_remove, pathway_genes))
                        for pathway, pathway_genes in genes.items()}
        else:
            genes = list(filter(lambda gene: gene not in genes_to_remove, genes))

        pool = mp.Pool(processes=n_cores)
        gene_combos = list(combinations(genes, num_per_combo))

        if pathways:
            mp_results = pool.starmap(filter_and_compute_samples,
                                      [({pathway: genes[pathway] for pathway in combo},
                                        key,
                                        pathways,
                                        print_info) for combo in gene_combos],
                                      chunksize=chunksize)
        else:
            mp_results = pool.starmap(filter_and_compute_samples,
                                      [(combo,
                                        key,
                                        pathways,
                                        print_info) for combo in gene_combos],
                                      chunksize=chunksize)

        counts = {combo: count_array
                  for result in mp_results
                  for combo, count_array in result.items()}


    if save_results:
        print("Saving results...")
        np.save(os.path.join(location_output,
                             samples_per_combination_files[key]),
                counts)

        print("done.")

    return counts


def are_all_fluxes_computable(samples):
    """Return True if we can compute all fluxes, that is, if there are
    patients in `data` with every possible combination of `genes`
    (except for the combination of all genes)

    There is a version of this function in the
    :mod:`count_combinations` module, but it this one is way faster if
    the samples per combination have already been computed.

    :type samples: numpy.ndarray
    :param samples: One dimensional array with the samples. It should
        be of the same size as S, and have in each entry the number of
        individuals that have the respective mutation combination.

    :rtype: bool
    :return: True if all fluxes are computable.

    """
    return np.all(samples[:-1] > 0)


def all_but_last_layer_computable(samples, M):
    S = build_S_as_array(M)
    necessary_pos_ind = np.where(np.sum(S, axis=1) < M-1)

    return np.all(samples[necessary_pos_ind] > 0)


def lambdas_from_samples(samples, max_bound_changes=4):
    """Estimate the fluxes from the samples, varying the bounds if
    necessary.

    :type samples: numpy.ndarray
    :param samples: One dimensional array with the samples. It should
        be of the same size as S, and have in each entry the number of
        individuals that have the respective mutation combination.

    :type max_bound_changes: int
    :param max_bound_changes: If the estimating algorithm does not
        converge using the default bounds ([0, 1]) for the uniform
        prior use for each of the lambdas, then those bounds are
        changed for the particular lambdas that we have identified as
        problematic with those bounds for certain models. We change
        the upper bound to half of the previous bound. This variable
        is the maximum number of times that we will try to change
        those bounds before concluding that we cannot compute the
        model for those samples (likely they are too few samples for
        the gene combination).

    :rtype: dict
    :return: Estimates of the fluxes as a dictionary with the MAP
        (also MLE, because the priors are uniform) estimates.

    """

    bounds = 1
    bound_maxes = np.array([bounds])

    print(f"Bounds for fluxes: {bounds}")
    draws = 1
    MLE = estimate_lambdas(samples, draws=draws,
                           upper_bound_prior=bounds,
                           kwargs={'return_raw':True})
    if len(samples) == 2 and draws == 1: # M=1 model with one draw
        print(f"MLE: {MLE['lambdas']}")
        return MLE
    else:
        print(f"MLE: {MLE[0]['lambdas']}")

    bound_changes = 0

    # conditions necessitating change of bounds
    if any(x in bound_maxes for x in MLE[0]['lambdas']) or MLE[1].fun < -1e+20:
        M = int(np.log2(len(samples)))
        num_fluxes = M*2**(M-1) # formula for number of fluxes
        bounds = np.array([bounds]*num_fluxes, dtype='float')
    # if any of the fluxes reach the upper bound, common for pathway analysis
    while any(x in bound_maxes for x in MLE[0]['lambdas']):
        increment = 1
        if bound_changes == max_bound_changes:
            print(f"We have changed the bounds {bound_changes} times already, "
                  "and the algorithm has not converged. Concluding that these "
                  "samples are incomputable.")
            return "incomputable"
        print("Upper bound hit for one or more fluxes, increasing bounds...")
        bound_changes += 1
        to_increase = [i for i in range(num_fluxes) if MLE[0]['lambdas'][i] in bound_maxes]
        bounds = [bounds[i]+increment if i in to_increase else bounds[i] for i in range(num_fluxes)]
        bound_maxes = np.append(bound_maxes, bound_maxes[-1]+increment)

        print(f"Proposed bounds: {bounds}")
        try:
            MLE = estimate_lambdas(samples, draws=draws,
                                upper_bound_prior=bounds,
                                kwargs={'return_raw':True})
        except SamplingError:
            return "incomputable"
        print(f"MLE: {MLE[0]['lambdas']}")

    while MLE[1].fun < -1e+20:
        if bound_changes == max_bound_changes:
            print(f"We have changed the bounds {bound_changes} times already, "
                  "and the algorithm has not converged. Concluding that these "
                  "samples are incomputable.")
            return "incomputable"
        print("Algorithm did not converge changing bounds...")
        bound_changes += 1

        # when the flux is below an arbitarily chosen value, lower the upper bound to optimize the search
        to_decrease = np.where(MLE[0]['lambdas'] < 1e-5)[0]
        for i in to_decrease:
            bounds[i] = bounds[i]/2
        # bounds = [bounds[i]/2 if i in to_decrease else bounds[i] for i in range(num_fluxes)]

        print(f"Proposed bounds: {bounds}")
        MLE = estimate_lambdas(samples, draws=draws,
                               upper_bound_prior=bounds,
                               kwargs={'return_raw':True})
        print(f"MLE: {MLE[0]['lambdas']}")

    print("")

    return MLE[0]


def compute_lambda_for_combo(combo, counts, flexible_last_layer):
    samples = counts[combo]

    if are_all_fluxes_computable(samples):
        print(f"Estimating fluxes for {combo}...")
        mle = lambdas_from_samples(samples)

        if(mle == 'incomputable'):
            with open(os.path.join(location_output,"incomputable_limited_selection_fluxes.txt"),'a', encoding='utf-8') as incomputable_output_file:
                incomputable_output_file.write(str(combo) + '\n')
            return None

        combo_lambda_mles = convert_lambdas_to_dict(mle)

        print("Estimating asymptotic confidence intervals...")
        cis = asymp_CI_lambdas(mle['lambdas'], samples)
        combo_lambda_cis = convert_lambdas_to_dict(cis)

        return combo, combo_lambda_mles, combo_lambda_cis

    elif flexible_last_layer and all_but_last_layer_computable(samples, len(combo)):
        print(f"Estimating a limited selection of fluxes for {combo}...")
        mle = lambdas_from_samples(samples)

        if(mle == 'incomputable'):
            with open(os.path.join(location_output,"incomputable_limited_selection_fluxes.txt"),'a', encoding='utf-8') as incomputable_output_file:
                incomputable_output_file.write(str(combo) + '\n')
            return None

        S = build_S_as_array(len(combo))
        # Remove fluxes from second-to-last layer (all but 1 mutation) to last layer (all mutations)
        states_with_zero = [tuple(x) for x in S[np.where(np.sum(S, axis=1) >= len(combo)-1)]]
        indices_with_zero = [order_pos_lambdas(S).index((x, tuple(y)))
                                for x, y in order_pos_lambdas(S)
                                if x in states_with_zero]
        for x in indices_with_zero:
            mle['lambdas'][x] = np.nan

        combo_lambda_mles = convert_lambdas_to_dict(mle)

        print("Estimating asymptotic confidence intervals...")
        cis = asymp_CI_lambdas(mle['lambdas'], samples)
        combo_lambda_cis = convert_lambdas_to_dict(cis)

        return combo, combo_lambda_mles, combo_lambda_cis

    else:
        print(f"Skipping estimation for combination {combo} "
                "because the fluxes of interest are not "
                "computable for that model.")
    print("")

    return None


def compute_all_lambdas(key, all_counts, flexible_last_layer=False, chunksize=None, save_results=True):
    """Compute all estimates of the fluxes for the data set `key`
    iterating over all genes in :const:`gene_list`.in

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking

    :type all_counts: dict
    :param all_counts: Dictionary with keys being the results keys,
        and values being dictionaries that contain the samples per
        mutation combination as return by :func:`compute_samples` for
        each of genes in :const:`gene_list`.

    :type save_results: bool
    :param save_results: If True (default) save results.

    :rtype: tuple
    :return: A tuple with the maximum likelihood estimations and the
        95% asymptomatic confidence intervals for the fluxes.

    """
    lambdas_mles = {}
    lambdas_cis = {}

    counts = all_counts[key]

    pool = mp.Pool(processes=n_cores)
    mp_results = pool.starmap(compute_lambda_for_combo,
                              [(combo, counts, flexible_last_layer) for combo in counts.keys()],
                              chunksize=chunksize)
    lambdas_mles = {result[0]: result[1] for result in mp_results if result is not None}
    lambdas_cis = {result[0]: result[2] for result in mp_results if result is not None}

    if save_results:
        np.save(os.path.join(location_output,
                                f"{key}_fluxes_mles.npy"),
                lambdas_mles)
        np.save(os.path.join(location_output,
                                f"{key}_fluxes_cis.npy"),
                lambdas_cis)

    return lambdas_mles, lambdas_cis


def compute_all_gammas(key, all_lambdas, mus, pathways=False, pathway_genes_dict=None, save_results=True):
    """Compute all estimates of the selection coefficient for the data
    set `key` iterating over all gene combinations that are present in
    the keys of :dict:`all_lambdas`

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking
        - smoking_plus
        - nonsmoking_plus

    :type all_lambdas: dict
    :param all_lambdas: Dictionary with the results for fluxes as
        obtained from :func:`compute_all_lambdas`. It is indexed by
        tuples containing the a key and either 'mles' or 'cis'

    :type mus: dict
    :param mus: Dictionary with the mutation rates for each gene.

    :type save_results: bool
    :param save_results: If True (default) save results.

    :rtype: tuple
    :return: A tuple with the maximum likelihood estimations and the
        95% asymptomatic confidence intervals for the fluxes.

    """
    gammas_mles = {}
    gammas_cis = {}

    lambdas_mles = all_lambdas[key, 'mles']
    lambdas_cis = all_lambdas[key, 'cis']

    total_num_combos = len(lambdas_mles.keys())

    for counter, combo in enumerate(lambdas_mles.keys()):

        print(f"Estimating gammas for combination {counter+1}/{total_num_combos}")

        M = len(combo)
        # Constructs dictionary of the form:
        # mu_combo = {combo:{(1, 0, 0): mus[combo[0]],
        #                    (0, 1, 0): mus[combo[1]],
        #                    (0, 0, 1): mus[combo[2]]}}
        # if M == 3

        if not pathways:
            mu_combo = {(i*(0,) + (1,) + (M-i-1)*(0,)): mus[combo[i]]
                for i in range(M)}
        else:
            subset_pathway_genes_dict = {pathway:pathway_genes_dict[pathway] for pathway in combo}
            gene_indices = [i for i, included in enumerate(subset_pathway_genes_dict.values()) if len(included) == 1]
            mu_combo = {(i*(0,) + (1,) + (M-i-1)*(0,)): mus[list(subset_pathway_genes_dict.values())[i][0]] if i in gene_indices
                        else sum([mus[gene] for gene in list(subset_pathway_genes_dict.values())[i] if gene in mus.keys()])
                    for i in range(M)}

        gammas_mles[combo] = compute_gammas(
            lambdas_mles[combo],
            mu_combo)

        gammas_cis[combo] = compute_CI_gamma(
            lambdas_cis[combo],
            mu_combo)

        if save_results:
            np.save(os.path.join(location_output,
                                 f"{key}_selections_mles.npy"),
                    gammas_mles)
            np.save(os.path.join(location_output,
                                 f"{key}_selections_cis.npy"),
                    gammas_cis)

    return gammas_mles, gammas_cis


def main(genes=None,
         num_per_combo=[1,2,3],
         keys=None,
         mu_method="variant",
         pathways=False,
         flexible_last_layer=False,
         recompute_samples_per_combination=False,
         print_info=True,
         save_results=True):
    """Main method for the estimation of all the fluxes.

    :type genes: list, dict, or NoneType
    :param genes: Either a list of genes from which gene-set combinations will
        be made or a dictionary with the genes or pathways as keys and the gene
        or the genes in the pathway, respectively, as values.

    :type num_per_combo: int
    :param num_per_combo: How many genes to include in each set of
        genes

    :type recompute_samples_per_combination: bool
    :param recompute_samples_per_combination: If True force
        recomputing the samples per combination of each
        gene. Otherwise (default) try load the respective file if
        available.

    :type save_results: bool
    :param save_results: If True (default) save results.

    :rtype: tuple
    :return: A tuple with three dictionaries that contain all the
        samples per combination, all the fluxes estimated and all the
        scaled selection coefficients. The keys of the dictionaries
        are the keys in :const:`results_keys`.

    """
    all_counts = {}
    all_lambdas = {}
    all_gammas = {}

    if genes is None:
        genes = gene_list
    if keys is None:
        keys = results_keys

    if isinstance(genes, list):
        if pathways:
            raise Warning("A list was detected even though `pathways` was set to true."
                          "If you are evaluating pathways, please provide a dictionary,"
                          "with the pathway names as keys and the genes in the pathways"
                          "as values")
        pathways = False
    elif isinstance(genes, dict):
        pathways = True
        if not all([isinstance(entry, list) for entry in genes.values()]):
            raise ValueError("When including pathways, all entries in "
                             "`genes` must be lists, not strings, "
                             "even if there is only gene for the pathway")
    else:
        raise IOError("`genes` must be either a list of genes "
                      "or a dictionary of genes and pathways.")

    if isinstance(num_per_combo, int):
        num_per_combo = [num_per_combo]

    chunksize = ceil(1/3 * sum([comb(len(genes),k) for k in num_per_combo]) / n_cores)

    for key in keys:
        print("")

        mus = load_results('mutations', mu_method)[key]

        if (recompute_samples_per_combination
            or not os.path.exists(samples_per_combination_files[key])):
            print(f"Computing number of samples per combination for {key}...")
            all_counts[key] = compute_samples_for_all_combinations(genes, key, num_per_combo, mus.keys(),
                                                                   pathways, chunksize, print_info, save_results)
        else:
            print(f"Loading counts per combination for {key}...")
            all_counts[key] = np.load(samples_per_combination_files[key],
                                      allow_pickle=True).item()
        print(f"done computing samples per combination for {key}.")
        print("")
        print("")

        print(f"Estimating all epistatic models for {key}...")
        print("")
        lambdas_mles, lambdas_cis = compute_all_lambdas(key, all_counts, flexible_last_layer, chunksize, save_results)
        all_lambdas[(key, 'mles')] = lambdas_mles
        all_lambdas[(key, 'cis')] = lambdas_cis
        print(f"done estimating all epistatic models for {key}.")
        print("")
        print("")

        print(f"Computing selection coefficients for {key}...")
        print("")
        if pathways:
            gammas_mles, gammas_cis = compute_all_gammas(key, all_lambdas, mus, pathways, genes, save_results)
        else: # if just genes
            gammas_mles, gammas_cis = compute_all_gammas(key, all_lambdas, mus, pathways, save_results)
        all_gammas[(key, 'mles')] = gammas_mles
        all_gammas[(key, 'cis')] = gammas_cis
        print(f"done computing selection coefficients for {key}.")
        print("")
        print("")


    return all_counts, all_lambdas, all_gammas



if __name__ == "__main__":
    print("Running main...")
    print("")
    main(recompute_samples_per_combination=True,
         flexible_last_layer=False,
         pathways=False,
         num_per_combo=[1, 2, 3],
         mu_method="variant",
         keys = results_keys)
    print("")
    print('Done running main.')
