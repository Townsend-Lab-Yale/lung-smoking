import os
import pandas as pd
import numpy as np
import multiprocessing as mp

from itertools import combinations

from count_combinations import updated_compute_samples
from count_combinations import convert_samples_to_dict

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

from filter_data import key_filtered_dbs
from filter_data import dbs_filtered_for_TP53_KRAS
from filter_data import filter_samples_for_genes



# gene_list = list(pd.read_csv(gene_list_file, header=None)[0])
# gene_list = [gene.upper() for gene in gene_list]
# gene_list = gene_list[:103]

gene_list = ["TP53","KRAS","PIK3CA","BRAF","RB1","STK11","KEAP1","EGFR","CDKN2A.p16INK4a"]

def filter_and_compute_samples(combo, key):
    print(combo)
    db = filter_samples_for_genes(combo, key_filtered_dbs[key], print_info=True)
    counts = updated_compute_samples(db,
                                    mutations = list(combo),
                                    print_info = True)
    return {combo: counts}
    

def compute_samples_for_all_combinations(genes = None, key=None, num_per_combo = 3, save_results=True):
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
    
    :type num_per_combo: int
    :param num_per_combo: How many genes to include in
        each set of genes

    :type save_results: bool
    :param save_results: If True (default) save results as a csv file.

    :rtype: dict
    :return: Dictionary with keys being the gene-sets, and values being
        arrays with number of samples per mutation combination as
        return by :func:`compute_samples`.

    """

    if key is None:
        key = "pan_data"
    if genes is None:
        raise IOError("Please input a list of genes.")
    if num_per_combo < 2 or not isinstance(num_per_combo, int):
        raise ValueError("The number of genes in each combinations must be an integer greater than 1.")
    print(f"Computing samples for all combinations of {num_per_combo} from {len(genes)} genes.")

    unrepresented_genes = [gene for gene in genes if gene not in key_filtered_dbs[key].columns]
    if len(unrepresented_genes) > 0:
        print(f"No sample availability information for the following genes: "
              f"{str(unrepresented_genes)}."
              "\nThis may be because there are not mutations in those genes in the data, "
              "or because the genes are not correctly represented in sequencing data")

    genes = list(filter(lambda gene: gene not in unrepresented_genes, genes))

    pool = mp.Pool(processes=8)
    gene_combos = list(combinations(genes, num_per_combo))
    mp_results = pool.starmap(filter_and_compute_samples, 
                              [(gene_combos[i], key) for i in range(len(gene_combos))], 
                              chunksize=8)
    counts = {combo: count_array for result in mp_results for combo, count_array in result.items()}
    # print(counts)
    # for combo in combinations(genes, num_per_combo):
    #     db = filter_samples_for_genes(combo, key_filtered_dbs[key], print_info=True)
    #     counts[combo] = updated_compute_samples(db,
    #                                             mutations = list(combo),
    #                                             print_info = True)

    if save_results:
        print("Saving results...")
        df = pd.DataFrame.from_dict(counts, orient='index',
                                    columns=[str(x) for x in build_S_with_tuples(3)])
        df.index.name = "gene combination"
        df.to_csv(samples_per_combination_files[key])

        print("done.")

    return counts

def compute_samples_for_TP53_KRAS_gene_model(key=None, save_results=True):
    """Compute number of patients with each combination of the TP53,
    KRAS and a third gene, for each third gene in the
    const:`gene_list`.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking

    :type save_results: bool
    :param save_results: If True (default) save results as a csv file.

    :rtype: dict
    :return: Dictionary with keys being the genes, and values being
        arrays with number of samples per mutation combination as
        return by :func:`compute_samples`.

    """

    if key is None:
        key = "pan_data"

    genes = gene_list[2:]

    unrepresented_genes = []

    counts = {}

    for gene in genes:
        db = filter_samples_for_genes(gene, dbs_filtered_for_TP53_KRAS[key],print_info=True)

        if gene in db.columns:
            counts[gene] = updated_compute_samples(db,
                                        mutations=['TP53', 'KRAS', gene])
        else:
            unrepresented_genes.append(gene)
            counts[gene] = np.repeat(0, 8)
    if len(unrepresented_genes) > 0:
        print(f"No sample availability information for the following genes: "
              f"{str(unrepresented_genes)}."
              "\nThis may be because there are not mutations in those genes in the data, "
              "or because the genes are not correctly represented in sequencing data")

    if save_results:
        print("Saving results...")
        df = pd.DataFrame.from_dict(counts, orient='index',
                                    columns=[str(x) for x in build_S_with_tuples(3)])
        df.index.name = "third gene"
        df.to_csv(samples_per_combination_files[key])

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


def at_least_000_to_001_and_110_to_111(samples):
    """Return True if the flux from wildtype to only the third gene and
    the flux from TP53 + KRAS to TP53 + KRAS + third gene are both computable"""
    return np.all(samples[[0,-2]] > 0)

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
    print(f"Bounds for fluxes: {bounds}")
    MLE = estimate_lambdas(samples, draws=1,
                           upper_bound_prior=bounds,
                           kwargs={'return_raw':True})
    print(f"MLE: {MLE[0]['lambdas']}")

    bound_changes = 0
    while MLE[1].fun < -1e+20:
        if bound_changes == max_bound_changes:
            print(f"We have changed the bounds {bound_changes} times already, "
                  "and the algorithm has not converged. Concluding that these "
                  "samples are incomputable.")
            return "incomputable"
        print(f"Algorithm did not converge changing bounds...")
        bound_changes += 1
        bounds = np.array(
            [1/(2**bound_changes),  # (0, 0, 0) -> (0, 0, 1)
             1,                     # (0, 0, 0) -> (0, 1, 0)
             1,                     # (0, 0, 1) -> (0, 1, 1)
             1/(2**bound_changes),  # (0, 1, 0) -> (0, 1, 1)
             1,                     # (0, 0, 0) -> (1, 0, 0)
             1,                     # (0, 0, 1) -> (1, 0, 1)
             1/(2**bound_changes),  # (1, 0, 0) -> (1, 0, 1)
             1,                     # (0, 1, 0) -> (1, 1, 0)
             1,                     # (1, 0, 0) -> (1, 1, 0)
             1,                     # (0, 1, 1) -> (1, 1, 1)
             1,                     # (1, 0, 1) -> (1, 1, 1)
             1/(2**bound_changes)]) # (1, 1, 0) -> (1, 1, 1)
        print(f"Proposed bounds: {bounds}")
        MLE = estimate_lambdas(samples, draws=1,
                               upper_bound_prior=bounds,
                               kwargs={'return_raw':True})
        print(f"MLE: {MLE[0]['lambdas']}")

    return MLE[0]

def compute_lambda_for_combo(combo, counts, flexible_last_layer):
    samples = counts[combo]

    if (are_all_fluxes_computable(samples)) or (flexible_last_layer and all_but_last_layer_computable(samples, len(combo))):
        print(f"Estimating fluxes for {combo}...")
        mle = lambdas_from_samples(samples)
        combo_lambda_mles = convert_lambdas_to_dict(mle)

        print("Estimating asymptotic confidence intervals...")
        cis = asymp_CI_lambdas(mle['lambdas'], samples)
        combo_lambda_cis = convert_lambdas_to_dict(cis)

        return combo, combo_lambda_mles, combo_lambda_cis

    else:
        print(f"Skipping estimation for combination {combo} "
                "because the fluxes of interest are not "
                "computable for that model.")

    
def compute_all_lambdas(key, all_counts, flexible_last_layer=False, save_results=True):
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

    pool = mp.Pool(processes=8)
    mp_results = pool.starmap(compute_lambda_for_combo, 
                              [(combo, counts, flexible_last_layer) for combo in counts.keys()], 
                              chunksize=8)
    print(mp_results)
    lambdas_mles = {result[0]: result[1] for result in mp_results if result is not None}
    lambdas_cis = {result[0]: result[2] for result in mp_results if result is not None}

    print(lambdas_mles)
    print(lambdas_cis)

    if save_results:
        np.save(os.path.join(location_output,
                                f"{key}_fluxes_mles.npy"),
                lambdas_mles)
        np.save(os.path.join(location_output,
                                f"{key}_fluxes_cis.npy"),
                lambdas_cis)




    # for i, combo in enumerate(counts.keys()):
    #     samples = counts[combo]
    #     if are_all_fluxes_computable(samples):
    #         print(f"Running model with {combo} "
    #               f"(combo number {i+1}/{len(counts)} "
    #               f"for {key})")

    #         print("Estimating fluxes...")
    #         mle = lambdas_from_samples(samples)
    #         lambdas_mles[combo] = convert_lambdas_to_dict(mle)
    #         if save_results:
    #             np.save(os.path.join(location_output,
    #                                     f"{key}_fluxes_mles.npy"),
    #                     lambdas_mles)

    #         print("Estimating asymptotic confidence intervals...")
    #         cis = asymp_CI_lambdas(mle['lambdas'], samples)
    #         lambdas_cis[combo] = convert_lambdas_to_dict(cis)
    #         if save_results:
    #             np.save(os.path.join(location_output,
    #                                     f"{key}_fluxes_cis.npy"),
    #                     lambdas_cis)

    #     else:
    #         print(f"Skipping estimation for combination {combo} "
    #                 "because the fluxes of interest are not "
    #                 "computable for that model.")

    #     print("")

    return lambdas_mles, lambdas_cis


def compute_all_lambdas_for_TP53_KRAS_gene_model(key, all_counts, save_results=True):
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

    S_3 = build_S_as_array(3)

    for i, gene in enumerate(gene_list[2:]):
        samples = counts[gene]
        if are_all_fluxes_computable(samples):
            print(f"Running model with third gene {gene} "
                  f"(gene number {i+1}/{len(gene_list[2:])}"
                  f"for {key})")

            print("Estimating fluxes...")
            mle = lambdas_from_samples(samples)
            lambdas_mles[gene] = convert_lambdas_to_dict(mle)
            if save_results:
                np.save(os.path.join(location_output,
                                        f"{key}_fluxes_mles.npy"),
                        lambdas_mles)

            print("Estimating asymptotic confidence intervals...")
            cis = asymp_CI_lambdas(mle['lambdas'], samples)
            lambdas_cis[gene] = convert_lambdas_to_dict(cis)
            if save_results:
                np.save(os.path.join(location_output,
                                        f"{key}_fluxes_cis.npy"),
                        lambdas_cis)

        elif at_least_000_to_001_and_110_to_111(samples):
            print(f"Running model with third gene {gene} "
                    f"(gene number {i+1}/{len(gene_list[2:])})")

            print("Estimating a limited selection of fluxes...")
            mle = lambdas_from_samples(samples)
            if(mle == 'incomputable'):
                with open(os.path.join(location_output,"incomputable_limited_selection_fluxes.txt"),'a') as incomputable_output_file:
                    incomputable_output_file.write(gene + '\n')
                continue

            states_with_zero = [tuple(x) for x in S_3[samples == 0]]
            indices_with_zero = [order_pos_lambdas(S_3).index((x, tuple(y)))
                                 for x, y in order_pos_lambdas(S_3)
                                 if x in states_with_zero]
            for x in indices_with_zero: mle['lambdas'][x] = np.nan

            lambdas_mles[gene] = convert_lambdas_to_dict(mle)

            if save_results:
                np.save(os.path.join(location_output,
                                        f"{key}_fluxes_mles.npy"),
                        lambdas_mles)

            print("Estimating asymptotic confidence intervals...")
            cis = asymp_CI_lambdas(mle['lambdas'], samples)
            lambdas_cis[gene] = convert_lambdas_to_dict(cis)
            if save_results:
                np.save(os.path.join(location_output,
                                        f"{key}_fluxes_cis.npy"),
                        lambdas_cis)

        else:
            print(f"Skipping estimation for gene {gene} "
                    "because the fluxes of interest are not "
                    "computable for that model.")

        print("")

    return lambdas_mles, lambdas_cis


def compute_TP53_KRAS_model(key='pan_data', save_results=True):
    """Compute a model for M=2 that only includes TP53 and KRAS.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking

    :type save_results: bool
    :param save_results: If True (default) save results.

    :rtype: dict
    :return: A dictionary with all results (lambdas, mus and gammas).

    """
    results = {}

    samples = pd.read_csv(samples_per_combination_files[key],
                          index_col='third gene')

    samples = samples.loc['EGFR'] ## these ones are for 3 gene model

    samples_00 = samples["(0, 0, 0)"] + samples["(0, 0, 1)"]
    samples_01 = samples["(0, 1, 0)"] + samples["(0, 1, 1)"]
    samples_10 = samples["(1, 0, 0)"] + samples["(1, 0, 1)"]
    samples_11 = samples["(1, 1, 0)"] + samples["(1, 1, 1)"]

    samples = [samples_00, samples_01, samples_10, samples_11]

    results['samples'] = {(0, 0):samples_00,
                          (0, 1):samples_01,
                          (1, 0):samples_10,
                          (1, 1):samples_11}

    bounds = np.repeat(1,4)
    MLE = estimate_lambdas(samples, draws=1,
                           upper_bound_prior=bounds,
                           kwargs={'return_raw':True})[0]

    results[('lambdas', 'mle')] = convert_lambdas_to_dict(MLE)

    results[('lambdas', 'cis')] = convert_lambdas_to_dict(
        asymp_CI_lambdas(MLE['lambdas'], samples))

    if key[-4:] == 'plus':
        final_part_key = 'w_panel'
    else:
        final_part_key = key[-4:]
    mus = pd.read_csv(
        os.path.join(location_output,
                     f"{key[:-4] + final_part_key}_mutation_rates.txt"),
        index_col='gene')

    if key == 'pan_data':
        mus = {(1, 0): float(mus.loc['TP53']),
               (0, 1): float(mus.loc['KRAS'])}
    else:
        mus = {(1, 0): mus.loc['TP53', 'rate_grp_1'],
               (0, 1): mus.loc['KRAS', 'rate_grp_1']}
    results['mus'] = mus

    results[('gammas', 'mle')] = compute_gammas(results[('lambdas', 'mle')], mus)
    results[('gammas', 'cis')] = compute_CI_gamma(results[('lambdas', 'cis')], mus)

    if save_results:
        print("Saving results...")

        np.save(os.path.join(location_output,
                             f"results_TP53_KRAS_model_{key}.npy"),
                results)

    return results


def compute_TP53_KRAS_gene_model(gene, key='pan_data', compute_CIs=False, save_results=True):
    """Compute a model for M=3 that includes TP53, KRAS and gene.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking

    :type save_results: bool
    :param save_results: If True (default) save results.

    :rtype: dict
    :return: A dictionary with all results (lambdas, mus and gammas).

    """
    results = {}

    print("Loading samples...")
    samples = pd.read_csv(samples_per_combination_files[key],
                          index_col='third gene')

    samples = samples.loc[gene]

    results['samples'] = convert_samples_to_dict(samples)
    print("...done.")


    print("Computing fluxes MLE...")
    bounds = np.repeat(1, 12)
    MLE = estimate_lambdas(samples, draws=1,
                           upper_bound_prior=bounds,
                           kwargs={'return_raw':True})[0]

    results[('lambdas', 'mle')] = convert_lambdas_to_dict(MLE)
    print("...done.")


    if compute_CIs:
        print("Computing fluxes CI...")
        results[('lambdas', 'cis')] = convert_lambdas_to_dict(
            asymp_CI_lambdas(MLE['lambdas'], samples))
        print("...done.")

    print("Loading mutation rates...")
    if key[-4:] == 'plus':
        final_part_key = 'w_panel'
    else:
        final_part_key = key[-4:]
    mus = pd.read_csv(
        os.path.join(location_output,
                     f"{key[:-4] + final_part_key}_mutation_rates.txt"),
        index_col='gene')

    if key == 'pan_data':
        mus = {(1, 0, 0): float(mus.loc['TP53']),
               (0, 1, 0): float(mus.loc['KRAS']),
               (0, 0, 1): float(mus.loc[gene])}
    else:
        mus = {(1, 0, 0): mus.loc['TP53', 'rate_grp_1'],
               (0, 1, 0): mus.loc['KRAS', 'rate_grp_1'],
               (0, 0, 1): mus.loc[gene, 'rate_grp_1']}
    results['mus'] = mus
    print("...done.")


    print("Computing selection coefficients...")
    results[('gammas', 'mle')] = compute_gammas(results[('lambdas', 'mle')], mus)
    if compute_CIs:
        results[('gammas', 'cis')] = compute_CI_gamma(results[('lambdas', 'cis')], mus)
    print("...done.")


    if save_results:
        print("Saving results...")

        np.save(os.path.join(location_output,
                             f"results_TP53_KRAS_{gene}_model_{key}.npy"),
                results)
        print("...done.")

    return results


def load_mutation_rates(key, method): 
    """Load the mutation rates in a format that :func:`compute_gammas` uses.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking
        - smoking_plus
        - nonsmoking_plus
    
    :type method: str or NoneType
    :param method: Which mutation rate calculation method to use. Can be one of:
        - variant (sum of mutation rate for each variant site in gene)
        - cesR (gene mutation rate directly from cancereffectsizeR)

    :rtype: dict
    :return: A dictionary with the genes as keys and the mutation rates as values.
        - gene: mutation rate
    """
    # TODO: Account for when mutation rates are not available for a gene
    
    if(method == "variant"):
        variant_based_mutation_rates = pd.read_csv(os.path.join(location_output,
                                                        'variant_based_mutation_rates.txt'),
                                                        index_col=0)
        genes_available = set.intersection(set(gene_list),
                                           set(variant_based_mutation_rates.index))
        variant_based_mutation_rates.columns = "pan_data","smoking","nonsmoking"
        variant_based_mutation_rates = variant_based_mutation_rates.to_dict()

        if key[-4:] == 'plus':
            key = key[:-5]

        mus = variant_based_mutation_rates[key]

    elif(method == "cesR"):
        if key[-4:] == 'plus':
            final_part_key = 'w_panel'
        else:
            final_part_key = key[-4:]
        mus_df = pd.read_csv(
            os.path.join(location_output,
                        f"{key[:-4] + final_part_key}_mutation_rates.txt"),
            index_col='gene')
        
        mus = mus_df["rate_grp_1"].to_dict()
        
    else:
        raise ValueError("The only options for gene mutation rates are variant-sum-based (variant) or directly from cancereffectsizeR (cesR)")

    
    return mus

def load_mutation_rates_for_TP53_KRAS_gene_model(key, method):
    """Load the mutation rates in a format that :func:`compute_gammas_for_TP53_KRAS_gene_model`
    uses, where we will assume an M=3 model including TP53, KRAS and a third gene.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking
        - smoking_plus
        - nonsmoking_plus
    
    :type method: str or NoneType
    :param method: Which mutation rate calculation method to use. Can be one of:
        - variant (sum of mutation rate for each variant site in gene)
        - cesR (gene mutation rate directly from cancereffectsizeR)

    :rtype: dict
    :return: A dictionary of dictionaries indexed by the third gene,
        where the mutation rates are represented by:
        - (1, 0, 0): TP53
        - (0, 1, 0): KRAS
        - (0, 0, 1): the third gene in the model.
    """
    if(method == "variant"):
        variant_based_mutation_rates = pd.read_csv(os.path.join(location_output,
                                                        'variant_based_mutation_rates.txt'),
                                                        index_col=0)
        genes_available = set.intersection(set(gene_list),
                                           set(variant_based_mutation_rates.index))
        variant_based_mutation_rates.columns = "pan_data","smoking","nonsmoking"
        variant_based_mutation_rates = variant_based_mutation_rates.to_dict()

        if key[-4:] == 'plus':
            key = key[:-5]

        variant_based_mutation_rates = variant_based_mutation_rates[key]

        mus = {gene:{(1, 0, 0): variant_based_mutation_rates['TP53'],
                     (0, 1, 0): variant_based_mutation_rates['KRAS'],
                     (0, 0, 1): variant_based_mutation_rates[gene]}
           for gene in genes_available
           if gene not in ['TP53', 'KRAS']}

    elif(method == "cesR"):
        if key[-4:] == 'plus':
            final_part_key = 'w_panel'
        else:
            final_part_key = key[-4:]
        mus_df = pd.read_csv(
            os.path.join(location_output,
                        f"{key[:-4] + final_part_key}_mutation_rates.txt"),
            index_col='gene')
        genes_available = set.intersection(set(gene_list),
                                        set(mus_df.index))
        mus = {gene:{(1, 0, 0): mus_df.loc['TP53', 'rate_grp_1'],
                    (0, 1, 0): mus_df.loc['KRAS', 'rate_grp_1'],
                    (0, 0, 1): mus_df.loc[gene, 'rate_grp_1']}
            for gene in genes_available
            if gene not in ['TP53', 'KRAS']}
        
    else:
        raise ValueError("The only options for gene mutation rates are variant-sum-based (variant) or directly from cancereffectsizeR (cesR)")

    
    return mus

def compute_all_gammas(key, all_lambdas, mus, save_results=True):
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
    :param mus: Dictionary with the mutation rates for each gene in
        data set `key` as obtained from :func:`load_mutation_rates`.

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

    for i, combo in enumerate(lambdas_mles.keys()):

        combo_length = len(combo)
        # Constructs dictionary of the form:
        # mu_combo = {combo:{(1, 0, 0): mus[combo[0]],
        #                    (0, 1, 0): mus[combo[1]],
        #                    (0, 0, 1): mus[combo[2]]}}
        # if combo_length == 3
        mu_combo = {tuple([0]*i + [1] + [0]*(combo_length-i-1)):
                        mus[combo[i]]
                     for i in range(0,combo_length)}

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


def compute_all_gammas_for_TP53_KRAS_gene_model(key, all_lambdas, mus, save_results=True):
    """Compute all estimates of the selection coefficient for the data
    set `key` iterating over all genes in the intersection of
    :const:`gene_list` and the genes available from the results of
    CES.

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
    :param mus: Dictionary with the results of mutation rates for
        `key` as obtained from :func:`load_mutation_rates`.

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

    genes_available = set.intersection(set(gene_list),
                                       set(lambdas_mles.keys()),
                                       set(mus.keys()))

    for i, gene in enumerate(genes_available):

        gammas_mles[gene] = compute_gammas(
            lambdas_mles[gene],
            mus[gene])

        gammas_cis[gene] = compute_CI_gamma(
            lambdas_cis[gene],
            mus[gene])

        if save_results:
            np.save(os.path.join(location_output,
                                 f"{key}_selections_mles.npy"),
                    gammas_mles)
            np.save(os.path.join(location_output,
                                 f"{key}_selections_cis.npy"),
                    gammas_cis)

    return gammas_mles, gammas_cis


def main(genes = gene_list, num_per_combo=3, flexible_last_layer=False, recompute_samples_per_combination=False, save_results=True):
    """Main method for the estimation of all the fluxes.

    :type genes: list or NoneType
    :param genes: List of genes from which gene-set combinations will
        be made.

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
    for key in results_keys:
        print("")
        if (recompute_samples_per_combination
            or not os.path.exists(samples_per_combination_files[key])):
            print(f"Computing number of samples per combination for {key}...")
            all_counts[key] = compute_samples_for_all_combinations(genes, key, num_per_combo, save_results)
        else:
            print(f"Loading counts per combination for {key}...")
            df = pd.read_csv(samples_per_combination_files[key], index_col='gene combination')
            all_counts[key] = {combo:np.array(df.loc[combo]) for combo in combinations(genes, num_per_combo)}
        print(f"done computing samples per combination for {key}.")
        print("")
        print("")

        print(f"Estimating all epistatic models for {key}...")
        print("")
        lambdas_mles, lambdas_cis = compute_all_lambdas(key, all_counts, flexible_last_layer, save_results)
        all_lambdas[(key, 'mles')] = lambdas_mles
        all_lambdas[(key, 'cis')] = lambdas_cis
        print(f"done estimating all epistatic models for {key}.")
        print("")
        print("")

        print(f"Computing selection coefficients for {key}...")
        print("")
        mus = load_mutation_rates(key, method = "variant")
        gammas_mles, gammas_cis = compute_all_gammas(key, all_lambdas, mus, save_results)
        all_gammas[(key, 'mles')] = gammas_mles
        all_gammas[(key, 'cis')] = gammas_cis
        print(f"done computing selection coefficients for {key}.")
        print("")
        print("")


    return all_counts, all_lambdas, all_gammas

def estimate_fluxes_TP53_KRAS_gene_model(recompute_samples_per_combination=False, save_results=True):
    """Main method for the estimation of all the fluxes.

    :type recompute_samples_per_combination: bool
    :param recompute_samples_per_combination: If True force
        recomputing the samples per combination of each
        gene. Otherwise (default) try load the respective file if
        available.

    :type save_results: bool
    :param save_results: If True (default) save results.

    :rtype: tuple
    :return: A tuple with two dictionaries that contain all the
        samples per combination and all the fluxes estimated. The keys
        of the dictionaries are the keys in :const:`results_keys`.

    """
    all_counts = {}
    all_lambdas = {}
    all_gammas = {}
    for key in results_keys:
        print("")
        if (recompute_samples_per_combination
            or not os.path.exists(samples_per_combination_files[key])):
            print(f"Computing number of samples per combination for {key}...")
            all_counts[key] = compute_samples_for_TP53_KRAS_gene_model(key, save_results)
        else:
            print(f"Loading counts per combination for {key}...")
            df = pd.read_csv(samples_per_combination_files[key], index_col='third gene')
            all_counts[key] = {gene:np.array(df.loc[gene]) for gene in gene_list[2:]}
        print(f"done computing samples per combination for {key}.")
        print("")
        print("")

        print(f"Estimating all epistatic models for {key}...")
        print("")
        lambdas_mles, lambdas_cis = compute_all_lambdas(key, all_counts, save_results)
        all_lambdas[(key, 'mles')] = lambdas_mles
        all_lambdas[(key, 'cis')] = lambdas_cis
        print(f"done estimating all epistatic models for {key}.")
        print("")
        print("")

        print(f"Computing selection coefficients for {key}...")
        print("")
        mus = load_mutation_rates(key, method = "variant")
        gammas_mles, gammas_cis = compute_all_gammas(key, all_lambdas, mus, save_results)
        all_gammas[(key, 'mles')] = gammas_mles
        all_gammas[(key, 'cis')] = gammas_cis
        print(f"done computing selection coefficients for {key}.")
        print("")
        print("")


    return all_counts, all_lambdas, all_gammas

if __name__ == "__main__":
    main(recompute_samples_per_combination = True)



# import matplotlib.pyplot as plt
# from locations import location_figures

# default_genes = ['TP53', 'KRAS', 'EGFR']

# df_smoking = pd.read_csv(samples_per_combination_files['smoking_plus'], index_col='third gene')
# counts_smoking = {gene:np.array(df_smoking.loc[gene]) for gene in gene_list[2:]}['EGFR']
# counts_smoking = {x:y for x, y in zip(build_S_with_tuples(3), counts_smoking)}

# df = pd.DataFrame(columns=default_genes)

# for mutation, count in counts_smoking.items():
#     temp_df = pd.DataFrame([mutation]*count*4, columns=df.columns)
#     df = pd.concat([temp_df, df])

# df = df.replace(0, 0.25)

# # Convert DataFrame to numpy array for plotting
# data_array = df.T.astype(float)

# # for col in data_array.columns:
# #     for _ in range(2):
# #         data_array = pd.concat([data_array, data_array[col]], axis=1)

# total_smoking = sum(counts_smoking.values())

# for i in range(4, total_smoking*5, 5):
#     data_array.insert(i, f'Z{i}', np.zeros(len(data_array)))


# # Plot the heatmap
# plt.figure(figsize=(10, 1.8))
# plt.imshow(data_array, cmap='gray_r', interpolation='none', aspect='auto')
# # plt.colorbar(label='Mutation status')
# percs = np.round(np.array([sum([value for x, value in counts_smoking.items()
#                                 if x[i] == 1])
#                            for i in range(3)])/total_smoking*100).astype(int)
# percs = ["47", "36",  " 8"]
# plt.yticks(ticks=np.arange(3), labels=[f"$\it{y}$    {perc}%"
#                                        for y, perc in zip(df.columns,
#                                                           percs)])#
# # np.round(np.array([sum([value for x, value in counts_smoking.items() if x[i] == 1]) for i in range(3)])/total_smoking*100)

# plt.xticks([])
# # plt.xlabel('Patients')

# # plt.ylabel('Genes')
# plt.title(f"Lung adenocarcinoma smokers (multiple sources) $n={total_smoking}$")
# plt.savefig(os.path.join(location_figures,
#                          "heatmap_smoking.png"), dpi=600)



# df_nonsmoking = pd.read_csv(samples_per_combination_files['nonsmoking_plus'], index_col='third gene')
# counts_nonsmoking = {gene:np.array(df_nonsmoking.loc[gene]) for gene in gene_list[2:]}['EGFR']
# counts_nonsmoking = {x:y for x, y in zip(build_S_with_tuples(3), counts_nonsmoking)}


# df = pd.DataFrame(columns=default_genes)

# for mutation, count in counts_nonsmoking.items():
#     temp_df = pd.DataFrame([mutation]*count*4, columns=df.columns)
#     df = pd.concat([temp_df, df])

# df = df.replace(0, 0.25)

# # Convert DataFrame to numpy array for plotting
# data_array = df.T.astype(float)

# # for col in data_array.columns:
# #     for _ in range(2):
# #         data_array = pd.concat([data_array, data_array[col]], axis=1)

# total_nonsmoking = sum(counts_nonsmoking.values())

# for i in range(4, total_nonsmoking*5, 5):
#     data_array.insert(i, f'Z{i}', np.zeros(len(data_array)))


# # Plot the heatmap
# plt.figure(figsize=(10, 1.8))
# plt.imshow(data_array, cmap='gray_r', interpolation='none', aspect='auto')
# # plt.colorbar(label='Mutation status')

# percs = np.round(np.array([sum([value for x, value in counts_nonsmoking.items()
#                                 if x[i] == 1])
#                            for i in range(3)])/total_nonsmoking*100).astype(int)
# percs = ["37", " 8", "26"]
# plt.yticks(ticks=np.arange(3), labels=[f"$\it{y}$    {perc}%"
#                                        for y, perc in zip(df.columns,
#                                                           percs)])#

# plt.xticks([])
# # plt.xlabel('Patients')

# # plt.ylabel('Genes')
# plt.title(f"Lung adenocarcinoma non-smokers (multiple sources) $n={total_nonsmoking}$")
# plt.savefig(os.path.join(location_figures,
#                          "heatmap_nonsmoking.png"), dpi=600)
