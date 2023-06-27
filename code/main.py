import os
import pandas as pd
import numpy as np

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

from filter_data import dbs_filtered_for_TP53_KRAS
from filter_data import filter_samples_for_genes


gene_list = list(pd.read_csv(gene_list_file, header=None)[0])
gene_list = [gene.upper() for gene in gene_list]
gene_list = gene_list[:103]


def compute_samples_for_all_genes(key=None, save_results=True):
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

    number_of_genes = len(genes)

    unrepresented_genes = []

    counts = {}

    for i, gene in enumerate(genes):
        if (i+1) % 100 == 0 or i == 0:
            print(f"Computing gene {i+1}/{number_of_genes} "
                  f"({round((i+1)/number_of_genes, 2)*100}% done)")
        db = filter_samples_for_genes(gene, dbs_filtered_for_TP53_KRAS[key])

        if gene in db.columns:
            counts[gene] = updated_compute_samples(db,
                                        mutations=['TP53', 'KRAS', gene])
        else:
            unrepresented_genes.append(gene)
            counts[gene] = np.repeat(0, 8)
    if len(unrepresented_genes) > 0:
        print(f"No sample availability information for the following genes: {str(unrepresented_genes)}. This may be because there it is not mutated in any tumors or because it is not correctly represented in sequencing data")

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


def compute_all_lambdas(key, all_counts, save_results=True):
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


def main(recompute_samples_per_combination=False, save_results=True):
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
    for key in results_keys:
        print("")
        if (recompute_samples_per_combination
            or not os.path.exists(samples_per_combination_files[key])):
            print(f"Computing number of samples per combination for {key}...")
            all_counts[key] = compute_samples_for_all_genes(key, save_results)
        else:
            f"Loading counts per combination for {key}..."
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

    return all_counts, all_lambdas

if __name__ == "__main__":
    main(recompute_samples_per_combination = True)
