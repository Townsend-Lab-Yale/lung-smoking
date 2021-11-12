import os
import pandas as pd
import numpy as np
from count_combinations import compute_samples

from cancer_epistasis import estimate_lambdas
from cancer_epistasis import asymp_CI_lambdas
from cancer_epistasis import convert_lambdas_to_dict

from theory import build_S_with_tuples

from locations import gene_list_file
from locations import location_output
from locations import results_keys
from locations import samples_per_combination_files

from filter_data import prefiltered_dbs
from filter_data import filter_db_for_gene

# from figures import plot_lambdas_gammas


gene_list = list(pd.read_csv(gene_list_file, header=None)[0])


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
    :param save_results: If True (default) save results.

    """

    if key is None:
        key = "pan_data"

    genes = gene_list[2:]

    number_of_genes = len(genes)

    counts = {}

    for i, gene in enumerate(genes):
        if (i+1) % 100 == 0 or i == 0:
            print(f"Computing gene {i+1}/{number_of_genes} "
                  f"({round((i+1)/number_of_genes, 2)*100}% done)")
        subsetted_db = filter_db_for_gene(gene, prefiltered_dbs[key])
        counts[gene] = compute_samples(subsetted_db,
                                       mutations=['TP53', 'KRAS', gene])

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
    mod:`count_combinations` module, but it this one is way faster if
    the samples per combination have already been computed.

    """
    return np.all(samples[:-1] > 0)


def lambdas_from_samples(samples):
    bounds = 1
    print(f"Bounds for fluxes: {bounds}")
    MLE = estimate_lambdas(samples, draws=bounds,
                           upper_bound_prior=1,
                           kwargs={'return_raw':True})
    print(f"MLE: {MLE[0]['lambdas']}")

    bound_changes = 0
    while MLE[1].fun < -1e+20:
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
    """Compute all estimates of the fluxes for the data set `key`.

    """
    lambdas_mles = {}
    lambdas_cis = {}

    counts = all_counts[key]

    for i, gene in enumerate(gene_list[2:]):
        samples = counts[gene]
        if not are_all_fluxes_computable(samples):
            print(f"Skipping estimation for gene {gene} "
                  "because not all fluxes are computable for that model.")
        else:
            print(f"Running model with third gene {gene} "
                  f"(gene number {i+1}/{len(gene_list[2:])})")

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

        print("")

    return lambdas_mles, lambdas_cis


def main(recompute_samples_per_combination=False, save_results=True):
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
    main()
