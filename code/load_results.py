import os
import pandas as pd
import numpy as np


from locations import location_output
from locations import results_keys


def load_mutation_rates(key, method, extension=None):
    """Load the mutation rates in a format that :func:`compute_gammas` uses.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking
        - smoking_plus
        - nonsmoking_plus

    :type method: str
    :param method: Which mutation rate calculation method to use. Can be one of:
        - variant (sum of mutation rate for each variant site in gene)
        - cesR (gene mutation rate directly from cancereffectsizeR)

    :type extension: str
    :param extension: Which mutation rate calculation method to use. Can be one of:
        - variant (sum of mutation rate for each variant site in gene)
        - cesR (gene mutation rate directly from cancereffectsizeR)

    :rtype: dict
    :return: A dictionary with the genes as keys and the mutation rates as values.
        - gene: mutation rate
    """

    if extension is not None:
        extended_location_output = os.path.join(location_output, extension)
    else:
        extended_location_output = location_output

    if key[-4:] == 'plus':
        accession_key = key[:-5]
    else:
        accession_key = key

    if(method == "variant"):
        variant_based_mutation_rates = pd.read_csv(os.path.join(extended_location_output,
                                                                'variant_based_mutation_rates.txt'),
                                                   index_col=0)
        variant_based_mutation_rates = variant_based_mutation_rates.to_dict()

        mus = variant_based_mutation_rates[accession_key]
        mus = {gene: mu for gene, mu in mus.items() if not np.isnan(mu)}

    elif(method == "cesR"):
        mus_df = pd.read_csv(
            os.path.join(extended_location_output,
                        f"{accession_key}_mutation_rates.txt"),
            index_col='gene')

        mus = mus_df["rate_grp_1"].to_dict()

    else:
        raise ValueError("The only options for gene mutation rates "
                         "are variant-sum-based (variant) or directly "
                         "from cancereffectsizeR (cesR)")


    return mus



def load_results(result_type, which=None, extension=None):
    """Load the results.

    :type result_type: str
    :param result_type: What type of result to load. Can be one of:
        - 'samples'
        - 'mutations'
        - 'fluxes'
        - 'selections'

    :type which: str or NoneType
    :param which: What estimates to use.
        - if `result_type` is 'fluxes' or 'selections, it can be one of:
            + 'mles' (default), to load the MLE
            + 'cis', to load the 95% confidence intervals
        - if `result_type` is 'mutations', it can be one of:
            + 'variant' (default), to load the MLE
            + 'cesR', to load the 95% confidence intervals
        - if `result_type` is 'samples' it does not have an effect


    :rtype: dict
    :return: A dictionary the fluxes results indexed by the keys as in
        :const:`results_keys`.

    """

    if result_type == 'mutations':
        if which is None:
            which = 'variant'
        results = {key: load_mutation_rates(key, which, extension)
                   for key in results_keys}
    else:
        if which is None:
            which = 'mles'

        file_names = {key: f"{key}_{result_type}" +
                           (f"_{which}" if result_type != 'samples'
                            else "") +
                           ".npy"
                      for key in results_keys}

        if extension is not None:
            extended_location_output = os.path.join(location_output,extension)
        else:
            extended_location_output = location_output

        results = {key: np.load(os.path.join(extended_location_output,
                                             file_name),
                                allow_pickle=True).item()
                   for key,file_name in file_names.items()}

    return results
