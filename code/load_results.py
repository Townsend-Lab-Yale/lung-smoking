import os
import pandas as pd
import numpy as np


from locations import location_output
from locations import results_keys
from locations import samples_per_combination_files


def load_mutation_rates(key, method):
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

    :rtype: dict
    :return: A dictionary with the genes as keys and the mutation rates as values.
        - gene: mutation rate
    """

    if(method == "variant"):
        variant_based_mutation_rates = pd.read_csv(os.path.join(location_output,
                                                        'variant_based_mutation_rates.txt'),
                                                        index_col=0)
        ## TODO: Fix the availability of pan_data (and rerun main.R)
        # variant_based_mutation_rates.columns = "pan_data","smoking","nonsmoking" pan_data is not available
        variant_based_mutation_rates.columns = "smoking", "nonsmoking"
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
        raise ValueError("The only options for gene mutation rates "
                         "are variant-sum-based (variant) or directly "
                         "from cancereffectsizeR (cesR)")


    return mus



def load_results(result_type, which=None):
    """Load the results.

    :type result_type: str
    :param result_type: What type of result to load. Can be one of:
        - 'samples'
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
        results = {key: load_mutation_rates(key, which)
                   for key in results_keys}
    else:
        if which is None:
            which = 'mles'

        file_names = {key: f"{key}_{result_type}" +
                           (f"_{which}" if result_type != 'samples'
                            else "") +
                           ".npy"
                      for key in results_keys}

        results = {key: np.load(os.path.join(location_output,
                                             file_name),
                                allow_pickle=True).item()
                   for key,file_name in file_names.items()}

    return results
