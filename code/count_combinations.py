import numpy as np
import pandas as pd

from theory import build_S_as_array

from locations import gene_coordinates_file

genes = pd.read_csv(gene_coordinates_file, index_col='gene')

def updated_compute_samples(data,
                    mutations,
                    print_info=False):
    """Compute samples numbers for each mutation combinations in S.

    :type data: pandas.core.frame.DataFrame
    :param data: Table in which each row is a sample and the
        columns are genes. The value of each sample/gene combination
        is a binary representation of whether the sample has a mutation
        in that gene or not.

    :type mutations: list
    :param mutations: List with the mutation names to extract
        information about. This list determines the total number of
        mutations, M, and thus the space of mutation combinations,
        S. The order in the list is also important as the vectors in S
        correspond to that order. Each item on the list can be a
        string, or a list of strings. In the latter case (NOT
        IMPLEMENTED YET), those mutations will be aggregated into a
        single category.

    """

    pts_per_combination = data.groupby(mutations)['Sample ID'].count()
    levels = [[0,1]] * len(mutations)
    new_index = pd.MultiIndex.from_product(levels, names=pts_per_combination.index.names)
    pts_per_combination = pts_per_combination.reindex(new_index, fill_value=0)

    if print_info:
        print(pts_per_combination.reset_index().rename(columns={'Sample ID':'Sample Count'}))
    return np.array(pts_per_combination)


def convert_samples_to_dict(samples):
    """Convert a samples array to a dictionary.

    :type samples: list
    :param samples: Number of patients in each mutation combination
        (as returned by :func:`count_pts_per_combination`).

    :rtype: dict
    :return: Dictionary with the samples, indexed by tuples of 1's and
        0's representing whether the mutation occur fot the gene or
        not.

    """

    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    S = build_S_as_array(M)

    results_as_dict = {tuple(x):value
                       for x, value in zip(S, samples)}

    return results_as_dict
