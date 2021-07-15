import numpy as np
import pandas as pd

def compute_samples(data,
                    mutations=None,
                    tumor_col_name=None,
                    print_info=True,
                    patient_id_col_name=None):
    """Compute samples numbers for each mutation combinations in S.

    :type data: pandas.core.frame.DataFrame
    :param data: Data frame with the mutation data.

    :type mutations: list or int or NoneType
    :param mutations: List with the mutation names to extract
        information about. This list determines the total number of
        mutations, M, and thus the space of mutation combinations,
        S. The order in the list is also important as the vectors in S
        correspond to that order. Each item on the list can be a
        string, or a list of strings. In the latter case (NOT
        IMPLEMENTED YET), those mutations will be aggregated into a
        single category. If instead of a list an int M is passed, pick
        the M most common mutations from the database, and if None
        (default), pick M to be equal to :const:`default_M`.

    :type tumor_col_name: str
    :param tumor_col_name: Name of the column in the data that
        contains the tumor allele. If None, try to infer it from the
        data.

    :type patient_id_col_name: str
    :param patient_id_col_name: Name of the column in the data that
        contains the patient identification. If None, try to infer it
        from the data.

    :rtype: tuple
    :return: A tuple with a one dimensional array of size 2^M with the
        computed sample numbers per mutation combination (ordered as
        items in S), M as in integer and S as a numpy array.

    """

    if mutations is None:
        mutations = default_M

    if not isinstance(mutations, list):
        mutations = list(most_common_mutations(data, mutations).index)

    M = len(mutations)

    S = build_S_as_array(M)

    if print_info:
        print("Mutation labels")
        for k, mutation in enumerate(mutations):
            print(f"  {np.eye(1, M, k, dtype=np.int)[0]} = {mutation}")

    pts_per_mutation = [set(data[
        data['Start_Position'] >= genes.loc[mutation, 'start']][
            data['Start_Position'] <= genes.loc[mutation, 'end']][
                data['Chromosome'] == str(genes.loc[mutation, 'chromosome'])][
                    'Patient ID']) if mutation in genes.index else
        set(data[data['Mutation'] == mutation]['Patient ID'])
                        for mutation in mutations]

    pts_per_combination = (
        [len(data['Patient ID'].unique())]
        + [len(set.intersection(*[pts_per_mutation[i]
                                  for i in indices]))
           for indices in [[i for i in range(M)
                            if Sj[i] == 1]
                           for Sj in S[1:]]])

    pts_per_combination = np.array(pts_per_combination)

    for m in range(M, 0, -1):
        for k, y in enumerate(S):
            if np.sum(y) == m:
                in_y_indices = np.array(
                    [True
                     if np.sum(y >= S[i]) == M
                     and np.sum(y - S[i]) >= 1
                     else False
                     for i in range(len(S))])
                pts_per_combination[in_y_indices] = (
                    pts_per_combination[in_y_indices]
                    - pts_per_combination[k])

    if print_info:
        print("")
        print("Counts per combination")
        for x, pts in zip(S, pts_per_combination):
            print(f"  {x} : {pts}")

    return pts_per_combination
