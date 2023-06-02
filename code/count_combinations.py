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


def compute_samples(data,
                    mutations=None,
                    print_info=False,
                    save_info=False):
    """Compute samples numbers for each mutation combinations in S.

    :type data: pandas.core.frame.DataFrame
    :param data: Data frame with the mutation data.

    :type mutations: list
    :param mutations: List with the mutation names to extract
        information about. This list determines the total number of
        mutations, M, and thus the space of mutation combinations,
        S. The order in the list is also important as the vectors in S
        correspond to that order. Each item on the list can be a
        string, or a list of strings. In the latter case (NOT
        IMPLEMENTED YET), those mutations will be aggregated into a
        single category.

    :type tumor_col_name: str
    :param tumor_col_name: Name of the column in the data that
        contains the tumor allele. If None, try to infer it from the
        data.

    :type patient_id_col_name: str
    :param patient_id_col_name: Name of the column in the data that
        contains the patient identification. If None, try to infer it
        from the data.

    :rtype: numpy.ndarray
    :return: A one dimensional array of size 2^M with the computed
        sample numbers per mutation combination (ordered as items in
        an S given as the result of function func:`build_S_as_array`),
        where M is the number of `mutations`.

    """
    if save_info:
        temp_file = open('mut_sets.txt','w')
        temp_file.close()

    M = len(mutations)

    S = build_S_as_array(M)

    if print_info:
        print("Mutation labels")
        for k, mutation in enumerate(mutations):
            print(f"  {np.eye(1, M, k, dtype=np.int)[0]} = {mutation}")

    '''creates list of sets containing sample IDs in which a mutation is present, index of set corresponds to index of gene in mutation list'''

    pts_per_mutation = [set(
        data[
            (data['Start_Position'] >= genes.loc[mutation, 'start']) &
            (data['Start_Position'] <= genes.loc[mutation, 'end']) &
            (data['Chromosome'] == genes.loc[mutation, 'chromosome'])]
        ['Sample ID'])
                        for mutation in mutations]

    pts_per_combination = (

        [len(data['Sample ID'].unique())] #for no mutations
        + [len(set.intersection(*[pts_per_mutation[i] #performs set intersection across all sets containing a mutation in the gene currently selected in the for loop looping across genes mutated in row Sj
                                for i in indices]))
        for indices in [[i for i in range(M) #for genes that are mutated in that row (Sj)
                            if Sj[i] == 1]
                        for Sj in S[1:]]])#for row in possible combinations of mutations

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
                '''at this point, (0,0,1,1,1) would contain a subset of patients in (0,0,0,1,1); this removes those.'''
                pts_per_combination[in_y_indices] = (
                    pts_per_combination[in_y_indices]
                    - pts_per_combination[k])

    if print_info:
        print("")
        print("Counts per combination")
        for x, pts in zip(S, pts_per_combination):
            print(f"  {x} : {pts}")
    if save_info:
        if 0 not in pts_per_combination[:-1]:
            with open ('mut_sets.txt','a') as out_file:
                out_file.write(str(mutations))
                out_file.write("\n")
            with open ('combinations.txt','w') as out_file2:
                out_file2.write(str(pts_per_combination))
                out_file2.write("\n")
    return pts_per_combination


def are_all_fluxes_computable(data, mutations):
    """Return True if we can compute all fluxes, that is, if there are
    patients in `data` with every possible combination of `genes`
    (except for the combination of all genes).

    """
    M = len(mutations)
    S = build_S_as_array(M)
    pts_per_mutation = [set(
        data[
            (data['Start_Position'] >= genes.loc[mutation, 'start']) &
            (data['Start_Position'] <= genes.loc[mutation, 'end']) &
            (data['Chromosome'] == str(genes.loc[mutation, 'chromosome']))]
        ['Sample ID'])
                        for mutation in mutations]
    in_all = len(set.intersection(*pts_per_mutation))
    all_but_one = [len(set.intersection(
        *(pts_per_mutation[0:i] + pts_per_mutation[i+1:]))) - in_all
                   for i in range(M)]
    if 0 not in all_but_one:
        return 0 not in compute_samples(data, mutations)[:-1]
    else:
        return False


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
