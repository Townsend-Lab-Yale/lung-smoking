
import numpy as np
import pandas as pd
from itertools import product, combinations
from time import time
from ranked_list import produce_ranked_list

from locations import gene_coordinates_file
from locations import merged_luad_file_name

genes = pd.read_csv(gene_coordinates_file, index_col='gene')

def build_S_with_tuples(M):
    """Build entire space of mutations combinations S.

    It will represented by a list of size 2^M with items being tuples
    of size M with 0s and 1s.

    :type M: int
    :param M: Number of mutations.

    :rtype: list
    :return: S as an ordered list of tuples.

    """
    return [x for x in product([0, 1], repeat=M)]

def build_S_as_array(M):
    """Build entire space of mutations combinations S.

    It will represented by a numpy array of shape (2^M, M).

    :type M: int
    :param M: Number of mutations.

    :rtype: numpy.ndarray
    :return: S as a numpy array of 0s and 1s.

    """
    return np.array(build_S_with_tuples(M))

def number_mutations(data, mutations = None):
    numbers_per_mutation = {}

    numbers_per_mutation = {mutation: len(data[
            (data['Start_Position'] >= genes.loc[mutation, 'start']) &
                (data['Start_Position'] <= genes.loc[mutation, 'end']) &
                (data['Chromosome'] == genes.loc[mutation, 'chromosome'])]
            ['Sample ID'].unique())
                        for mutation in mutations}

    result = [(k,v) for k, v in sorted(numbers_per_mutation.items(), key=lambda item: item[1])][::-1]

    with open('ranked_mutations.txt', 'w') as fp:
        fp.write('\n'.join('%s %s' % x for x in result))
        fp.write('\n')

    return True

def compute_samples(data,
                    mutations=None,
                    tumor_col_name=None,
                    print_info=False,
                    patient_id_col_name=None,
                    N=10,
                    R=6):
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

    :rtype: tuple
    :return: A tuple with a one dimensional array of size 2^M with the
        computed sample numbers per mutation combination (ordered as
        items in S), M as in integer and S as a numpy array.

    """
    temp_file = open('mut_sets.txt','w')
    temp_file.close()

    if mutations == None:

        results = []
        working_mutation_sets = []

        ranked_list = produce_ranked_list(set_length_=R)
        combination_mutation_sets = []

        for indices in combinations(range(N),R):
            combination_mutation_sets.append([ranked_list[ind] for ind in indices])

        for mutations in combination_mutation_sets:
            M = len(mutations)

            S = build_S_as_array(M)

            if print_info:
                print("Mutation labels")
                for k, mutation in enumerate(mutations):
                    print(f"  {np.eye(1, M, k, dtype=np.int)[0]} = {mutation}")

            '''creates list of sets containing sample IDs in which a mutation is present, index of set corresponds to index of gene in mutation list'''

            pts_per_mutation = [set(data[
                data['Start_Position'] >= genes.loc[mutation, 'start']][
                    data['Start_Position'] <= genes.loc[mutation, 'end']][
                        data['Chromosome'] == str(genes.loc[mutation, 'chromosome'])][
                            'Sample ID'])
                                for mutation in mutations]

            pts_per_combination = (

                [len(data['Sample ID'].unique())] #for no mutations
                + [len(set.intersection(*[pts_per_mutation[i] #performs set intersection across all sets containing a mutation in the gene currently selected in the for loop looping across genes mutated in row Sj
                #this is where I need to check if len = 0
                                        for i in indices]))
                for indices in [[i for i in range(M) #for genes that are mutated in that row (Sj)
                                    if Sj[i] == 1]
                                for Sj in S[1:]]])# if (len(set.intersection(*[pts_per_mutation[i] for i in indices])) > 0)]) #for row in possible combinations of mutations

            '''
                for indices in [[i for i in range(M)
                                    if Sj[i] == 1]
                                for Sj in S[1:]]:
                    if len(set.intersection(*[pts_per_mutation[i] for i in indices])) > 0:
                        [len(data['Sample ID'].unique())] + len(set.intersection(*pts_per_mutation[i])))
            '''

            pts_per_combination = np.array(pts_per_combination)

            #if 0 in pts_per_combination:
            #    continue

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
            #print(pts_per_combination)
            if 0 not in pts_per_combination[:-1]:
                results.append(pts_per_combination)
                working_mutation_sets.append(mutations)
                with open ('mut_sets.txt','a') as out_file:
                    out_file.write(str(mutations) + ', ')
                print(pts_per_combination)
        #return results
        print(results)
        print(working_mutation_sets)

    else:
        temp1234 = open('combinations.txt','w')
        temp1234.close()

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
            #this is where I need to check if len = 0
                                    for i in indices]))
            for indices in [[i for i in range(M) #for genes that are mutated in that row (Sj)
                                if Sj[i] == 1]
                            for Sj in S[1:]]])# if (len(set.intersection(*[pts_per_mutation[i] for i in indices])) > 0)]) #for row in possible combinations of mutations

        '''
            for indices in [[i for i in range(M)
                                if Sj[i] == 1]
                            for Sj in S[1:]]:
                if len(set.intersection(*[pts_per_mutation[i] for i in indices])) > 0:
                    [len(data['Sample ID'].unique())] + len(set.intersection(*pts_per_mutation[i])))
        '''

        pts_per_combination = np.array(pts_per_combination)

        #if 0 in pts_per_combination:
        #    continue

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
        if 0 not in pts_per_combination[:-1]:
            with open ('mut_sets.txt','a') as out_file:
                out_file.write(str(mutations))
                out_file.write("\n")
            with open ('combinations.txt','a') as out_file2:
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
        ['Our Sample ID'])
                        for mutation in mutations]
    in_all = len(set.intersection(*pts_per_mutation))
    all_but_one = [len(set.intersection(
        *(pts_per_mutation[0:i] + pts_per_mutation[i+1:]))) - in_all
                   for i in range(M)]
    if 0 not in all_but_one:
        return 0 not in compute_samples(data, mutations)[:-1]
    else:
        return False



maf_file = pd.read_csv(merged_luad_file_name)
maf_file = maf_file[maf_file['Variant_Classification'] != 'Silent']
all_pts_ordered = list(maf_file["Sample ID"].unique())
#maf_file["Our Sample ID"] = maf_file["Sample ID"].apply( #this apply function takes 8-9 seconds to run
#    lambda x: all_pts_ordered.index(x))



#number_mutations(maf_file, mutations = genes.index)
t3 = time()
print(len(genes.index))
if number_mutations(maf_file, mutations=genes.index):
    print('Successful')
print(time() - t3)
'''
possible = []
for five_genes in twenty_choose_five:
    for added_gene in ranked_list_test5[35:45]:
        if combination_works(maf_file,combination = five_genes + [added_gene]):
            possible.append(five_genes + [added_gene])
print(possible)
'''

'''
possible_sets = open('possible_sets.txt','w')
possible_sets.close()

#possible = []
N = 30
R = 6

ranked_list = produce_ranked_list(set_length_=R)
combination_mutation_sets = []

for indices in combinations(range(N),R):
    combination_mutation_sets.append([ranked_list[ind] for ind in indices])
for mut_set in combination_mutation_sets:
    if combination_works(maf_file, combination = mut_set):
        with open('possible_sets.txt','a') as possible_sets:
            possible_sets.write(str(mut_set))
        #possible.append(mut_set)
#with open('possible_sets.txt','a') as possible_sets:
#    possible_sets.write(str(possible))

t4 = time()
with open('timer.txt','w') as time_file:
    time_file.write(str(t4-t3))

possible_sets_all5 = open('possible_sets_all5.txt','w')
possible_sets_all5.close()

R_2 = 5
ranked_list_2 = produce_ranked_list(set_length_=R_2)
N_2 = len(ranked_list_2)

combination_mutation_sets_2 = []
for indices in combinations(range(N_2),R_2):
    combination_mutation_sets_2.append([ranked_list_2[ind] for ind in indices])
for mut_set in combination_mutation_sets_2:
    if combination_works(maf_file, combination = mut_set):
        with open('possible_sets_all5.txt','a') as possible_sets_all5:
            possible_sets_all5.write(str(mut_set))

t5 = time()
with open('timer.txt','a') as time_file:
    time_file.write(str(t5-t4))
'''