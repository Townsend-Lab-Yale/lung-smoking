import numpy as np
from theory import epistatic_comparisons


def epistatic_ratios(results, M, results_cis=None):

    """Compute ratios between values of the `results' in a somatic
    genotype vs another one with one less mutated gene.

    :type results: dict
    :param results: Dictionary with the results for which we will
        compute epistatic ratios. It should be indexed by tuples
        of genes that were included in the model. Here we can put
        either selections (gammas), fluxes (lambdas, but Jeff doesn't
        like this approach), or later even mutation rates (mus).

    :type M: int
    :param M: Number of genes to consider in the models

    :type results_cis: dict or None
    :param results_cis: Dictionary with the results confidence
       intervals. Results of ratios are set to 1 if they are not
       statistically significant. Statistical significance is
       conservatively determined by non overlapping confidence
       intervals. If None is provided (default), then do consider
       statistical significance for the ratios.


    :rtype: dict
    :return: A dictionary with the epistatic ratios, with keys being
       the genes in each model.
    """

    comparisons = epistatic_comparisons(M)

    if results_cis is None:
        ratios = {genes:{comparison:(
            results[genes][comparison[1]]/results[genes][comparison[0]])
                         for comparison in comparisons}
                  for genes in results.keys() if len(genes) == M}
    else:
        ratios = {genes:{comparison:(
            results[genes][
                comparison[1]]/results[genes][comparison[0]])
                         if (max(results_cis[genes][comparison[0]][0],
                                 results_cis[genes][comparison[1]][0]) >
                             min(results_cis[genes][comparison[0]][1],
                                 results_cis[genes][comparison[1]][1]))
                         else 1
                         for comparison in comparisons}
                  for genes in results.keys() if len(genes) == M}

    return ratios


def epistatic_ratios_2_matrix(results, genes_ordered_list, results_cis=None):
    """Order results of the pair-wise epistatic ratios in a matrix.

    This function calls :func:`epistatic_ratios` for M=2 to obtain all
    epistatic ratios: the ratio of `results' (say selection) for a
    gene in a somatic genotype with another gene mutated over
    `results' for the same gene in a normal somatic genotype (the
    other gene not mutated).

    The results are provided as len(genes) x len(genes) matrix where
    the rows represent the common gene (say, the gene selected), and
    the columns the other gene for which the somatic genotype context
    that ratio compares (with or without that other gene).

    :type results: dict
    :param results: Dictionary with the results for which we will
        compute epistatic ratios. It should be indexed by tuples
        of genes that were included in the model. Here we can put
        either selections (gammas), fluxes (lambdas, but Jeff doesn't
        like this approach), or later even mutation rates (mus).

    :type M: int
    :param M: Number of genes to consider in the models

    :type results_cis: dict or None
    :param results_cis: Dictionary with the results confidence
       intervals. Results of ratios are set to 1 if they are not
       statistically significant. Statistical significance is
       conservatively determined by non overlapping confidence
       intervals. If None is provided (default), then do consider
       statistical significance for the ratios.


    :rtype: dict
    :return: A dictionary with the epistatic ratios, with keys being
       the genes in each model.

    """

    ratios_dict = epistatic_ratios(results, 2, results_cis)

    n = len(genes_ordered_list)

    matrix = np.array(n*[n*[np.nan]])

    for genes, values in ratios_dict.items():
        for comparison, value in values.items():
            gene_selected = genes[comparison[0][1].index(1)]
            context_gene = genes[comparison[1][0].index(1)]
            matrix[genes_ordered_list.index(gene_selected)][
                genes_ordered_list.index(context_gene)] = value

    return matrix



def epistatic_ratios_3rd_gene_effects(results, results_cis):
    """When does the presence of a third gene in a somatic genotype,
    alter the epistatic effect of the selection of a gene in the
    presence of another.

    Results are a dictionary where the keys are of the form

        ('gene_always_previously_mutated_in_comparison', 'gene_selected', 'third_gene')

    and the values are the gene ratios of

        gene_always_there+third_gene -> gene_always_there+third_gene+gene_selected

    and

       gene_always_there -> gene_always_there+gene_selected

    """

    ratios_dict = epistatic_ratios(results, 3, results_cis)

    effects = {}

    for genes, values in ratios_dict.items():
        for comparison, value in values.items():
            if value != 1:
                from_gene = genes[(comparison[0][0]).index(1)]
                to_gene_tuple = tuple(np.array(comparison[0][1]) -
                                      np.array(comparison[0][0]))
                to_gene = genes[to_gene_tuple.index(1)]
                gene_3rd_tuple = tuple(np.array(comparison[1][0]) -
                                       np.array(comparison[0][0]))
                gene_3rd = genes[gene_3rd_tuple.index(1)]
                effects[(from_gene, to_gene, gene_3rd)] = value


    return effects
