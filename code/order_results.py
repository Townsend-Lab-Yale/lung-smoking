## Ordering of results

def find_model(genes, main_gene_list):
    """If results where produced from `main_gene_list' then for any
    non-ordered tuple or list `genes', return the rightly ordered
    tuple that will index the results.

    For example, if main_gene_list = ['TP53', 'KRAS', 'EGFR', 'BRAF', ...]

    and ('TP53', 'BRAF', 'EGFR') is provided as genes

    it will return

    ('TP53', 'EGFR', 'BRAF')

    """

    gene_positions = {gene: index for index, gene in enumerate(main_gene_list)}

    # Sort the genes based on their positions in the main_gene_list
    sorted_genes = sorted(genes, key=lambda gene: gene_positions[gene])

    # Convert the sorted list back to a tuple and return
    return tuple(sorted_genes)



def order_genes_by_result_values(results):
    """Give a list of the genes ordered by values.

    :type results: dict
    :param results: Dictionary with the results for which will be
        ordered by result values. It should be indexed by tuples of
        genes that were included in the model, and where here we will
        only take the M=1 (no epistasis) results for the ordering. Use
        for example as input: load_results('selections')['pan_data']
        to order by selection, where load_results is the function from
        load_results.py

    :rtype: list
    :return: A list with the genes, ordered by values in `results'.

    """

    gene_list = {gene[0]:max(value.values()) # max doesn't
                                             # matter there
                                             # is only one
                                             # value because
                                             # M = 1
                 for gene, value in results.items()
                 if len(gene) == 1}
    gene_list = sorted(gene_list,
                       key=lambda k: gene_list[k],
                       reverse=True)

    return gene_list
