import pandas as pd
from locations import gene_coordinates_file
from locations import merged_maf_file_name
from locations import pts_by_mutation_file

genes = pd.read_csv(gene_coordinates_file, index_col='gene')

db = pd.read_csv(merged_maf_file_name)
db = db[db['Variant_Classification'] != 'Silent']


def compute_pts_per_mutation(mutations=None):
    numbers_per_mutation = {}

    if mutations is None:
        mutations = list(genes.index)

    numbers_per_mutation = {mutation: len(db[
        (db['Start_Position'] >= genes.loc[mutation, 'start']) &
        (db['Start_Position'] <= genes.loc[mutation, 'end']) &
        (db['Chromosome'] == genes.loc[mutation, 'chromosome'])]
                                          ['Sample ID'].unique())
                        for mutation in mutations}

    pts_per_mutation = {gene:number for gene, number in sorted(
        numbers_per_mutation.items(),
        key=lambda item: item[1], reverse=True)}

    return pts_per_mutation


def main(save_result=True):
    pts_per_mutation = compute_pts_per_mutation()
    if save_result:
        pd.Series(pts_per_mutation).to_csv(pts_by_mutation_file,
                                           index_label='gene',
                                           header=['number'])

    return pts_per_mutation
