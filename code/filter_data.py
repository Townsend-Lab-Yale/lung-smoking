"""filter_data -- Jorge A. Alfaro-Murillo

This module filters the merged data to remove panels that do not
include the genes of interest (the genes that would be used for the
epistasis model). This is a necessary step because otherwise the count
for the "normal" genotype (that without any mutations in any of the
three genes) would be overestimated as we do not know if in panels
without a gene the gene is mutated or not.

"""

import pandas as pd


from locations import results_keys
from locations import all_panel_genes_file_name
from locations import all_panel_samples_file_name
from locations import nonsmoking_sample_ids_file
from locations import smoking_sample_ids_file
from locations import panel_nonsmoking_sample_ids_file
from locations import panel_smoking_sample_ids_file
from locations import genes_per_sample_file_name

all_panel_genes = pd.read_csv(all_panel_genes_file_name)
all_panel_samples = pd.read_csv(all_panel_samples_file_name)
all_panels = pd.unique(all_panel_genes['SEQ_ASSAY_ID'])


smoking_sample_ids = pd.read_csv(smoking_sample_ids_file,
                                 header=None)[0].to_list()
nonsmoking_sample_ids = pd.read_csv(nonsmoking_sample_ids_file,
                                    header=None)[0].to_list()
panel_smoking_sample_ids = pd.read_csv(panel_smoking_sample_ids_file,
                                 header=None)[0].to_list()
panel_nonsmoking_sample_ids = pd.read_csv(panel_nonsmoking_sample_ids_file,
                                    header=None)[0].to_list()


def filter_samples_for_genes(genes, db, print_info=False):
    """Remove for the database `db` all patients from panels that do not
    include the `genes`.

    Return the filtered database.

    """
    if type(genes) == str:
        genes = [genes]

    panels_to_remove = [
        panel for panel in all_panels
        if any(gene not in all_panel_genes[
            all_panel_genes['SEQ_ASSAY_ID'] == panel]['Hugo_Symbol'].tolist()
            for gene in genes)]

    samples_to_remove = all_panel_samples[all_panel_samples['Panel'].isin(panels_to_remove)]['Sample ID']
    if print_info:
        print("Panels excluded because they did not sequence at least one of "
              + ' or '.join(genes) + ': ' + ' '.join(panels_to_remove))
    return db[~db['Sample ID'].isin(samples_to_remove)]


def filter_db_for_key(key, db):
    """Filter the samples of the database `db` to the result `key`:

        - If the `key` is 'pan_data' do not filter.
        - If the `key` is 'smoking' filter to patients with the
          smoking signature.
        - If the `key` is 'nonsmoking' filter to patients without the
          smoking signature.

    """
    if key == 'pan_data':
        return db[db['Sample ID'].isin(smoking_sample_ids + panel_smoking_sample_ids + nonsmoking_sample_ids + panel_nonsmoking_sample_ids)]
    elif key == 'smoking':
        return db[db['Sample ID'].isin(smoking_sample_ids)]
    elif key == 'nonsmoking':
        return db[db['Sample ID'].isin(nonsmoking_sample_ids)]
    elif key == 'smoking_plus':
        return db[db['Sample ID'].isin(smoking_sample_ids + panel_smoking_sample_ids)]
    elif key == 'nonsmoking_plus':
        return db[db['Sample ID'].isin(nonsmoking_sample_ids + panel_nonsmoking_sample_ids)]

genes_per_sample = pd.read_csv(genes_per_sample_file_name)

key_filtered_dbs = {key: filter_db_for_key(key, genes_per_sample)
                   for key in results_keys}
