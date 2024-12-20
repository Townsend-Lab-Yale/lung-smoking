import os

from string import ascii_uppercase as alphabet

if '__file__' not in globals():
    __file__ = '.'

location_figures = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "figures/"))
"""Location of directory that contains figures for the project."""

if not os.path.exists(location_figures):
    os.mkdir(location_figures)

location_data = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "data/"))
"""Location of directory that contains data for the model."""


location_output = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "output/"))
"""Location of directory that contains output for the model."""

if not os.path.exists(location_output):
    os.mkdir(location_output)

data_sets_directories = {
    'TSP':'luad_tsp',
    'OncoSG':'luad_oncosg_2020',
    'MSK2015':'luad_mskcc_2015',
    'Broad':'luad_broad',
    'MSK2018':'nsclc_pd1_msk_2018',
    'MSK2017':'lung_msk_2017',
    'TCGA':'luad_tcga',
    'TracerX':'nsclc_tracerx_2017',
    'Genie':'genie_9',
    'FM-AD':'luad_fm-ad',
    'NCI':'lung_nci_2022',
    'CPTAC':'luad_cptac_2020',
    'Yale':'yale_luad'}
"""Location of the directories that contain each data set in the
`location_data' directory.

"""


location_gene_panels = os.path.join(location_data,
                                    "gene_panels")
"""Location of directory that contains the gene panels."""

all_panel_genes_file_name = os.path.join(location_gene_panels,
                                         'all_panel_genes.txt')

all_panel_samples_file_name = os.path.join(location_gene_panels,
                                           'all_panel_samples.txt')



full_maf_file_names = {
    db:os.path.join(location_data,
                    directory,
                    "data_mutations_extended.txt")
    for db, directory in data_sets_directories.items()}

gene_list_file = os.path.join(location_data, "genes_list.txt")


pts_by_mutation_file = os.path.join(location_output, "pts_by_mutation.csv")


merged_maf_file_name = os.path.join(location_output, 'merged_luad_maf.txt')
merged_clinical_file_name = os.path.join(location_output, 'merged_luad_clinical.txt')
merged_maf_clinical_file_name = os.path.join(location_output, 'merged_final.txt')

cesR_filtered_maf_file_name = os.path.join(location_output, 'cesR_maf_for_epistasis_analysis.txt')

genes_per_sample_file_name = os.path.join(location_output, 'genes_per_sample.txt')

results_keys = ["pan_data", "smoking", "nonsmoking", "smoking_plus", "nonsmoking_plus"]
"""List with possible keys to classify the data.

"""

full_mutation_rate_file_names = {
    key:os.path.join(location_output,
                    f"{key}_mutation_rates.txt")
    for key in results_keys[:3]
}

full_flux_mle_file_names = {
    key:os.path.join(location_output,
                    f"{key}_fluxes_mles.npy")
    for key in results_keys
}

full_flux_ci_file_names = {
    key:os.path.join(location_output,
                    f"{key}_fluxes_cis.npy")
    for key in results_keys
}

full_selection_mle_file_names = {
    key:os.path.join(location_output,
                    f"{key}_selections_mles.npy")
    for key in results_keys
}

full_selection_ci_file_names = {
    key:os.path.join(location_output,
                    f"{key}_selections_cis.npy")
    for key in results_keys
}

smoking_sample_ids_file = os.path.join(location_data,
                                       'smoking_sample_ids.txt')
nonsmoking_sample_ids_file = os.path.join(location_data,
                                          'nonsmoking_sample_ids.txt')
panel_smoking_sample_ids_file = os.path.join(location_data,
                                       'panel_smoking_sample_ids.txt')
panel_nonsmoking_sample_ids_file = os.path.join(location_data,
                                          'panel_nonsmoking_sample_ids.txt')


default_mutation_names = list(alphabet)
"""Default mutation names, in case `mutation_names` is not provided
for the plotting functions. The default is A, B, C, etc.

"""
