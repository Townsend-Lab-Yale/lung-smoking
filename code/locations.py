import os

if '__file__' not in globals():
    __file__ = '.'

location_figures = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "figures/"))
"""Location of directory that contains figures for the project."""


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
    'FM-AD':'luad_fm-ad'}
"""Location of the directories that contain each data set in the
`location_data' directory.

"""


full_maf_file_names = {
    db:os.path.join(location_data,
                    directory,
                    "data_mutations_extended.txt")
    for db, directory in data_sets_directories.items()}


full_maf_file_names_lifted = {
    db:os.path.join(location_data,
                    directory,
                    "data_mutations_extended_lifted.txt")
    for db, directory in data_sets_directories.items()}


gene_list_file = os.path.join(location_data, "genes_list.txt")
gene_coordinates_file = os.path.join(location_data, "gene_coordinates.csv")
mutation_rates_file = os.path.join(location_data, "mutation_rates.txt")


merged_maf_file_name = os.path.join(location_output, 'merged_luad_maf.txt')
merged_clinical_file_name = os.path.join(location_output, 'merged_luad_clinical.txt')
merged_maf_clinical_file_name = os.path.join(location_output, 'merged_final.txt')
