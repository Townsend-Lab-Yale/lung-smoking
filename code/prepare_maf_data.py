import os
import pandas as pd

from filter_clinical_data import dup_sample_ids_1718
from filter_clinical_data import multi_sample_ids_2017
from filter_clinical_data import keep_tracer_samples
from filter_clinical_data import dup_sample_ids_gen18
from filter_clinical_data import multi_sample_ids_genie
from filter_clinical_data import dup_sample_ids_gen17
from filter_clinical_data import metastatic_sample_ids_2017
from filter_clinical_data import metastatic_sample_ids_tracer
from filter_clinical_data import metastatic_sample_ids_genie
from filter_clinical_data import fmad_non_primary
from filter_clinical_data import non_luad_sample_ids_fmad
from filter_clinical_data import non_luad_sample_ids_msk2015
from filter_clinical_data import non_luad_sample_ids_msk2018
from filter_clinical_data import non_luad_sample_ids_tracer
from filter_clinical_data import non_luad_sample_ids_genie
from filter_clinical_data import non_luad_sample_ids_nci
from filter_clinical_data import non_yale_sample_ids



from locations import full_maf_file_names as maf_file_names
from locations import location_gene_panels
from locations import merged_maf_file_name


data_sets_sample_id_col_names ={key:'case_id' if key in ['TCGA', 'FM-AD']
                                else 'tumor_name' if key == 'Yale'
                                else None
                                for key in maf_file_names.keys()}



def import_maf_data(db,
                    tumor_col_name=None,
                    sample_id_col_name=None,
                    clear_silent=True):
    """Build a data frame with relevant information including sample IDs
    and mutations.

    :type db: str or list
    :param db: Name of the data set or list containing multiple data
        set names. Data files names available thus far are:

            -'TSP'
            -'OncoSG'
            -'MSK2015'
            -'Broad'
            -'MSK2018'
            -'MSK2017'
            -'TCGA'
            -'TracerX'
            -'Genie'
            -'NCI'
            -'FM-AD'
            -'Yale'

    :type tumor_col_name: str
    :param tumor_col_name: Name of the column in the data that
        contains the tumor allele. If None, try to infer it from the
        data.

    :type sample_id_col_name: str
    :param sample_id_col_name: Name of the column in the data that
        contains the sample identification. If None, try to infer it
        from the data.

    :type clear_silent: bool
    :param clear_silent: Remove silent mutations from the data. Leave
        as True (default), to later provide ranges for a gene.

    :rtype: pandas.core.frame.DataFrame
    :return: Pandas DataFrame with sample IDs and mutations.

    """

    dfs = []
    if isinstance(db, list):
        for x in db:
            file_path = maf_file_names[x]
            if x == 'FM-AD' and not os.path.exists(file_path):
                # File does not exist, create an empty DataFrame
                print(f"File for '{x}' not found. Creating an empty DataFrame.")
                df = pd.DataFrame(columns=[
                    'Sample ID', 'Chromosome', 'Start_Position', 'Mutation',
                    'Reference_Allele', 'Tumor_Seq_Allele2', 'Source', 'Variant_Classification'
                ])
            else:
                # File exists, read the data
                df = pd.read_csv(file_path, sep="\t", comment="#")
                df = df.assign(Source=x)

            dfs.append(df)

        # Concatenate all DataFrames
        data = pd.concat(dfs, ignore_index=True, sort=False)

    else:
        file_path = maf_file_names[db]
        if db == 'FM-AD' and not os.path.exists(file_path):
            # File does not exist, create an empty DataFrame
            print(f"File for '{db}' not found. Creating an empty DataFrame.")

            data = pd.DataFrame(columns=[
                'Sample ID', 'Chromosome', 'Start_Position', 'Mutation',
                'Reference_Allele', 'Tumor_Seq_Allele2', 'Source',
                'Variant_Classification'])
        else:
            # File exists, read the data
            data = pd.read_csv(file_path, sep="\t", comment="#")
            data = data.assign(Source=db)

    if not data.empty:
        if clear_silent and 'Variant_Classification' in data.columns:
            data = data[data['Variant_Classification'] != 'Silent']

        if tumor_col_name is None:
            if 'Tumor_Seq_Allele2' in data.columns:
                tumor_col_name = 'Tumor_Seq_Allele2'
            elif 'Tumor_Allele' in data.columns:
                tumor_col_name = 'Tumor_Allele'
            else:
                raise Exception("Unknown tumor allele. "
                                "Provide variable 'tumor_col_name'")

        if sample_id_col_name is None:
            if 'Tumor_Sample_Barcode' in data.columns:
                sample_id_col_name = 'Tumor_Sample_Barcode'
            elif 'Unique_Patient_Identifier' in data.columns:
                sample_id_col_name = 'Unique_Patient_Identifier'
            else:
                raise Exception("Unknown sample identifier. "
                                "Provide variable 'sample_id_col_name'")

        # Remove 'chr' prefix if present
        if 'Chromosome' in data.columns:
            data['Chromosome'] = data['Chromosome'].astype(str).str.replace('chr', '', regex=False)

        # Create a unique mutation description
        data['Mutation'] = data.apply(
            lambda x: "{}:{} {}>{}".format(
                x['Chromosome'],
                x['Start_Position'],
                x['Reference_Allele'],
                x[tumor_col_name]),
            axis=1)

        # Prepare the list of columns to keep
        cols_except_id = ['Chromosome', 'Start_Position', 'Mutation',
                          'Reference_Allele', tumor_col_name, 'Source']
        if not clear_silent and 'Variant_Classification' in data.columns:
            cols_except_id.append('Variant_Classification')

        # Select and rename columns
        data = data[[sample_id_col_name] + cols_except_id]
        data.columns = ['Sample ID'] + cols_except_id

    else:
        # Handle the case where data is empty
        if tumor_col_name is None:
            tumor_col_name = 'Tumor_Seq_Allele2'  # Default column name
        if sample_id_col_name is None:
            sample_id_col_name = 'Tumor_Sample_Barcode'  # Default column name

        # Create empty DataFrame with expected columns
        cols_except_id = ['Chromosome', 'Start_Position', 'Mutation',
                          'Reference_Allele', tumor_col_name, 'Source']
        if not clear_silent:
            cols_except_id.append('Variant_Classification')

        data = pd.DataFrame(columns=['Sample ID'] + cols_except_id)

    return data


dfs = {
    db:import_maf_data(
        db,
        sample_id_col_name=data_sets_sample_id_col_names[db],
        clear_silent=False)
    for db in maf_file_names.keys()}


## Filter data sets
dfs['Broad'] = dfs['Broad'][dfs['Broad']['Sample ID'] != 'LU-A08-43']

dfs['OncoSG'] = dfs['OncoSG'][~(dfs['OncoSG']['Chromosome'].isin(['GL000230.1', 'hs37d5', 'GL000211.1','MT', 'GL000192.1', 'GL000214.1', 'GL000241.1', 'GL000220.1', 'GL000212.1','GL000205.1', 'GL000195.1', 'GL000218.1', 'GL000216.1', 'GL000226.1','GL000224.1', 'GL000231.1', 'GL000221.1', 'GL000234.1', 'GL000219.1','GL000191.1', 'GL000229.1', 'GL000238.1']))]

dfs['MSK2015'] = dfs['MSK2015'][~dfs['MSK2015']['Sample ID'].isin(non_luad_sample_ids_msk2015)]

dfs['MSK2018'] = dfs['MSK2018'][~dfs['MSK2018']['Sample ID'].isin(non_luad_sample_ids_msk2018)]

dfs['MSK2017'] = dfs['MSK2017'][~dfs['MSK2017']['Sample ID'].isin(metastatic_sample_ids_2017)]
dfs['MSK2017'] = dfs['MSK2017'][~dfs['MSK2017']['Sample ID'].isin(multi_sample_ids_2017)]
dfs['MSK2017'] = dfs['MSK2017'][~dfs['MSK2017']['Sample ID'].isin(dup_sample_ids_1718)]

dfs['TracerX'] = dfs['TracerX'][~dfs['TracerX']['Sample ID'].isin(metastatic_sample_ids_tracer)]
dfs['TracerX'] = dfs['TracerX'][~dfs['TracerX']['Sample ID'].isin(non_luad_sample_ids_tracer)]
dfs['TracerX'] = dfs['TracerX'][dfs['TracerX']['Sample ID'].isin(keep_tracer_samples)]

dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(non_luad_sample_ids_genie)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(metastatic_sample_ids_genie)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(multi_sample_ids_genie)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(dup_sample_ids_gen18)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(dup_sample_ids_gen17)]

dfs['NCI'] = dfs['NCI'][~dfs['NCI']['Sample ID'].isin(non_luad_sample_ids_nci)]

dfs['FM-AD'] = dfs['FM-AD'][~dfs['FM-AD']['Sample ID'].isin(non_luad_sample_ids_fmad)]
dfs['FM-AD'] = dfs['FM-AD'][~dfs['FM-AD']['Sample ID'].isin(fmad_non_primary)]

dfs['Yale'] = dfs['Yale'][~dfs['Yale']['Sample ID'].isin(non_yale_sample_ids)]

## Merge data

final_df = pd.concat([df for df in dfs.values()])
final_df = final_df.reset_index(drop=True)


# ## Add panel information

panels_used = {key:pd.read_csv(
    os.path.join(location_gene_panels,
                 f"{key.lower()}_panels_used.txt")).set_index(
                     "Sample Identifier").to_dict()
    for key in ['Genie', 'MSK2017', 'MSK2018']}

for key in ['Genie', 'MSK2017', 'MSK2018']:
    final_df.loc[final_df['Source'] == key, 'Panel'] = final_df['Sample ID'].map(
        panels_used[key]['Sequence Assay ID' if key == 'Genie'
                         else 'Gene Panel'])

final_df.loc[final_df['Source'] == 'TSP', 'Panel'] = 'TSP'
final_df.loc[final_df['Source'] == 'FM-AD', 'Panel'] = 'FoundationOne'


## Save merged data to file
final_df.to_csv(merged_maf_file_name)
