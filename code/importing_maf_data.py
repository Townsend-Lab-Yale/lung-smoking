import os
import pandas as pd
from importing_clinical_data import dup_sample_ids_1718, multi_sample_ids_2017, keep_tracer_samples, genie_non_luad_id, dup_sample_ids_gen18, multi_sample_ids_genie, dup_sample_ids_gen17, metastatic_sample_ids_2017, metastatic_sample_ids_tracer, metastatic_sample_ids_genie, fmad_non_luad, fmad_non_primary

from locations import full_maf_file_names
from locations import location_output


data_sets_sample_id_col_names ={key:'case_id' if key in ['TCGA', 'FM-AD']
                                else None
                                for key in full_maf_file_names.keys()}



def import_maf_data(db,
                    tumor_col_name=None,
                    sample_id_col_name=None,
                    clear_silent=True):
    """Build a data frame with relevant information including sample IDs
and mutations.

    :typed db: str or list
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
            -'FM-AD'

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

    if isinstance(db, list):
        dfs = [pd.read_csv(full_maf_file_names[x],
                           sep="\t",
                           comment="#")
               for x in db]
        dfs = [df.assign(Source=x)
               for df, x in zip(dfs, db)]
        data = pd.concat(dfs, ignore_index=True)
    else:
        data = pd.read_csv(full_maf_file_names[db],
                           sep="\t",
                           comment="#")
        data = data.assign(Source=db)

    if clear_silent:
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

    # removes 'chr' when there:
    data['Chromosome'] = data.apply(
        lambda x: str(x['Chromosome']).split("chr")[-1], axis=1)
    # includes new column with unique description of mutation
    data['Mutation'] = data.apply(
        lambda x: "{}:{} {}>{}".format(
            x['Chromosome'],
            x['Start_Position'],
            x['Reference_Allele'],
            x[tumor_col_name]),
        axis=1)

    cols_except_id =  (['Chromosome', 'Start_Position', 'Mutation', 'Source']
                       + (['Variant_Classification'] if not clear_silent else []))
    data = data[[sample_id_col_name] + cols_except_id]
    data.columns = ['Sample ID'] + cols_except_id

    return data


dfs = {
    db:import_maf_data(
        db,
        sample_id_col_name=data_sets_sample_id_col_names[db])
    for db in full_maf_file_names.keys()}


dfs['Broad'] = dfs['Broad'][dfs['Broad']['Sample ID'] != 'LU-A08-43']

dfs['MSK2017'] = dfs['MSK2017'][~dfs['MSK2017']['Sample ID'].isin(metastatic_sample_ids_2017)]
dfs['MSK2017'] = dfs['MSK2017'][~dfs['MSK2017']['Sample ID'].isin(multi_sample_ids_2017)]
dfs['MSK2017'] = dfs['MSK2017'][~dfs['MSK2017']['Sample ID'].isin(dup_sample_ids_1718)]

dfs['TracerX'] = dfs['TracerX'][~dfs['TracerX']['Sample ID'].isin(metastatic_sample_ids_tracer)]
dfs['TracerX'] = dfs['TracerX'][dfs['TracerX']['Sample ID'].isin(keep_tracer_samples)]

dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(genie_non_luad_id)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(metastatic_sample_ids_genie)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(multi_sample_ids_genie)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(dup_sample_ids_gen18)]
dfs['Genie'] = dfs['Genie'][~dfs['Genie']['Sample ID'].isin(dup_sample_ids_gen17)]

dfs['FM-AD'] = dfs['FM-AD'][~dfs['FM-AD']['Sample ID'].isin(fmad_non_luad)]
dfs['FM-AD'] = dfs['FM-AD'][~dfs['FM-AD']['Sample ID'].isin(fmad_non_primary)]


final_df = pd.concat([df for df in dfs.values()])
final_df.to_csv(os.path.join(location_output, 'merged_luad_maf.txt'))
