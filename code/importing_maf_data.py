import os
import pandas as pd
from importing_clinical_data import dup_sample_ids_1718, multi_sample_ids_2017, keep_tracer_samples, genie_non_luad_id, dup_sample_ids_gen18, multi_sample_ids_genie, dup_sample_ids_gen17, metastatic_sample_ids_2017, metastatic_sample_ids_tracer, metastatic_sample_ids_genie, fmad_non_luad, fmad_non_primary

if '__file__' not in globals():
    __file__ = '.'

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


default_db_file_name = os.path.join(location_data, "tcga.luad.maf.txt")
"""Default file name to use as data source."""


def db_full_name(db_file):
    """Return full path file name of `db_file'.

    :typed db_file: str
    :param db_file: Name of the MAF file containing the data. It
        should be located in the :const:`location_data`. Data files
        names available thus far are:

            - 'cca_ihc.txt' alias 'ihc'
            - 'melanoma.maf.txt' alias 'melanoma'
            - 'rectal_adenocarcinoma.maf.txt' alias 'rectal'
            - 'prostate_adenocarcinoma.maf.txt' alias 'prostate'
            - 'ucec.maf.txt' alias 'ucec'
            - 'tcga.luad.maf.txt' alias 'luad'
            - 'tcga.lusc.maf.txt' alias 'lusc'

    """
    db_file_name = db_file.lower()
    '''
    if db_file_name in "ihc":
        db_file_name = "cca_ihc.txt"
    if db_file_name == "melanoma":
        db_file_name = 'melanoma.maf.txt'
    if db_file_name == "rectal":
        db_file_name = 'rectal_adenocarcinoma.maf.txt'
    if db_file_name == "prostate":
        db_file_name = 'prostate_adenocarcinoma.maf.txt'
    if db_file_name == "ucec":
        db_file_name = 'ucec.maf.txt'
    if db_file_name == 'luad':
        db_file_name = 'tcga.luad.maf.txt'
    if db_file_name == 'lusc':
        db_file_name = 'tcga.lusc.maf.txt'

    if "/" not in db_file_name:
    '''
    db_file_name = os.path.join(location_data, db_file_name)

    return db_file_name


def filter_db_by_mutation(db=default_db_file_name,
                          tumor_col_name=None,
                          sample_id_col_name=None,
                          clear_silent=True):
    """Build a data frame with sample IDs and mutations.

    :type db: str or list
    :param db: Name of the MAF file (or files if a db is a list)
        containing the data. It should be located in the
        :const:`location_data`. Data files names available thus far
        are:

            - 'cca_ihc.txt' alias 'ihc'
            - 'melanoma.maf.txt' alias 'melanoma'
            - 'rectal_adenocarcinoma.maf.txt' alias 'rectal'
            - 'prostate_adenocarcinoma.maf.txt' alias 'prostate'
            - 'ucec.maf.txt' alias 'ucec'
            - 'tcga.luad.maf.txt' alias 'luad'
            - 'tcga.lusc.maf.txt' alias 'lusc'

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
        db_full = [db_full_name(x) for x in db]
        dfs = [pd.read_csv(x, sep="\t", comment="#")
               for x in db_full]
        dfs = [df.assign(Source=db_name)
               for df, db_name in zip(dfs, db)]
        data = pd.concat(dfs, ignore_index=True)
    else:
        db_full = db_full_name(db)
        data = pd.read_csv(db_full, sep="\t", comment="#")
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

files = {
    'TSP':['luad_tsp/data_mutations_extended.txt',''],
    'OncoSG':['luad_oncosg_2020/data_mutations_extended.txt',''],
    'MSK2015':['luad_mskcc_2015/data_mutations_extended.txt',''],
    'Broad':['luad_broad/data_mutations_extended.txt',''],
    'MSK2018':['nsclc_pd1_msk_2018/data_mutations_extended.txt',''],
    'MSK2017':['lung_msk_2017/data_mutations_extended.txt',''],
    'TCGA':['luad_tcga/data_mutations_extended.txt','case_id'],
    'TracerX':['nsclc_tracerx_2017/data_mutations_extended.txt',''],
    'Genie':['genie_9/data_mutations_extended.txt',''],
    'FM-AD':['luad_FM-AD/data_mutations_extended.txt','case_id']
}

for key, value in files.items():
    if value[1] == '':
        files[key].append(filter_db_by_mutation(db = value[0]))
    else:
        files[key].append(filter_db_by_mutation(db = value[0], sample_id_col_name=value[1]))

files['Broad'][2] = files['Broad'][2][files['Broad'][2]['Sample ID'] != 'LU-A08-43']

files['MSK2017'][2] = files['MSK2017'][2][~files['MSK2017'][2]['Sample ID'].isin(metastatic_sample_ids_2017)]
files['MSK2017'][2] = files['MSK2017'][2][~files['MSK2017'][2]['Sample ID'].isin(multi_sample_ids_2017)]
files['MSK2017'][2] = files['MSK2017'][2][~files['MSK2017'][2]['Sample ID'].isin(dup_sample_ids_1718)]

files['TracerX'][2] = files['TracerX'][2][~files['TracerX'][2]['Sample ID'].isin(metastatic_sample_ids_tracer)]
files['TracerX'][2] = files['TracerX'][2][files['TracerX'][2]['Sample ID'].isin(keep_tracer_samples)]

files['Genie'][2] = files['Genie'][2][~files['Genie'][2]['Sample ID'].isin(genie_non_luad_id)]
files['Genie'][2] = files['Genie'][2][~files['Genie'][2]['Sample ID'].isin(metastatic_sample_ids_genie)]
files['Genie'][2] = files['Genie'][2][~files['Genie'][2]['Sample ID'].isin(multi_sample_ids_genie)]
files['Genie'][2] = files['Genie'][2][~files['Genie'][2]['Sample ID'].isin(dup_sample_ids_gen18)]
files['Genie'][2] = files['Genie'][2][~files['Genie'][2]['Sample ID'].isin(dup_sample_ids_gen17)]

files['FM-AD'][2] = files['FM-AD'][2][~files['FM-AD'][2]['Sample ID'].isin(fmad_non_luad)]
files['FM-AD'][2] = files['FM-AD'][2][~files['FM-AD'][2]['Sample ID'].isin(fmad_non_primary)]

final_file = pd.concat([value[2] for value in files.values()])
final_file.to_csv(os.path.join(location_output, 'merged_luad_maf.txt'))
