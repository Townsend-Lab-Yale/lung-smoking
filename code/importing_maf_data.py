import os
import pandas as pd
from importing_clinical_data import dup_sample_ids_1718, multi_sample_ids_2017, keep_tracer_samples, genie_non_luad_id, dup_sample_ids_gen18, multi_sample_ids_genie, dup_sample_ids_gen17, metastatic_sample_ids_2017, metastatic_sample_ids_tracer, metastatic_sample_ids_genie

if '__file__' not in globals():
    __file__ = '.'

location_data = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "data/"))
"""Location of directory that contains data for the model."""

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
        db_file_name = os.path.join(location_data, db_file_name)

    return db_file_name


def filter_db_by_mutation(db=default_db_file_name,
                          tumor_col_name=None,
                          patient_id_col_name=None,
                          clear_silent=True):
    """Build a data frame with patient IDs and mutations.

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

    :type patient_id_col_name: str
    :param patient_id_col_name: Name of the column in the data that
        contains the patient identification. If None, try to infer it
        from the data.

    :type clear_silent: bool
    :param clear_silent: Remove silent mutations from the data. Leave
        as True (default), to later provide ranges for a gene.

    :rtype: pandas.core.frame.DataFrame
    :return: Pandas DataFrame with patient IDs and mutations.

    """

    if isinstance(db, list):
        db = [db_full_name(x) for x in db]
        dfs = [pd.read_csv(x, sep="\t", comment="#")
               for x in db]
        data = pd.concat(dfs, ignore_index=True)
    else:
        db = db_full_name(db)
        data = pd.read_csv(db, sep="\t", comment="#")

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

    if patient_id_col_name is None:
        if 'Tumor_Sample_Barcode' in data.columns:
            patient_id_col_name = 'Tumor_Sample_Barcode'
        elif 'Unique_Patient_Identifier' in data.columns:
            patient_id_col_name = 'Unique_Patient_Identifier'
        else:
            raise Exception("Unknown patient identifier. "
                            "Provide variable 'patient_id_col_name'")

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

    data = data[[patient_id_col_name, 'Chromosome', 'Start_Position', 'Mutation', 'Variant_Classification']]
    data.columns = ['Sample ID', 'Chromosome', 'Start_Position', 'Mutation', 'Variant_Classification']

    return data

files1 = ["tsp.luad.maf.txt", "oncosg.luad.maf.txt","mskcc.2015.luad.maf.txt", "broad.luad.maf.txt", "mskcc.2018.nsclc.maf.txt"]
result1 = filter_db_by_mutation(db = files1)
#removing metastatic sample from broad dataset
result1 = result1[result1['Sample ID'] != 'LU-A08-43']

#has to be done separately so samples from MSK 2018 with same sample IDs aren't accidentally removed
result2 = filter_db_by_mutation(db = "mskcc.2017.luad.maf.txt")
result2 = result2[~result2['Sample ID'].isin(metastatic_sample_ids_2017)]
#multi sample ids refers to sample IDs that should be removed as they are extra samples from the same patient
result2 = result2[~result2['Sample ID'].isin(multi_sample_ids_2017)]
#dup sample ids refers to IDs repeated between MSK 2017 and 2018 and then removed from 2017.
result2 = result2[~result2['Sample ID'].isin(dup_sample_ids_1718)]

result3 = filter_db_by_mutation(db = 'tcga.luad.maf.txt', patient_id_col_name="case_id")

result4 = filter_db_by_mutation(db = "tracerx.nsclc.maf.txt")
result4 = result4[~result4['Sample ID'].isin(metastatic_sample_ids_tracer)]
result4 = result4[result4['Sample ID'].isin(keep_tracer_samples)]

result5 = filter_db_by_mutation(db = 'genie.maf.txt')
#print(len(pd.unique(result5['Sample ID'])))
#output: 98465 sample IDs in maf file vs 110704 sample IDs in the clinical file
result5 = result5[~result5['Sample ID'].isin(genie_non_luad_id)]
#print(len(pd.unique(result5['Sample ID'])))
#output: 11774 non-LUAD sample IDs in maf file vs 12926 non-LUAD sample IDs in the clinical file
result5 = result5[~result5['Sample ID'].isin(metastatic_sample_ids_genie)]
result5 = result5[~result5['Sample ID'].isin(multi_sample_ids_genie)]
result5 = result5[~result5['Sample ID'].isin(dup_sample_ids_gen18)]
result5 = result5[~result5['Sample ID'].isin(dup_sample_ids_gen17)]

final = pd.concat([result1, result2, result3, result4, result5])
final.to_csv('output/merged_luad_maf.txt')
