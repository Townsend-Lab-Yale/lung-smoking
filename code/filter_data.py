"""filter_data -- Jorge A. Alfaro-Murillo

This module filters the merged data to remove panels that do not
include the genes of interest (the genes that would be used for the
epistasis model). This is a necessary step because otherwise the count
for the "normal" genotype (that without any mutations in any of the
three genes) would be overestimated as we do not know if in panels
without a gene the gene is mutated or not.

"""

import os

import pandas as pd


from locations import results_keys
from locations import all_panel_genes_file_name
from locations import all_panel_samples_file_name
from locations import nonsmoking_sample_ids_file
from locations import smoking_sample_ids_file
from locations import panel_nonsmoking_sample_ids_file
from locations import panel_smoking_sample_ids_file
from locations import genes_per_sample_file_name
from locations import indels_per_sample_file_name
from locations import location_output

all_panel_genes = pd.read_csv(all_panel_genes_file_name)
all_panel_samples = pd.read_csv(all_panel_samples_file_name)
all_panels = pd.unique(all_panel_genes['SEQ_ASSAY_ID'])

sample_id_override_dir = os.path.join(location_output, "sample_id_overrides")
"""Directory with optional task-specific cohort overrides.

Files are keyed by entries in `results_keys`. For each key, the following
plain-text files are recognized when present:
- `<key>.txt`: replace the default sample IDs for that key
- `<key>__add.txt`: add sample IDs after replacement/default selection
- `<key>__drop.txt`: remove sample IDs after replacement/default selection

Assumption:
- Each file contains one sample ID per line and is created before the modeling
  process starts. This keeps the override mechanism reproducible under
  multiprocessing because every process reads the same files from disk.
"""

sample_ids_by_key = {}
"""Resolved sample IDs used to define each analysis cohort."""

active_sample_id_overrides = {}
"""Metadata describing which task-specific overrides were applied."""

key_filtered_dbs = {}
"""Per-key filtered versions of `genes_per_sample` built from current files."""

key_filtered_indel_dbs = {}
"""Per-key filtered versions of `indels_per_sample` when present."""

genes_per_sample = None
"""Merged genotype table used as the substrate for cohort filtering."""

indels_per_sample = None
"""Merged frameshift-indel table aligned to `genes_per_sample` when present."""


def _ordered_unique(values):
    """Return values with duplicates removed while preserving first appearance."""
    return list(dict.fromkeys(values))


def _read_sample_ids(file_name):
    """Read one-sample-per-line identifiers from `file_name`.

    Inputs:
    - file_name: path to a plain-text sample ID file

    Outputs:
    - list[str]: sample IDs with order preserved and duplicates removed

    Assumptions:
    - Missing files are handled by the caller.
    - Sample IDs do not contain commas or embedded separators.
    """
    return _ordered_unique(
        pd.read_csv(file_name, header=None)[0].astype(str).str.strip().tolist()
    )


def _build_default_sample_ids_by_key():
    """Construct the default sample-ID definitions for each result key."""
    smoking_sample_ids = _read_sample_ids(smoking_sample_ids_file)
    nonsmoking_sample_ids = _read_sample_ids(nonsmoking_sample_ids_file)
    panel_smoking_sample_ids = _read_sample_ids(panel_smoking_sample_ids_file)
    panel_nonsmoking_sample_ids = _read_sample_ids(panel_nonsmoking_sample_ids_file)

    return {
        'pan_data': _ordered_unique(smoking_sample_ids
                                    + panel_smoking_sample_ids
                                    + nonsmoking_sample_ids
                                    + panel_nonsmoking_sample_ids),
        'smoking': smoking_sample_ids,
        'nonsmoking': nonsmoking_sample_ids,
        'smoking_plus': _ordered_unique(smoking_sample_ids
                                        + panel_smoking_sample_ids),
        'nonsmoking_plus': _ordered_unique(nonsmoking_sample_ids
                                           + panel_nonsmoking_sample_ids),
    }


def _resolve_sample_ids_for_key(key, default_sample_ids):
    """Apply optional replacement/add/drop files for one analysis key.

    Inputs:
    - key: one of `results_keys`
    - default_sample_ids: default sample IDs for `key`

    Outputs:
    - tuple[list[str], dict]: final sample IDs and override metadata

    Assumptions:
    - Override files, when present, are stored in `sample_id_override_dir`.
    - Replacement is applied before add/drop.
    """
    replace_file = os.path.join(sample_id_override_dir, f"{key}.txt")
    add_file = os.path.join(sample_id_override_dir, f"{key}__add.txt")
    drop_file = os.path.join(sample_id_override_dir, f"{key}__drop.txt")

    if os.path.exists(replace_file):
        sample_ids = _read_sample_ids(replace_file)
        mode = 'replace'
    else:
        sample_ids = list(default_sample_ids)
        mode = 'default'

    added_sample_ids = _read_sample_ids(add_file) if os.path.exists(add_file) else []
    dropped_sample_ids = _read_sample_ids(drop_file) if os.path.exists(drop_file) else []

    sample_ids = _ordered_unique(sample_ids + added_sample_ids)
    if dropped_sample_ids:
        dropped_sample_ids = set(dropped_sample_ids)
        sample_ids = [sample_id for sample_id in sample_ids
                      if sample_id not in dropped_sample_ids]

    override_metadata = {
        'mode': mode,
        'base_n': len(default_sample_ids),
        'final_n': len(sample_ids),
        'replace_file': replace_file if os.path.exists(replace_file) else None,
        'add_file': add_file if os.path.exists(add_file) else None,
        'drop_file': drop_file if os.path.exists(drop_file) else None,
        'n_added': len(added_sample_ids),
        'n_dropped': len(dropped_sample_ids) if isinstance(dropped_sample_ids, set) else 0,
    }

    return sample_ids, override_metadata


def get_active_sample_id_overrides():
    """Return a copy of currently active sample-ID override metadata."""
    return {
        key: metadata.copy()
        for key, metadata in active_sample_id_overrides.items()
    }


def frameshift_table_is_available():
    """Return True when an `indels_per_sample` table exists for this run."""
    return indels_per_sample is not None


def get_key_filtered_indel_db(key):
    """Return the cohort-filtered indel table for `key`.

    Inputs:
    - key: analysis cohort key from `results_keys`

    Outputs:
    - pandas.DataFrame: indel table filtered to the cohort

    Assumptions:
    - `refresh_key_filtered_dbs()` has already been called in the current
      process.
    """
    if key not in key_filtered_indel_dbs:
        raise ValueError(f"No indel table is available for result key: {key}")
    return key_filtered_indel_dbs[key]


def _validate_indel_table_alignment(indel_table, gene_table):
    """Ensure that the indel table matches the genes table sample universe.

    Inputs:
    - indel_table: `indels_per_sample` DataFrame
    - gene_table: `genes_per_sample` DataFrame

    Outputs:
    - None. Raises `ValueError` if alignment assumptions are violated.

    Assumptions:
    - Row identity is defined by `Sample ID`, `Source`, and `Panel`.
    """
    identity_cols = ['Sample ID', 'Source', 'Panel']

    gene_identity = gene_table[identity_cols].astype(str)
    indel_identity = indel_table[identity_cols].astype(str)

    if len(gene_identity) != len(indel_identity):
        raise ValueError(
            "`indels_per_sample` and `genes_per_sample` do not contain the "
            "same number of samples."
        )

    gene_identity_index = pd.MultiIndex.from_frame(gene_identity)
    indel_identity_index = pd.MultiIndex.from_frame(indel_identity)

    if not gene_identity_index.equals(indel_identity_index):
        raise ValueError(
            "`indels_per_sample` and `genes_per_sample` do not have the same "
            "sample rows in the same order."
        )

    missing_genes = [
        gene for gene in gene_table.columns
        if gene not in identity_cols and gene not in indel_table.columns
    ]
    if missing_genes:
        raise ValueError(
            "`indels_per_sample` is missing gene columns present in "
            f"`genes_per_sample`: {missing_genes}"
        )


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
    if key not in sample_ids_by_key:
        raise ValueError(f"Unrecognized result key: {key}")

    return db[db['Sample ID'].isin(sample_ids_by_key[key])]


def refresh_key_filtered_dbs():
    """Rebuild cohort definitions and filtered genotype tables from disk.

    Outputs:
    - dict[str, pandas.DataFrame]: filtered genotype tables for each result key

    Assumptions:
    - `genes_per_sample_file_name` exists in the current run root.
    - Optional cohort overrides are stored in `sample_id_override_dir`.
    """
    global genes_per_sample
    global indels_per_sample

    default_sample_ids_by_key = _build_default_sample_ids_by_key()
    resolved_sample_ids = {}
    override_summary = {}

    for key in results_keys:
        sample_ids, metadata = _resolve_sample_ids_for_key(
            key,
            default_sample_ids_by_key[key])
        resolved_sample_ids[key] = sample_ids
        if any(metadata[file_field] is not None
               for file_field in ['replace_file', 'add_file', 'drop_file']):
            override_summary[key] = metadata

    genes_per_sample = pd.read_csv(genes_per_sample_file_name, low_memory=False)

    if os.path.exists(indels_per_sample_file_name):
        indels_per_sample = pd.read_csv(indels_per_sample_file_name, low_memory=False)
        _validate_indel_table_alignment(indels_per_sample, genes_per_sample)
    else:
        indels_per_sample = None

    sample_ids_by_key.clear()
    sample_ids_by_key.update(resolved_sample_ids)

    active_sample_id_overrides.clear()
    active_sample_id_overrides.update(override_summary)

    key_filtered_dbs.clear()
    key_filtered_dbs.update({
        key: filter_db_for_key(key, genes_per_sample)
        for key in results_keys
    })

    key_filtered_indel_dbs.clear()
    if indels_per_sample is not None:
        key_filtered_indel_dbs.update({
            key: filter_db_for_key(key, indels_per_sample)
            for key in results_keys
        })

    return key_filtered_dbs


refresh_key_filtered_dbs()
