"""
This module compiles and processes metagenomic data
It combines headers, checkM2 and GTDBTK data into a single dataframe
Creates quality categories based on Contamination and Completeness
labels contigs that are both large and circular
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# create an absolute path to folder, in which this .py file is located
HERE = Path(__file__).resolve().parent

# append path to folder into sys.path, so that parsers.py can be found
sys.path.append(str(HERE))

# import parsers.py module
from parsers import parse_myloasm_headers, parse_metamdbg_headers, parse_checkm2, parse_gtdbtk

def compile_data(paths: dict):
    """
    compiles all available data into a single dataframe
    adjusts Contamination and Completeness columns
    creates quality categories based on contamination and completeness
    analyses which contigs are large circular

    :param paths: dictionary of paths to all required files
    :return compiled and adjusted dataframe of all data
    """

    # load assembler headers
    print("Parsing assembler headers...")
    df_mylo = parse_myloasm_headers(paths['headers_mylo'])
    df_meta = parse_metamdbg_headers(paths['headers_meta'])

    # concatenate headers for both myloasm and metamdbg, resetting index
    df_main = pd.concat([df_mylo, df_meta], ignore_index=True)

    # load CheckM2 reports
    print("Parsing CheckM2 reports...")
    cm2_mylo = parse_checkm2(paths['checkm2_mylo'])
    cm2_meta = parse_checkm2(paths['checkm2_meta'])

    # concatenate reports for both myloasm and metamdbg, resetting index
    df_checkm = pd.concat([cm2_mylo, cm2_meta], ignore_index=True)

    # load GTDBTK reports for both bacteria and archaea
    print("Parsing GTDBTK reports...")
    gtdb_mylo = parse_gtdbtk(paths['gtdb_mylo_ar'], paths['gtdb_mylo_bac'])
    gtdb_meta = parse_gtdbtk(paths['gtdb_meta_ar'], paths['gtdb_meta_bac'])

    # concatenate reports for both myloasm and metamdbg, resetting index
    df_gtdb = pd.concat([gtdb_mylo, gtdb_meta], ignore_index=True)

    # merge data based on contig_id; how=left ensures that no data is lost
    print("Merging data into a single dataframe...")
    df_merged = df_main.merge(df_checkm, on='contig_id', how='left')
    df_merged = df_merged.merge(df_gtdb, on='contig_id', how='left')

    # replaces NA values in Completeness and Contamination
    completeness = df_merged['Completeness'].fillna(0)
    contamination = df_merged['Contamination'].fillna(100)

    # list of conditions and choices to determine quality
    conditions = [
        (completeness > 90) & (contamination < 5),
        (completeness > 50) & (contamination < 10)]
    choices = ['High', 'Medium']

    # apply conditions and choices on data and create a new column
    df_merged['quality_category'] = np.select(conditions, choices, default='Low')

    # create a new column based on length and circularity
    df_merged['is_large_circular'] = (
            (df_merged['length'] > 500_000) &
            (df_merged['is_circular'] == True))

    # replace any NaN values in Phylum with Unknown
    df_merged['Phylum'] = df_merged['Phylum'].fillna('Unknown')

    return df_merged

if __name__ == "__main__":
    results_dir = HERE.parent / 'results'

    # define all file paths in a dictionary
    file_paths = {
        # headers
        'headers_mylo': results_dir / 'myloasm_assembly_headers.txt',
        'headers_meta': results_dir / 'metamdbg_assembly_headers.txt',

        # CheckM2
        'checkm2_mylo': results_dir / 'checkm2/myloasm/quality_report.tsv',
        'checkm2_meta': results_dir / 'checkm2/metamdbg/quality_report.tsv',

        # GTDBTK
        'gtdb_mylo_ar': results_dir / 'gtdbtk/myloasm/classify/gtdbtk.ar53.summary.tsv',
        'gtdb_mylo_bac': results_dir / 'gtdbtk/myloasm/classify/gtdbtk.bac120.summary.tsv',
        'gtdb_meta_ar': results_dir / 'gtdbtk/metamdbg/classify/gtdbtk.ar53.summary.tsv',
        'gtdb_meta_bac': results_dir / 'gtdbtk/metamdbg/classify/gtdbtk.bac120.summary.tsv',}

    df_final = compile_data(file_paths)

    #quick look at the data as a validation
    print(df_final.head())
    print("\nCounts based on quality categories:")
    print(df_final['quality_category'].value_counts())
    print(df_final[["coverage", "length"]])
    print(f"Number of rows: {len(df_final)}")
    print("Counts based on assembler:")
    print(df_final['assembler'].value_counts())
    print("\nNumber of circular contigs:")
    print(df_final.groupby(['assembler', 'is_large_circular']).size())
    print(df_final["Phylum"].value_counts())
