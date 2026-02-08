"""
This module contains different parsing functions used to extract data from .txt and .tsv
files into dataframes
"""

import pandas as pd

def parse_myloasm_headers(file_path):
    """
    reads myloasm header and extracts all data into a dataframe
    :param file_path: system path to desired .txt file
    :return: dataframe containing extracted columns of interest
    """
    with open(file_path, 'r', encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.startswith(">")]

    # insert each row into dataframe with single column raw_header
    df = pd.DataFrame(lines, columns=["raw_header"])

    # regex pattern for data extraction from raw_header
    pattern = \
        r'^>(?P<contig_id>.+?)_len-(?P<length>\d+)_circular-(?P<circular>\w+)_depth-(?P<coverage>[\d\.]+)'

    # create new dataframe of extracted data from raw_header
    extracted = df["raw_header"].str.extract(pattern)

    # convert extracted length into its own column and a numeric type
    extracted["length"] = pd.to_numeric(extracted["length"])

    # convert extracted length into its own column and a numeric type
    extracted["coverage"] = pd.to_numeric(extracted["coverage"])

    # convert circularity from "yes" and "possible" to True or False
    extracted["is_circular"] = extracted["circular"].isin(["yes", "possible"])

    # name the assembler
    extracted["assembler"] = "myloasm"

    return extracted[["contig_id", "assembler", "length", "coverage", "is_circular"]]

def parse_metamdbg_headers(file_path):
    """
    reads metamdbg header and extracts all data into a dataframe
    :param file_path: system path to desired .txt file
    :return: dataframe containing extracted columns of interest
    """
    with open(file_path, 'r', encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.startswith(">")]

    # insert each row into dataframe with single column raw_header
    df = pd.DataFrame(lines, columns=["raw_header"])

    # extract contig_id using regex (stuff after ">")
    df['contig_id'] = df['raw_header'].str.extract(r'^>(\S+)')

    # extract length as float using regex
    df['length'] = df['raw_header'].str.extract(r'(?:length)[=_:]([\d]+)').astype(float)

    # extract coverage as float using regex
    df['coverage'] = df['raw_header'].str.extract(r'(?:coverage)[=_:]([\d\.]+)').astype(float)

    # extract circularity using regex
    circ_extracted = df['raw_header'].str.extract(r'(?:circular)[=_:](\w+)')

    # convert "yes" and "no" into True and False values
    df["is_circular"] = circ_extracted[0].str.lower().isin(["yes"])

    # name the assembler column
    df["assembler"] = "metaMDBG"

    return df[['contig_id', 'assembler', 'length', 'coverage', 'is_circular']]

def parse_checkm2(file_path):
    """
    finds and reads all checkm2 files - then extracts data of interest
    :param file_path: system path to desired .tsv file
    :return: dataframe containing extracted columns of interest
    """
    # save file contents into dataframe
    df = pd.read_csv(file_path, sep="\t")

    # rename Name so its consistent
    df = df.rename(columns={"Name": "contig_id"})

    return df[["contig_id", "Completeness", "Contamination"]]

def parse_gtdbtk(file_bac, file_ar):
    """
    finds and reads all gtdbtk files - then extracts data of interest
    uses a regex to extract Phylum information and ignores all Unclassified
    :param file_bac: system path to desired .tsv file for bacteria
    :param file_ar: system path to desired .tsv file for archaea
    :return: dataframe containing extracted columns of interest
    """
    # save contents of the files into dataframe
    df_bac = pd.read_csv(file_bac, sep="\t")
    df_ar = pd.read_csv(file_ar, sep="\t")

    #concatenate dataframe together, resetting indexes
    df = pd.concat([df_bac, df_ar], ignore_index=True)

    #rename user_genome so its consistent
    df = df.rename(columns={"user_genome": "contig_id"})

    #from classification extract Phylum using regex
    df["Phylum"] = df["classification"].str.extract(r"p__([^;]+)")

    #remove any Phylum values that are either NaN or Unclassified
    df = df.dropna(subset=["Phylum"])
    df = df[~df['Phylum'].str.contains('Unclassified')]

    return df[['contig_id', 'Phylum']]
