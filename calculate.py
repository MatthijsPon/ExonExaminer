#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Perform calculations on the parsed exons
"""
import collections

import pandas as pd
from plots import *


def exon_trans_gene_amount(exon_df):
    """Calculate the amount of unique exons, transcripts and genes.

    :param exon_df: pandas dataframe object
    :return: int, no. of exons, transcripts and genes
    """
    # The no. of rows reflect the amount of unique exons
    exons = len(exon_df.index)

    # For the transcripts, all the unique values within the transcript lists have to be taken
    # Adapted from: https://towardsdatascience.com/dealing-with-list-values-in-pandas-dataframes-a177e534f173
    transcripts = pd.Series([x for trans_list in exon_df["transcriptID"] for x in trans_list]).nunique()

    # Gather all unique gene names
    genes = exon_df["geneID"].nunique()

    return exons, transcripts, genes


def prime_dataframes(exon_df):
    """Return four different dataframes, one for each combination of the 5' and 3' variables.

    :param exon_df: pandas dataframe object
    :return: dictionary of dataframes, (bool, bool): dataframes
    """
    df_dict = {}
    for idx, df in exon_df.groupby(["five_prime", "three_prime"]):
        df_dict[idx] = df
    return df_dict


def bucket_exons(exon_df, bucket_list, index):
    """Bucket exons based on a certain size.

    :param exon_df: pandas dataframe object
    :param bucket_list: list of ints, upper border of each bucket
    :param index: pandas index of group
    :return: exon dataframe with bucket column added
    """
    name_dict = {
        (False, False): "internal exons",
        (False, True): "3' exons",
        (True, False): "5' exons",
        (True, True): "double-prime exons",
    }
    # Check which exon is the biggest
    max_size = exon_df["size"].max()
    print("Biggest exon of {}: {} bp".format(name_dict[index], max_size))
    # Check if max exon fits in buckets
    if bucket_list[-1] < max_size:
        raise ValueError("Max bucket size does not fit largest exon.")

    exon_df["bucket"] = pd.cut(exon_df["size"], bucket_list)
    return exon_df


def bucket_exons_and_plot(prime_df_dict, exon_buckets):
    """Use the bucket_exons() function and create bar plots

    :param prime_df_dict: pandas df, exon df containing the prime indicators
    :param exon_buckets: list of ints, bucket borders
    :return: dict of series, (5' bool, 3' bool) coupled to bucketed exon counts
    """

    plot_dict = {}
    return_dict = {}
    for idx, value in prime_df_dict.items():
        bucket_df = bucket_exons(value, exon_buckets, idx)
        total = len(bucket_df.index)
        bucket_df = bucket_df.groupby(["bucket"]).size()
        return_dict[idx] = bucket_df
        perc_buckets = pd.Series(bucket_df/total)
        plot_dict[idx] = perc_buckets

    bucket_barplot(plot_dict, "./images/bucket_distribution")
    return return_dict


def largest_exon_genes(exon_df, quantity):
    """Return the gene names of the genes that contain the N largest exons in descending order.

    :param exon_df: pandas dataframe, exon dataframe from parse function
    :param quantity: amount of genes to return
    :return: list, genes ordered on largest exon size
    """
    subset = exon_df.nlargest(quantity, "size", keep="all")
    gene_list = subset["geneID"].tolist()
    return gene_list


def query_largest_exon_genes(gene_list, expression_file, filename):
    """"""
    # Determine which items are duplicates
    duplicates = [(item, count) for item, count in collections.Counter(gene_list).items() if count > 1]
    # Read in gene expression data
    expression_df = pd.read_csv(expression_file, skiprows=2, header=0, sep="\t")
    # Take a subset containing only the genes from the [gene_list]
    gene_expr = expression_df.loc[expression_df["Description"].isin(gene_list)]
    # This determines which genes are not included in the expression data

    excluded = list(set(gene_expr["Description"]).symmetric_difference(set(gene_list)))
    expression_heatmap(gene_expr, filename)
    return duplicates, excluded





