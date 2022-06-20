#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Parse a gff3 file into a usable dataframe.
"""
import numpy as np
import pandas as pd

from helpers import yield_single_gene
import argparse as arg


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Parse a gff3 file into a pandas dataframe (.pickle) file.")
    parser.add_argument("gff3_file", type=str,
                        help="gff3 file to parse")
    parser.add_argument("out_file", type=str,
                        help="output filename")
    args = parser.parse_args()
    return args.gff3_file, args.out_file


def determine_primes(gene_df):
    """Create a list of all internal exons in the gene dataframe

    :param gene_df: pandas dataframe, gene dataframe
    :return: set of str, exon IDs of internal exons (otherwise None)

    Select exons per transcript and determine which are internal. An internal exon has a non-prime position in at
    least one transcript. Returns all unique exon IDs.
    """
    three_prime = set()
    five_prime = set()
    only_exons = gene_df[gene_df["type"] == "exon"]

    for idx, group in only_exons.groupby("parent"):
        status = group["strand"].iloc[0]
        lowest = group.loc[group["start"] == group["start"].min(), "id"].iloc[0]
        highest = group.loc[group["stop"] == group["stop"].max(), "id"].iloc[0]
        if status == "+":
            five_prime.add(lowest)
            three_prime.add(highest)
        else:
            five_prime.add(highest)
            three_prime.add(lowest)
    # Add boolean for five and three prime to dataframe
    gene_df["three_prime"] = np.where(gene_df["id"].isin(three_prime), True, False)
    gene_df["five_prime"] = np.where(gene_df["id"].isin(five_prime), True, False)
    return gene_df


def merge_exons(gff_df):
    """Merge exons in a dataframe.

    :param gff_df: pandas dataframe object, gff3 dataframe
    :return: gff_df, with exons merged
    """
    columns = list(gff_df.columns)
    columns.remove("parent")  # only parent has to be removed, since equal exons have equal IDs
    gff_df = gff_df.groupby(columns, dropna=False)["parent"].apply(list).reset_index()
    return gff_df


def parse_gff_file(file_object):
    """Parse a gff file into a pandas dataframe.

    :param file_object: iterable, open gff3 file object
    :return: pandas dataframe object, parsed gff3 file
    """
    gene_list = []
    for gene in yield_single_gene(file_object):
        gene = determine_primes(gene)
        gene_list.append(gene)
    gff_df = pd.concat(gene_list, ignore_index=True)
    gff_df = merge_exons(gff_df)  # structure is lost here, due to the index resetting
    return gff_df


def main():
    """Main function"""
    gff, out_file = parse_arguments()
    with open(gff) as file:
        gff = parse_gff_file(file)
    gff.to_pickle(out_file)
    return None


if __name__ == "__main__":
    main()
