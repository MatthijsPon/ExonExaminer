#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from helpers import yield_single_gene


def determine_primes(gene_df):
    """Create a list of all internal exons in the gene dataframe

    :param gene_df: pandas dataframe, gene dataframe
    :param remove: str, exon ID (default: None)
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
    """"""
    columns = list(gff_df.columns)
    columns.remove("parent")
    gff_df = gff_df.groupby(columns, dropna=False)["parent"].apply(list).reset_index()
    return gff_df


def parse_gff_file(file_object):
    """"""
    gene_list = []
    for gene in yield_single_gene(file_object):
        gene = determine_primes(gene)
        gene_list.append(gene)
    gff_df = pd.concat(gene_list, ignore_index=True)
    gff_df = merge_exons(gff_df)  # structure is lost here, due to the index resetting
    return gff_df


def row_statistics(gff_df, row, out_folder, out_file):
    """"""


def exon_len_histograms(exon_df, out_folder, fig_name, bins, title=None, max_height=0):
    """Make histograms based on length

    :param exon_df: pandas dataframe, exon dataframe
    :param out_folder: str, path to output folder (excl. trailing "/")
    :param fig_name: str, name of output file (excl. ".png")
    :param bins: list of ints, bin borders
    :param max_height: int, determine max height of Y-axis
    :return: None, function creates .png file
    """
    if not title:
        title = fig_name

    # Create the plot and set name and axis titles
    exon_df["size"].plot.hist(bins=bins, title=title)
    plt.xlabel("Exon size (nt)")
    if max_height > 0:
        plt.ylim([0, max_height])

    plt.savefig("{}/{}.png".format(out_folder, fig_name))
    # Close figure, otherwise object keeps existing
    plt.close()


def main():
    """Main function."""
    # Replace with parser if needed:
    OUTDIR = "data/out/gff_statistics"
    GFF = "data/Homo_sapiens.GRCh38.106.chr.gff3"
    TEMP = "data/temp"

    # Actual code
    if not os.path.exists("{}/parsed_human_gff.pickle".format(TEMP)):
        with open(GFF) as file:
            gff = parse_gff_file(file)
        gff.to_pickle("{}/parsed_human_gff.pickle".format(TEMP))
    else:
        gff = pd.read_pickle("{}/parsed_human_gff.pickle".format(TEMP))

    # Create histogram of transcripts/gene and exons/transcript
    option = {"transcript": "transcripts per gene",
              "exon": "exons per transcripts"}
    output = {}
    for item in ["transcript", "exon"]:
        temp_data = gff.loc[gff["type"] == item]
        # Calculate
        temp_data = pd.Series([item for sublist in temp_data.parent for item in sublist])
        temp_data = temp_data.value_counts()
        # Add to dict for printing to file
        output[item] = temp_data.describe()
        # create histogram
        bins = [i for i in range(0, 40, 1)]
        plt.hist(temp_data, bins=bins)
        plt.xlabel("{}".format(option[item]))
        plt.savefig("{}/histogram_{}.png".format(OUTDIR, option[item]))
        plt.close()

    with open("{}/transcriptgene_exontranscript_information.txt".format(OUTDIR), "w+", newline="") as file:
        for item, value in option.items():
            file.write("{}\n".format(value))
            output[item].to_csv(file, sep="\t")
            file.write("\n")


    # Exon length distribution for all types
    bins = [i for i in range(0, 1000, 10)]
    exon_len_histograms(gff.loc[gff["type"] == "exon"], OUTDIR, "histogram_exonlen_0-1000bp", bins)
    bins = [i for i in range(1000, 20000, 10)]
    exon_len_histograms(gff.loc[gff["type"] == "exon"], OUTDIR, "histogram_exonlen_1000-20000bp", bins)
    # Exon length distribution per type
    name_dict = {
        (False, False): "Internal exons",
        (False, True): "3' exons",
        (True, False): "5' exons",
        (True, True): "Single exon genes"
    }
    for idx, group in gff.loc[gff["type"] == "exon"].groupby(["five_prime", "three_prime"]):
        bins = [i for i in range(0, 1000)]
        exon_len_histograms(group, OUTDIR, "histogram_exonlen_0-1000bp_{}".format(name_dict[idx]),
                            bins, title=name_dict[idx])
        bins = [i for i in range(0, 10000, 10)]
        exon_len_histograms(group, OUTDIR, "histogram_exonlen_0-10000bp_{}".format(name_dict[idx]),
                            bins, title=name_dict[idx])

    # Statistics
    stat_list = [gff.loc[gff["type"] == "exon"].describe(percentiles=[0.25, 0.5, 0.75, 0.9]).transpose()]
    stat_list[0]["exon_type"] = "All"
    for idx, group in gff.loc[gff["type"] == "exon"].groupby(["five_prime", "three_prime"]):
        stat_list.append(group.describe(percentiles=[0.25, 0.5, 0.75, 0.9]).transpose())
        stat_list[-1]["exon_type"] = name_dict[idx]
    # Save statistics to txt file
    info_pd = pd.concat(stat_list).set_index("exon_type")
    info_pd.to_csv("{}/exon_statistics.txt".format(OUTDIR), sep="\t", mode="w+")


if __name__ == '__main__':
    main()
