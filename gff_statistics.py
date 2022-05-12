#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

import os
import time

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def yield_gff_line(file_object):
    """Iterator, Yields a single non-comment line of a tsv gff file

    :param file_object: iterable, open file object or list of lines
    :return: yields a list, a gff line split on \t
    """
    for line in file_object:
        if line.startswith("#"):
            # Skip the comment line, which is not of interest
            continue
        line = line.strip()
        yield line.split("\t")


def determine_gff(gff_8, identifier, further_split=None):
    """Select a specific object from the information column (column 8 in 0-count)

    :param gff_8: list of str, information column, split
    :param identifier: str, identifier of column of interest (e.g. ID=). Everything
                       identifier will be removed from return string
    :param further_split: str, seperator on which the column should be split on
                          after removal of the identifier (default: None)
    :return: str, information following the identifier, otherwise None
    """
    for item in gff_8.split(";"):
        if item.startswith(identifier):
            if further_split:
                return item[len(identifier):].split(further_split)
            return item[len(identifier):]
    return None


def yield_single_gene(file_object):
    """Yield a single gene from a gff3 file, including related protein coding transcripts and exons.

    :param file_object: open file object, gff3 file
    :return: pandas dataframe, gene information containing type, id, parent, size and alias columns
    """
    gene_dict = {
        "type": [],
        "id": [],
        "chromosome": [],
        "start": [],
        "stop": [],
        "strand": [],
        "parent": [],
        "size": [],
    }

    # A dict needs all three to be complete
    gene_id = ""
    transcript = []
    exon = []

    for line in yield_gff_line(file_object):
        if line[2] == "gene":
            if gene_id and transcript and exon:
                yield pd.DataFrame.from_dict(gene_dict)
            # Reset stats
            gene_dict = {
                "type": [],
                "id": [],
                "chromosome": [],
                "start": [],
                "stop": [],
                "strand": [],
                "parent": [],
                "size": [],
            }
            gene_id = determine_gff(line[8], "ID=gene:")
            transcript = []
            exon = []
            # Add data to dict
            gene_dict["type"].append("gene")
            gene_dict["id"].append(gene_id)
            gene_dict["chromosome"].append(line[0])
            gene_dict["start"].append(line[3])
            gene_dict["stop"].append(line[4])
            gene_dict["strand"].append(line[6])
            gene_dict["parent"].append(None)
            gene_dict["size"].append(None)
        elif line[2] == "mRNA":
            # only add if there is a parent gene and the mRNA is protein_coding
            if determine_gff(line[8], "Parent=gene:") == gene_id and \
                    determine_gff(line[8], "biotype=") == "protein_coding":
                trans_id = determine_gff(line[8], "ID=transcript:")
                transcript.append(trans_id)
                # Add data to dict
                gene_dict["type"].append("transcript")
                gene_dict["id"].append(trans_id)
                gene_dict["chromosome"].append(line[0])
                gene_dict["start"].append(line[3])
                gene_dict["stop"].append(line[4])
                gene_dict["strand"].append(line[6])
                gene_dict["parent"].append(gene_id)
                gene_dict["size"].append(None)
        elif line[2] == "exon":
            if determine_gff(line[8], "Parent=transcript:") in transcript:
                exon_id = determine_gff(line[8], "Name=")
                exon.append(exon_id)
                gene_dict["type"].append("exon")
                gene_dict["id"].append(exon_id)
                gene_dict["chromosome"].append(line[0])
                gene_dict["start"].append(line[3])
                gene_dict["stop"].append(line[4])
                gene_dict["strand"].append(line[6])
                gene_dict["parent"].append(determine_gff(line[8], "Parent=transcript:"))
                gene_dict["size"].append((int(line[4]) - int(line[3])))

    if gene_id and transcript and exon:
        yield pd.DataFrame.from_dict(gene_dict)


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


def exon_len_histograms(exon_df, out_folder, fig_name, bins, max_height=0):
    """Make histograms based on length

    :param exon_df: pandas dataframe, exon dataframe
    :param out_folder: str, path to output folder (excl. trailing "/")
    :param fig_name: str, name of output file (excl. ".png")
    :param bins: list of ints, bin borders
    :param max_height: int, determine max height of Y-axis
    :return: None, function creates .png file
    """

    # Create the plot and set name and axis titles
    exon_df["size"].plot.hist(bins=bins, title=fig_name)
    plt.xlabel("exon size (bp)")
    if max_height > 0:
        plt.ylim([0, max_height])

    plt.savefig("{}/{}.png".format(out_folder, fig_name))
    # Close figure, otherwise object keeps existing and bins are added cumulatively
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

    # Exon length distribution for all types
    bins = [i for i in range(0, 1000, 10)]
    exon_len_histograms(gff.loc[gff["type"] == "exon"], OUTDIR, "histogram_exonlen_0-1000bp", bins)
    bins = [i for i in range(1000, 20000, 10)]
    exon_len_histograms(gff.loc[gff["type"] == "exon"], OUTDIR, "histogram_exonlen_1000-20000bp", bins)
    # Exon length distribution per type
    bins = [i for i in range(0, 10000, 10)]
    name_dict = {
        (False, False): "Internal exons",
        (False, True): "3' exons",
        (True, False): "5' exons",
        (True, True): "Single exon genes"
    }
    for idx, group in gff.loc[gff["type"] == "exon"].groupby(["five_prime", "three_prime"]):
        exon_len_histograms(group, OUTDIR, "histogram_exonlen_0-10000bp_{}".format(name_dict[idx]), bins)

    # Statistics
    stat_list = [gff.loc[gff["type"] == "exon"].describe(percentiles=[0.25, 0.5, 0.75, 0.9]).transpose()]
    stat_list[0]["exon_type"] = "All"
    for idx, group in gff.loc[gff["type"] == "exon"].groupby(["five_prime", "three_prime"]):
        stat_list.append(group.describe(percentiles=[0.25, 0.5, 0.75, 0.9]).transpose())
        stat_list[-1]["exon_type"] = name_dict[idx]

    info_pd = pd.concat(stat_list).set_index("exon_type")
    info_pd.to_csv("{}/exon_statistics.txt".format(OUTDIR), sep="\t", mode="w+")


if __name__ == '__main__':
    main()