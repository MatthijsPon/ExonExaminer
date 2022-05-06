#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: This script parses an ENSEMBL gff3 file and determines the average
             inclusion of exons compared to their length.
"""

from transcript_expression import determine_target_exon, yield_single_gene
from matplotlib import pyplot as plt

import argparse as arg
import pandas as pd
import scipy.stats as stats

import statistics
import os


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine the average inclusion of exons compared to their length.")
    parser.add_argument("gff3_file", type=str,
                        help="gff3 file to parse")
    parser.add_argument("out_dir", type=str,
                        help="output directory (incl. trailing \"/\"), default: ./data/out/")
    parser.add_argument("sizes_of_interest", type=int, nargs="+",
                        help="size cutoff(s) for selecting exon sizes")
    args = parser.parse_args()
    return args.gff3_file, args.out_dir, args.sizes_of_interest


def get_internal_exons(gene_df, remove=None):
    """Create a list of all internal exons in the gene dataframe

    :param gene_df: pandas dataframe, gene dataframe
    :param remove: str, exon ID (default: None)
    :return: set of str, exon IDs of internal exons (otherwise None)

    Select exons per transcript and determine which are internal. An internal exon has a non-prime position in at
    least one transcript. Returns all unique exon IDs.
    """
    non_prime = set()
    only_exons = gene_df[gene_df["type"] == "exon"]

    for idx, group in only_exons.groupby("parent"):
        exons_transcript = list(group["id"])
        exons_transcript.pop(0)  # Remove the first exon
        if len(exons_transcript) >= 1:  # Cannot pop twice if there is only one exon
            exons_transcript.pop(-1)  # Remove the last exon
        for exon in exons_transcript:
            non_prime.add(exon)
    # Only return a value when there are internal exons
    if len(non_prime) >= 1:
        if remove:
            non_prime.remove(remove)
            if len(non_prime) == 0:
                return None
        return non_prime
    return None


def exon_transcript_presence(gene_df, target_exon):
    """Determine in which transcripts an exon is present

    :param gene_df: pandas dataframe, gene dataframe
    :param target_exon: str, exon ID
    :return: numpy list, IDs of transcripts containing target exon
    """
    return gene_df.loc[gene_df["id"] == target_exon, "parent"].nunique()


def ratio_transcripts_exon(gene_df, target_exon=None):
    """Calculate the ratio of transcripts containing exon vs all transcripts

    :param gene_df: pandas dataframe, gene dataframe
    :param target_exon: str, exon ID (default: None; will target largest internal exon)
    :return: float, ratio of transcripts containing exon
    """
    if not target_exon:
        target_exon = determine_target_exon(gene_df, target="largest_internal")
        if not target_exon:
            return None
    # Determine how many transcripts are present
    all_trans = gene_df.loc[gene_df["type"] == "transcript", "id"].nunique()
    # Determine how many transcripts contain exon
    target_trans = gene_df.loc[gene_df["id"] == target_exon, "parent"].nunique()
    target_ratio = float(target_trans) / all_trans
    # Get all non-target internal exons
    internal_exons = get_internal_exons(gene_df)
    if not internal_exons:
        return None
    internal_ratios = []
    for internal_exon in internal_exons:
        int_ex_pres = exon_transcript_presence(gene_df, internal_exon)
        internal_ratios.append(float(int_ex_pres) / all_trans)
    internal_ratio = statistics.fmean(internal_ratios)
    return target_ratio - internal_ratio


def investigate_normality_plots(data, column_of_interest, out_dir, filename):
    """Create a histogram and QQplot vs. normal distribution for investigating normality on the data

    :param data: pandas dataframe, dataframe containing data
    :param column_of_interest: str, column of interest in the dataframe
    :param out_dir: str, output directory
    :param filename: str, start of filename for images to be saved to
    :return: None, images are created
    """
    # Histogram
    plt.rcParams.update({"font.size": 20, "figure.figsize": (19.2, 10.8)})
    data[column_of_interest].hist()  # Histogram
    plt.savefig("{}/{}_histogram.png".format(out_dir, filename))
    plt.close()
    # QQ plot
    plt.rcParams.update({"font.size": 20, "figure.figsize": (19.2, 10.8)})
    stats.probplot(data[column_of_interest], dist="norm", plot=plt)  # QQ plot
    plt.savefig("{}/{}_QQ_plot_vs_normal_distribution.png".format(out_dir, filename))
    plt.close()


def statistical_information(result_df, out_file, out_dir, size):
    """Calculate Welsch T-test if possible and output some statistical information

    :param result_df: pandas dataframe, ratios and sizes
    :param out_file: str, file to output information to
    :param out_dir: str, output directory path
    :param size: int, size value to split the data on
    :return: None, all information is outputted into outfile and figure
    """
    total = None
    out_file.write("----- Exons of interest: {}+ bp -----\n".format(size))
    subset = result_df.loc[result_df["exon_size"] >= size]
    if size > 0:
        total = result_df[~result_df["exon_size"].isin(subset["exon_size"])]

    out_file.write("Total average expression difference {}: {:.2f}\n"
                   "".format(size, statistics.fmean(list(subset["ratio"]))))
    out_file.write("No. of internal exons examined: {}\n".format(len(subset.index)))

    # Get mean, median, mode, quartiles, etc. for target group
    info_sub = subset.describe().transpose()
    out_file.write("\nInformation target group:\n")
    out_file.write(info_sub.to_string())
    out_file.write("\n")

    if total is not None:
        # Get mean, median, mode, quartiles, etc. for leftover group
        info_total = total.describe().transpose()
        out_file.write("\nInformation rest group:\n")
        out_file.write(info_total.to_string())
        out_file.write("\n")

        # Perform Welch T-test
        welch_t, p_val = stats.ttest_ind(total["ratio"], subset["ratio"], equal_var=False)
        out_file.write("Welch t-statistic: {}, p-value: {}\n".format(welch_t, p_val))

    out_file.write("\n\n")
    plt.rcParams.update({"font.size": 20, "figure.figsize": (19.2, 10.8)})
    subset.plot(x="exon_size", y="ratio", style="o")
    plt.savefig("{}/internal_exon_ratio_{}+.png".format(out_dir, size))
    plt.close()


def main():
    """Main function."""
    gff_file, out_dir, sizes_of_interest = parse_arguments()

    exon_sizes = []
    ratios = []
    with open(gff_file) as file:
        for gene in yield_single_gene(file):
            internal_exons = get_internal_exons(gene)
            if internal_exons is None:
                continue  # If there are no internal exons, skip the ratios
            for exon in internal_exons:
                ratio = ratio_transcripts_exon(gene, target_exon=exon)
                if ratio:
                    ratios.append(ratio)
                    size = gene.loc[gene["id"] == exon, "size"].iloc[0]
                    exon_sizes.append(int(size))
    # Create DF from exon size and ratio
    results = pd.DataFrame({"exon_size": exon_sizes, "ratio": ratios})

    investigate_normality_plots(results, "ratio", out_dir, "ratio")

    # Output Welch's T-test and dataframe distributions to file
    with open("{}/statistical_information.txt".format(out_dir), "w+") as out_file:
        for size in sizes_of_interest:
            statistical_information(results, out_file, out_dir, size)


if __name__ == '__main__':
    main()
