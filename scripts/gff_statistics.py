#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

import os
import pandas as pd
from matplotlib import pyplot as plt
import argparse as arg


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Use a parsed gff3 file to calculate gff statistics.")
    parser.add_argument("pickle", type=str,
                        help="gff dataframe pickle file")
    parser.add_argument("out_dir", type=str,
                        help="output directory")
    args = parser.parse_args()
    return args.pickle, args.out_dir


def exon_len_histograms(exon_df, out_folder, fig_name, bins, title=None, max_height=0):
    """Make histograms based on length

    :param exon_df: pandas dataframe, exon dataframe
    :param out_folder: str, path to output folder (excl. trailing "/")
    :param fig_name: str, name of output file (excl. ".png")
    :param bins: list of ints, bin borders
    :param title: title of the figure
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
    # Parse args
    pickle, out_dir = parse_arguments()

    # Check if pickle exists, otherwise
    if not os.path.exists(pickle):
        raise FileNotFoundError("Pickle file was not found.")
    gff = pd.read_pickle(pickle)

    gff = gff.loc[gff["chromosome"] != "MT"]
    # Get no. of genes, transcripts, exons
    genes = len(gff.loc[gff["type"] == "gene"].index)
    trans = len(gff.loc[gff["type"] == "transcript"].index)
    exons = len(gff.loc[gff["type"] == "exon"].index)
    with open("{}/exon_statistics.txt".format(out_dir), "w+") as outfile:
        outfile.write("# General information:\nNo. of genes: {}\nNo. of transcripts: {}\nNo. of exons: {}\n\n"
                      "".format(genes, trans, exons))

    # Create histogram of transcripts/gene and exons/transcript
    option = {"gene": "Genes per chromosome",
              "transcript": "Transcripts per gene",
              "exon": "Exons per transcripts"}
    output = {}
    for item in ["gene", "transcript", "exon"]:
        temp_data = gff.loc[gff["type"] == item]
        # Calculate
        if item == "gene":
            temp_data = pd.Series([item for item in temp_data.chromosome])
            temp_data = temp_data.value_counts()
        else:
            temp_data = pd.Series([item for sublist in temp_data.parent for item in sublist])
            temp_data = temp_data.value_counts()
        to_file = temp_data.describe()
        to_file["mode"] = temp_data.mode().iloc[0]
        to_file["median"] = temp_data.median()
        # Add to dict for printing to file
        output[item] = to_file
        # create histogram
        bins = [i for i in range(0, 40, 1)]
        plt.hist(temp_data, bins=bins)
        plt.xlabel("{}".format(option[item]))
        plt.ylabel("Count")
        plt.savefig("{}/histogram_{}.png".format(out_dir, option[item]))
        plt.close()

    with open("{}/transcriptgene_exontranscript_information.txt".format(out_dir), "w+", newline="") as file:
        for item, value in option.items():
            file.write("{}\n".format(value))
            output[item].to_csv(file, sep="\t")
            file.write("\n")

    # Exon length distribution for all types
    bins = [i for i in range(0, 1000, 10)]
    exon_len_histograms(gff.loc[gff["type"] == "exon"], out_dir, "histogram_exonlen_0-1000bp", bins)
    bins = [i for i in range(1000, 20000, 10)]
    exon_len_histograms(gff.loc[gff["type"] == "exon"], out_dir, "histogram_exonlen_1000-20000bp", bins)
    # Exon length distribution per type
    name_dict = {
        (False, False): "Internal exons",
        (False, True): "3' exons",
        (True, False): "5' exons",
        (True, True): "Single exon genes"
    }
    for idx, group in gff.loc[gff["type"] == "exon"].groupby(["five_prime", "three_prime"]):
        bins = [i for i in range(0, 1000)]
        exon_len_histograms(group, out_dir, "histogram_exonlen_0-1000bp_{}".format(name_dict[idx]),
                            bins, title=name_dict[idx])
        bins = [i for i in range(0, 10000, 10)]
        exon_len_histograms(group, out_dir, "histogram_exonlen_0-10000bp_{}".format(name_dict[idx]),
                            bins, title=name_dict[idx])

    # Statistics
    percentiles = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    stat_list = [gff.loc[gff["type"] == "exon"].describe(percentiles=percentiles).transpose()]
    stat_list[0]["exon_type"] = "All"
    for idx, group in gff.loc[gff["type"] == "exon"].groupby(["five_prime", "three_prime"]):
        stat_list.append(group.describe(percentiles=percentiles).transpose())
        stat_list[-1]["exon_type"] = name_dict[idx]
    # Save statistics to txt file
    info_pd = pd.concat(stat_list).set_index("exon_type")
    info_pd.to_csv("{}/exon_statistics.txt".format(out_dir), sep="\t", mode="a+")

    # Save the 5th and 95th quartile to file
    Q5, Q95 = info_pd.loc["Internal exons"].loc[["5%", "95%"]]
    Q5, Q95 = int(Q5), int(Q95)

    with open("{}/5_95_quartiles.txt".format(out_dir), "w+") as file:
        file.write("{}\t{}".format(Q5, Q95))
    return None


if __name__ == '__main__':
    main()
