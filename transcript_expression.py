#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import argparse as arg
import math
import os
import time
import pandas as pd
import gzip
import statistics as stats
from matplotlib import pyplot as plt

from helpers import get_internal_exons, yield_single_gene, determine_target_exon


def parse_arg():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Investigate transcript expression of transcripts containing long exons.")
    parser.add_argument("-o", "--out", default="./data/out/", type=str,
                        help="output directory (incl. trailing \"/\"), default: ./data/out/")
    parser.add_argument("gff3_file", type=str,
                        help="gff3 file to parse")
    parser.add_argument("expression_file", type=str,
                        help="GTEx transcript expression file")
    parser.add_argument("expression_companion", type=str,
                        help="GTEx transcript expression companion file, describing sampleIDs and tissues")
    parser.add_argument("-c", "--complexity", default=1, type=int, choices=[0, 1, 2],
                        help="determine the tissue depth (0: None (avg. exp per gene, 1: general tissue (e.g. brain), "
                             "2: specific tissue (e.g. brian - frontal cortex)")
    args = parser.parse_args()
    return args.gff3_file, args.expression_file, args.expression_companion, args.out, args.complexity


def create_averaged_expression(exp_file, exp_comp, out_dir, complexity=1):
    """Create an expression file with averaged expression over complexity 1 or 2

    :param exp_file: str, expression file path
    :param exp_comp: str, expression companion file path
    :param out_dir: str, output directory
    :param complexity: int, 1 or 2, complexity of the tissue (e.g. 1 - brain; 2 - brain, frontal cortext)
    :return: None, expression file is created
    """
    # Read in column information
    column_of_interest = "SMTS"
    if complexity == 2:
        column_of_interest = "SMTSD"

    companion_df = pd.read_csv(exp_comp, delimiter="\t", header=0)
    # Subset df to only relevant columns
    companion_df = companion_df[["SAMPID", column_of_interest]]
    options = tuple([item for item in companion_df[column_of_interest].unique() if item != "Bone Marrow"])
    option_columns = {x: [] for x in options}
    not_measured = set()
    # Open gzipped file and file to output data to
    with gzip.open(exp_file, "rt") as file, open("{}averaged_expression_complexity{}.gct"
                                                 "".format(out_dir, complexity), "w+") as out_f:
        count = 0
        for line in file:
            if count < 3:
                count += 1
                if count == 3:
                    # This is the header line:
                    columns = line.strip().split("\t")
                    for idx, item in enumerate(columns):
                        # idx 1 and 2 are transcript and gene ID respectively
                        if idx > 1:
                            # Couple column indices to cell types
                            col_value = companion_df.loc[companion_df["SAMPID"] == item, column_of_interest]
                            col_value = col_value.iloc[0]
                            option_columns[col_value].append(idx)
                    out_f.write("trans_id\t{}\n".format("\t".join(options)))
            else:
                count += 1
                if count % 1000 == 0:  # There are 199324 lines. This shows progress
                    print("Line {} of 199324".format(count))
                line = line.strip().split("\t")
                to_write = []
                for option in options:
                    # Get all the items per cell type from the line using the column indices
                    items = [float(line[i]) for i in option_columns[option]]
                    if items:
                        average = stats.fmean(items)
                    else:
                        not_measured.add(option)
                        average = 0.0
                    # Print the averages as floats with 2 decimals
                    to_write.append("{:.2f}".format(average))
                # Only take the transcript ID without the version number
                out_f.write("{}\t{}\n".format(str(line[0]).split(".")[0], "\t".join(to_write)))
    return None


def create_comp0_averaged_expression(exp_file, out_dir):
    """Create expression file for complexity 0

    :param exp_file: str, expression filepath
    :param out_dir: str, output directory
    :return: None, expression file is created
    """
    # Open gzipped file and file to output data to
    with gzip.open(exp_file, "rt") as file, open("{}averaged_expression_complexity0.gct"
                                                 "".format(out_dir), "w+") as out_f:
        count = 0
        for line in file:
            if count < 3:
                count += 1
                if count == 3:
                    # Replace header line
                    out_f.write("trans_id\taverage\n")
            else:
                line = line.strip().split("\t")
                average = stats.fmean([float(x) for x in line[2:]])
                # Only take the transcript ID without the version number
                out_f.write("{}\t{:.2f}\n".format(str(line[0]).split(".")[0], average))
    return None


def acquire_expression_dataframe(exp_file, exp_comp, out_dir, complexity=0):
    """Read in averaged expression file or create it

    :param exp_file: str, expression file path
    :param exp_comp: str, expression companion file path
    :param out_dir: str, output directory
    :param complexity: int, 1 or 2, complexity of the tissue (e.g. 1 - brain; 2 - brain, frontal cortext)
    :return: pandas dataframe, expression dataframe
    """
    save_location = "{}averaged_expression_complexity{}.gct".format(out_dir, complexity)
    if not os.path.exists(save_location):
        # Create averaged expression file
        if complexity == 0:
            create_comp0_averaged_expression(exp_file, out_dir)
        else:
            create_averaged_expression(exp_file, exp_comp, out_dir, complexity)
    # Read in file and return
    return pd.read_csv(save_location, header=0, sep="\t")


def ratio_expression_exon(gene_df, expression_df, target_exon=None):
    """Calculate the ratio between expression of transcripts including target exon and all transcripts.

    :param gene_df: pandas dataframe, dataframe of a single gene
    :param expression_df: pandas dataframe, expression dataframe
    :param target_exon: str, ID of target exon (default: None)
    :return: float, ratio
    """

    # If no target exon ID is given, determine which exon is the largest and take that ID
    if not target_exon:
        target_exon = determine_target_exon(gene_df, target="largest_internal")
        if not target_exon:
            return None
    # Determine in which transcripts this exon is present
    target_trans = list(gene_df.loc[gene_df["id"] == target_exon, "parent"])
    all_trans = list(gene_df.loc[gene_df["type"] == "transcript", "id"])
    # Couple expression to all transcripts
    subset_all = expression_df.loc[expression_df["trans_id"].isin(all_trans)]
    subset_target = subset_all.loc[subset_all["trans_id"].isin(target_trans)]
    subset_all = subset_all.set_index("trans_id")
    subset_target = subset_target.set_index("trans_id")
    total_exp_all = subset_all.sum(axis=1)
    # # Couple expression to target transcripts
    total_exp_target = subset_target.sum(axis=1)
    # total_exp_target = subset_target.loc[:, subset_target.columns != "trans_id"].sum()
    if len(total_exp_all) == 0:  # If the gene is not measured, discard the ratio
        return None
    total_all = stats.fmean(list(total_exp_all))
    if total_all == 0.0:  # If there is no expression measured at all, discard the ratio
        return None
    if len(total_exp_target) == 0:  # Target transcript was not measured, return None
        return None

    total_target = stats.fmean(list(total_exp_target))
    ratio = total_target / total_all
    return ratio


def create_full_scatterplot(data, out_dir, complexity):
    """Create a scatterplot of all the data

    :param data: pandas dataframe, containing data to plot
    :param out_dir: str, output directory
    :param complexity: int, complexity used for data creation
    :return: None, figures are created and written to disk
    """
    plt.rcParams.update({"font.size": 20, "figure.figsize": (30, 20)})
    data.plot(x="size", y="ratio", style="o")
    plt.savefig("{}full_scatterplot_complexity{}.png".format(out_dir, complexity))
    plt.close()


def create_zoom_scatterplot(data, out_dir, complexity):
    """Create a zoom in of the full scatterplot

    :param data: pandas dataframe, containing data to plot
    :param out_dir: str, output directory
    :param complexity: int, complexity used for data creation
    :return: None, figures are created and written to disk
    """
    plt.rcParams.update({"font.size": 20, "figure.figsize": (30, 20)})
    data = data.loc[(data["size"] >= 500) & (data["size"] <= 10000)]  # Select certain sizes
    data = data.loc[data["ratio"] <= 4]
    data.plot(x="size", y="ratio", style="o")
    plt.savefig("{}full_zoom_scatterplot_complexity{}.png".format(out_dir, complexity))
    plt.close()


def color_scatterplot_exonincorporation(data, out_dir, bigger_genes, smaller_genes):
    """"""
    data_high = data.loc[data["id"].isin(bigger_genes)].loc[data["size"] < 5000]
    data_low = data.loc[data["id"].isin(smaller_genes)].loc[data["size"] < 5000]
    plt.rcParams.update({"figure.figsize": (19.2, 16.8)})
    fig, ax = plt.subplots()
    ax.scatter(x=data_low["size"], y=data_low["ratio"], s=3, c="Blue", label="Exon usage < 0")
    ax.scatter(x=data_high["size"], y=data_high["ratio"], s=3,  c="Yellow", label="Exon usage >= 0")
    plt.xlim(right=5000)
    plt.ylim(top=3)
    ax.legend()
    plt.savefig("{}/scatterplot_split_on_usage.png".format(out_dir))
    plt.close()


def main():
    """Main function"""
    # Parse arguments
    gff, exp, exp_comp, out_dir, complexity = parse_arg()
    if not os.path.exists("./data/temp/expression.pickle"):
        # Get expression dataframe
        expression_df = acquire_expression_dataframe(exp, exp_comp, out_dir, complexity)
        data = {
            "id": [],
            "size": [],
            "ratio": []
        }
        with open(gff) as file:
            for gene_df in yield_single_gene(file):
                # For each gene determine the expression ratio
                internal_exons = get_internal_exons(gene_df)
                if internal_exons is None:
                    continue  # If there are no internal exons, skip the gene
                for exon in internal_exons:
                    ratio = ratio_expression_exon(gene_df, expression_df, target_exon=exon)
                    if ratio:
                        data["id"].append(exon)  # Number of transcripts in gene
                        data["size"].append(gene_df.loc[gene_df["id"] == exon, "size"].iloc[0])
                        data["ratio"].append(ratio)  # Ratio expression largest exon / total expression
                    else:  # If there is no ratio, skip the exon
                        continue
        data = pd.DataFrame(data)
        data.to_pickle("./data/temp/expression.pickle")
    else:
        data = pd.read_pickle("./data/temp/expression.pickle")

    create_full_scatterplot(data, out_dir, complexity)
    create_zoom_scatterplot(data, out_dir, complexity)

    # Create scatterplot coloured on exon incorporation
    if os.path.exists("./data/out/transcript_ratio/genes_ratio_biggerequal_0.csv") and \
            os.path.exists("./data/out/transcript_ratio/genes_ratio_smaller_0.csv"):
        with open("./data/out/transcript_ratio/genes_ratio_biggerequal_0.csv") as file:
            bigger_genes = file.read().strip().split(",")
        with open("./data/out/transcript_ratio/genes_ratio_smaller_0.csv") as file:
            smaller_genes = file.read().strip().split(",")
        color_scatterplot_exonincorporation(data, out_dir, bigger_genes, smaller_genes)

    # TODO Create boxplots
    return None


if __name__ == '__main__':
    main()
