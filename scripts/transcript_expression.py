#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import argparse as arg
import os
import pandas as pd
import gzip
import statistics as stats
import seaborn as sns
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
    parser.add_argument("smaller_genes", type=str,
                        help="csv file describing genes with a relative exon incorporation lower than 0")
    parser.add_argument("bigger_genes", type=str,
                        help="csv file describing genes with a relative exon incorporation higher than 0")
    parser.add_argument("-c", "--complexity", default=0, type=int, choices=[0, 1, 2],
                        help="determine the tissue depth (0: None (avg. exp per gene), 1: general tissue (e.g. brain), "
                             "2: specific tissue (e.g. brian - frontal cortex)")
    args = parser.parse_args()
    return args.gff3_file, args.expression_file, args.expression_companion, args.smaller_genes, args.bigger_genes,\
           args.out, args.complexity


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
    :return: (float, int), ratio and exit code

    exit code can be deciphered using the function: determine_error_code()
    """
    # Set exit code to 0 (nothing wrong)
    exit_code = 0

    # If no target exon ID is given, determine which exon is the largest and take that ID
    if not target_exon:
        target_exon = determine_target_exon(gene_df, target="largest_internal")
        if not target_exon:
            return None, None
    # Determine in which transcripts this exon is present
    target_trans = list(gene_df.loc[gene_df["id"] == target_exon, "parent"])
    all_trans = list(gene_df.loc[gene_df["type"] == "transcript", "id"])
    # Couple expression to all transcripts
    subset_all = expression_df.loc[expression_df["trans_id"].isin(all_trans)]
    subset_target = subset_all.loc[subset_all["trans_id"].isin(target_trans)]
    subset_all = subset_all.set_index("trans_id")
    subset_target = subset_target.set_index("trans_id")

    if len(subset_all.index) != len(all_trans):
        # Set exit code to 1 (not all transcripts measured within a gene)
        exit_code += 1

    # Determine total and target expression
    total_exp_all = subset_all.sum(axis=1)
    for item in total_exp_all:
        if item == 0.0 and exit_code < 2:
            # Add 2 to exit code (at least one transcript has no expression (expr = 0.0))
            exit_code += 2

    total_exp_target = subset_target.sum(axis=1)
    if len(total_exp_all) == 0:  # If the gene is not measured, discard the ratio
        return None, None
    total_all = stats.fmean(list(total_exp_all))
    if total_all == 0.0:  # If there is no expression measured at all, discard the ratio
        return None, None
    if len(total_exp_target) == 0:
        # Add 4 to exit code (target transcripts were not measured)
        exit_code += 4
        total_target = 0.0
    else:
        total_target = stats.fmean(list(total_exp_target))

    if total_target == 0.0:
        # Add 8 to exit code (target transcripts were not expressed)
        exit_code += 8

    if len(all_trans) == 1:
        # Add 16 to exit code (target exon is from single transcript gene)
        exit_code += 16

    ratio = total_target / total_all
    return ratio, exit_code


def determine_error_code(error_code):
    """Determine the meaning of an error code from ratio_expression_exon().

    :param error_code: int, error code
    :return: list of str, explaination of error code
    """
    codes = {
        1: "not all transcripts measured within the gene",
        2: "at least one transcript has no expression",
        4: "target transcript was not measured",
        8: "target transcripts were not expressed",
        16: "target exon is from single transcript gene",
    }
    keys = list(codes.keys())
    keys.sort(reverse=True)
    errors = {}
    for key in keys:
        if error_code // key > 0:
            error_code -= key
            errors[key] = codes[key]
    return errors


def sns_scatterplot(data, x, y, colour, out_dir, fig_name, x_max):
    """Create scatter plots coloured on error code.

    :param data: pandas dataframe object, data to plot
    :param x: str/index, index to access x values of data
    :param y: str/index, index to access y values of data
    :param colour: str, colour to paint the non-error points
    :param out_dir: str, output directory
    :param fig_name: str, figure name
    :param colour_subgroups: bool, if true, error codes will have different values
    :param x_max: int, max x value to plot
    :return: None, files are created
    """
    plt.rcParams.update({"figure.figsize": (19.2, 16.8), "font.size": 30})
    sns.scatterplot(data=data, x=x, y=y, s=10, color=colour)
    plt.xlim(0, x_max)
    plt.ylim(0.5, 15)
    plt.yscale("log")
    plt.xlabel("Exon size (nt)")
    plt.ylabel("Relative expression")
    plt.savefig("{}/{}.png".format(out_dir, fig_name))
    plt.close()
    return None


def color_scatterplot_exonincorporation(data, out_dir, bigger_genes, smaller_genes):
    """Colour a scatterplot based on the exon_incorperation of the biggest exon and error codes.

    :param data: pandas dataframe object, data to plot
    :param out_dir: str, output directory
    :param bigger_genes: genes with biggest exon incorporation >= 0
    :param smaller_genes: genes with biggest exon incorporation < 0
    :return:
    """
    # Exon usage ratio < 0
    data_low = data.loc[data["id"].isin(smaller_genes)]

    # Only take expression data with no error code
    data_low_sub = data_low.loc[data_low["error"] == 0]
    sns_scatterplot(data_low_sub, "size", "ratio", "Blue", out_dir,
                    "scatterplot_low_usage", x_max=2000)

    # Mean expression per size
    low_mean = data_low.groupby("size").mean()
    sns_scatterplot(low_mean, "size", "ratio", "Blue", out_dir,
                    "scatterplot_low_usage_avg_size", x_max=2000)

    # Exon usage ratio >= 0
    data_high = data.loc[data["id"].isin(bigger_genes)]

    # Only take expression data with no error code
    data_high_sub = data_high.loc[data_high["error"] == 0]
    sns_scatterplot(data_high_sub, "size", "ratio", "Red", out_dir,
                    "scatterplot_high_usage", x_max=2000)

    # Mean expression per size
    high_mean = data_high.groupby("size").mean()
    sns_scatterplot(high_mean, "size", "ratio", "Red", out_dir,
                    "scatterplot_high_usage_avg_size", x_max=2000)


    return None


def create_expression_df(gff, exp, exp_comp, out_dir, complexity):
    """Create an expression dataframe.

    :param gff: str, gff3 file path
    :param exp: str, expression file path
    :param exp_comp: str, expression companion file path
    :param out_dir: str, output directory
    :param complexity: int, [0, 1 or 2], complexity of the tissue (e.g. 1 - brain; 2 - brain, frontal cortext)
    :return: pandas dataframe, expression dataframe
    """
    expression_df = acquire_expression_dataframe(exp, exp_comp, out_dir, complexity)
    data = {
        "id": [],
        "size": [],
        "ratio": [],
        "error": []
    }
    with open(gff) as file:
        for gene_df in yield_single_gene(file):
            # For each gene determine the expression ratio
            internal_exons = get_internal_exons(gene_df)
            if internal_exons is None:
                continue  # If there are no internal exons, skip the gene
            for exon in internal_exons:
                ratio, error_code = ratio_expression_exon(gene_df, expression_df, target_exon=exon)
                if ratio:
                    data["id"].append(exon)  # Number of transcripts in gene
                    data["size"].append(gene_df.loc[gene_df["id"] == exon, "size"].iloc[0])
                    data["ratio"].append(ratio)  # Ratio expression largest exon / total expression
                    data["error"].append(error_code)
                else:  # If there is no ratio, skip the exon
                    continue
    return pd.DataFrame(data)


def main():
    """Main function"""
    # Parse arguments
    gff, exp, exp_comp, smaller_genes, bigger_genes, out_dir, complexity = parse_arg()
    data = create_expression_df(gff, exp, exp_comp, out_dir, complexity)

    # Create scatterplot coloured on exon incorporation
    color_scatterplot_exonincorporation(data, out_dir, bigger_genes, smaller_genes)
    return None


if __name__ == '__main__':
    main()
