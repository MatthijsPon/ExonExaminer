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

from parse import yield_gff_line


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
                "parent": [],
                "size": [],
            }
            gene_id = determine_gff(line[8], "ID=gene:")
            transcript = []
            exon = []
            # Add data to dict
            gene_dict["type"].append("gene")
            gene_dict["id"].append(gene_id)
            gene_dict["parent"].append(None)
            gene_dict["size"].append(None)
        elif line[2] in ["transcript", "mRNA", "primary_transcript"] and \
                determine_gff(line[8], "biotype=") == "protein_coding":
            if determine_gff(line[8], "Parent=gene:") == gene_id:
                trans_id = determine_gff(line[8], "ID=transcript:")
                transcript.append(trans_id)
                # Add data to dict
                gene_dict["type"].append("transcript")
                gene_dict["id"].append(trans_id)
                gene_dict["parent"].append(gene_id)
                gene_dict["size"].append(None)
        elif line[2] == "exon":
            if determine_gff(line[8], "Parent=transcript:") in transcript:
                exon_id = determine_gff(line[8], "Name=")
                exon.append(exon_id)
                gene_dict["type"].append("exon")
                gene_dict["id"].append(exon_id)
                gene_dict["parent"].append(determine_gff(line[8], "Parent=transcript:"))
                gene_dict["size"].append((int(line[4]) - int(line[3])))

    if gene_id and transcript and exon:
        yield pd.DataFrame.from_dict(gene_dict)


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
    with open("out.txt", "w+") as f:
        f.write(", ".join(list(companion_df["SAMPID"])))
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


def determine_target_exon(gene_df, target):
    """Determine which exon to target

    :param gene_df: pandas dataframe, gene dataframe
    :param target: str, "largest" or "largest_internal" to decide target
    :return: str, geneID if found, otherwise None
    """
    if target == "largest":
        return gene_df[gene_df["size"] == gene_df["size"].max()]["id"].iloc[0]  # Select largest exon
    elif target == "largest_internal":
        non_prime = set()
        only_exons = gene_df[gene_df["type"] == "exon"]
        # Select exons per transcript and apply True only if they are prime
        for idx, group in only_exons.groupby("parent"):
            exons_transcript = list(group["id"])
            exons_transcript.pop(0)  # Remove the first exon
            if len(exons_transcript) >= 1:  # Cannot pop twice if there is only one exon
                exons_transcript.pop(-1)  # Remove the last exon
            # All exons that have a non-prime function at least once are possible targets
            for exon in exons_transcript:
                non_prime.add(exon)
        if len(non_prime) >= 1:  # If there are internal exons return the largest, otherwise return none.
            subset_internal = gene_df[gene_df["id"].isin(non_prime)]  # Select only internal exons from gene_df
            # Select largest exon from the subset
            largest = subset_internal[subset_internal["size"] == subset_internal["size"].max()]["id"].iloc[0]
            return largest
        return None
    else:
        raise ValueError("Incorrect target given in determine_target_exon().")


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
    if len(total_exp_target) == 0:  # Target transcript was not measured, skip
        return None
    total_target = stats.fmean(list(total_exp_target))

    if total_all == 0.0:  # If there is no expression measured at all, discard the ratio
        return None
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
    data.plot(x="max_exon", y="ratio", style="o")
    plt.savefig("{}full_scatterplot_complexity{}.png".format(out_dir, complexity))
    plt.close()


def create_scatterplot_no_transcripts(data, out_dir, complexity):
    """Create a scatterplot split on No. of transcripts in gene

    :param data: pandas dataframe, containing data to plot
    :param out_dir: str, output directory
    :param complexity: int, complexity used for data creation
    :return: None, figures are created and written to disk
    """
    trans_numbers = data["no_trans"].unique()
    trans_numbers.sort()
    closest_root = math.ceil(math.sqrt(len(trans_numbers)))  # Determine minimum square size needed to plot all graphs
    # Plot stuff
    plt.rcParams.update({"font.size": 20, "figure.figsize": (30, 20)})
    # Plot each image separately
    for key in trans_numbers:
        data.loc[data["no_trans"] == key].plot(x="max_exon", y="ratio", style="o",
                                               title="{} transcripts in gene".format(key))
        plt.savefig("{}scatterplot_complexity{}_{}trans.png".format(out_dir, complexity, key))
        plt.close()

    """ Putting all plots in a single image
    fig, axes = plt.subplots(closest_root, closest_root)
    for idx, key in enumerate(trans_numbers):
        row = idx // closest_root  # row index
        col = idx % closest_root  # col index
        print(idx, col, key)
        data.loc[data["no_trans"] == key].plot(x="max_exon", y="ratio", style="o", ax=axes[(row, col)],
                                                title="{} transcripts in gene".format(key))
    """


def main():
    """Main function"""
    start_time = time.perf_counter()
    # Parse arguments
    gff, exp, exp_comp, out_dir, complexity = parse_arg()
    # Get expression dataframe
    expression_df = acquire_expression_dataframe(exp, exp_comp, out_dir, complexity)
    data = {
        "no_trans": [],
        "max_exon": [],
        "ratio": []
    }
    with open(gff) as file:
        for gene_df in yield_single_gene(file):
            # For each gene determine the expression ratio
            ratio = ratio_expression_exon(gene_df, expression_df)
            if ratio:
                data["no_trans"].append((gene_df.type.values == "transcript").sum())  # Number of transcripts in gene
                data["max_exon"].append(gene_df["size"].max())  # Size of the largest exon
                data["ratio"].append(ratio)  # Ratio expression largest exon / total expression
            else:  # If there is no ratio, skip the gene
                continue
    data = pd.DataFrame(data)
    create_full_scatterplot(data, out_dir, complexity)
    create_scatterplot_no_transcripts(data, out_dir, complexity)
    print("Done, time taken: {:.2f} seconds\n".format(time.perf_counter() - start_time))


if __name__ == '__main__':
    main()
