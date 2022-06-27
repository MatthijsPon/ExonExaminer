#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import argparse as arg
import statistics

import pandas as pd

from helpers import yield_single_gene, get_internal_exons, determine_target_exon


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Parse a gff3 file, calculate exon "
                                            "incorporation, and save as pandas pickle.")
    parser.add_argument("gff3_file", type=str,
                        help="gff3 file to parse")
    parser.add_argument("out_file", type=str,
                        help="output filename")
    args = parser.parse_args()
    return args.gff3_file, args.out_file


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
    target_trans = exon_transcript_presence(gene_df, target_exon)
    target_ratio = float(target_trans) / all_trans
    # Get all internal exons
    internal_exons = get_internal_exons(gene_df)
    if not internal_exons:
        return None
    internal_ratios = []
    for internal_exon in internal_exons:
        int_ex_pres = exon_transcript_presence(gene_df, internal_exon)
        internal_ratios.append(float(int_ex_pres) / all_trans)
    internal_ratio = statistics.fmean(internal_ratios)
    return target_ratio - internal_ratio


def exon_incorporation(infile, outfile):
    """"""
    ratio_dict = {
        "id": [],
        "exon_size": [],
        "ratio": [],
        "no_trans": [],
        "parent_gene": []
    }

    for gene in yield_single_gene(infile):
        internal_exons = get_internal_exons(gene)
        if internal_exons is None:
            continue  # If there are no internal exons, skip the ratios
        trans_amount = int(gene.loc[gene["type"] == "transcript", "id"].nunique())
        for exon in internal_exons:
            ratio = ratio_transcripts_exon(gene, target_exon=exon)
            if ratio is not None:
                ratio_dict["id"].append(exon)
                ratio_dict["ratio"].append(ratio)
                size = gene.loc[gene["id"] == exon, "size"].iloc[0]
                ratio_dict["exon_size"].append(int(size))
                ratio_dict["no_trans"].append(trans_amount)
                ratio_dict["parent_gene"].append(gene.loc[gene["type"] == "gene", "id"].iloc[0])
    # Create DF from exon size and ratio
    results = pd.DataFrame(ratio_dict)
    results.to_pickle(outfile)
    return None


def main():
    """Main function."""
    gff_file, output_file = parse_arguments()

    with open(gff_file) as file:
        exon_incorporation(file, output_file)
    return None


if __name__ == "__main__":
    main()
