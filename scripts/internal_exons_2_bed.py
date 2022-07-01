#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Convert the exons of a pandas dataframe (pickle) to a bed file.
"""
import pandas as pd
import argparse as arg


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Convert exons to a bed file.")
    parser.add_argument("pickle", type=str,
                        help="Pandas pickle (saved gff3 dataframe)")
    parser.add_argument("out_file", type=str,
                        help="output filename")
    args = parser.parse_args()
    return args.pickle, args.out_file


def write_to_bed(df, chrom, start, stop, item_id, out_file):
    """Write a gff3 dataframe to a bed file.

    :param df: pandas dataframe
    :param chrom: str/index, index to access chromosome information from df
    :param start: str/index, index to access start information from df
    :param stop: str/index, index to access stop information from df
    :param item_id: str/index, index to access ID information from df
    :param out_file: str, output file
    :return: None, file is created
    """
    columns = [chrom, start, stop, item_id]
    df[columns].to_csv(out_file, sep="\t", index=False, header=False)
    return None


def main():
    """Main function."""
    pickle, out_file = parse_arguments()

    df = pd.read_pickle(pickle)
    df = df.loc[df["type"] == "exon"]
    # Only take along internal exons
    df = df.loc[(df["three_prime"] == False) & (df["five_prime"] == False)]
    # BedTools sees the stop index as non-inclusive, but gff3 sees it as inclusive
    df["start"] = df["start"].apply(pd.to_numeric)
    df["start"] -= 1
    write_to_bed(df, "chromosome", "start", "stop", "id", out_file)
    return None


if __name__ == "__main__":
    main()
