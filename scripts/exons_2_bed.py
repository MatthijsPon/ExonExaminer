#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import pandas as pd


def write_to_bed(df, chrom, start, stop, id, out_file):
    """Write a gff3 dataframe to a bed file.

    :param df: pandas dataframe
    :param chrom: str/index, index to access chromosome information from df
    :param start: str/index, index to access start information from df
    :param stop: str/index, index to access stop information from df
    :param id: str/index, index to access ID information from df
    :param out_file: str, output file
    :return: None, file is created
    """
    columns = [chrom, start, stop, id]
    df[columns].to_csv(out_file, sep="\t", index=False)
    return None


def main():
    """Main function."""
    in_file = "data/temp/parsed_human_gff.pickle"
    out_file = "data/out/exons_2_bed/human_exons.bed"
    df = pd.read_pickle(in_file)
    df = df.loc[df["type"] == "exon"]
    # BedTools sees the stop index as non-inclusive, but gff3 sees it as inclusive
    df["stop"] += 1
    write_to_bed(df, "chromosome", "start", "stop", "id", out_file)
    return None


if __name__ == "__main__":
    main()
