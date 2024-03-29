#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

import pandas as pd
import argparse as arg
from helpers import yield_gff_line, determine_gff


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Convert exons to a bed file.")
    parser.add_argument("gff3", type=str,
                        help="gff3 file to parse")
    parser.add_argument("out_file", type=str,
                        help="output filename")
    parser.add_argument("-p", "--phase", action="store_true")
    args = parser.parse_args()
    return args.gff3, args.out_file, args.phase


def write_to_bed(df, chrom, start, stop, strand, parent, out_file):
    """Write a gff3 dataframe to a bed file.

    :param df: pandas dataframe
    :param chrom: str/index, index to access chromosome information from df
    :param start: str/index, index to access start information from df
    :param stop: str/index, index to access stop information from df
    :param strand: str/index, index to access strand information from df
    :param parent: str/index, index to access parent information from df
    :param out_file: str, output file
    :return: None, file is created
    """
    df["score"] = 0
    columns = [chrom, start, stop, parent, "score", strand]
    df[columns].to_csv(out_file, sep="\t", index=False, header=False)
    return None


def parse_cds_from_gff(file, phase_aware=False):
    """Parse all the CDS regions from a gff3 file.

    :param file: iterable, open file object
    :param phase_aware: bool, move the CDS start to the first codon in the exon (default: False)
    :return: pandas dataframe object, CDS information
    """
    data = {
        "chromosome": [],
        "start": [],
        "stop": [],
        "strand": [],
        "parent": []
    }

    for line in yield_gff_line(file):
        if line[2] == "CDS":
            data["chromosome"].append(line[0])
            if phase_aware:
                phase = int(line[7])
                data["start"].append(int(line[3]) + phase)
            else:
                data["start"].append(line[3])
            data["stop"].append(line[4])
            data["strand"].append(line[6])
            data["parent"].append(determine_gff(line[8], "Parent=transcript:"))
    df = pd.DataFrame(data)
    return df


def main():
    """Main function."""
    gff, out_file, phase = parse_arguments()

    with open(gff) as file:
        df = parse_cds_from_gff(file, phase)

    # BedTools sees the stop index as non-inclusive, but gff3 sees it as inclusive
    df["start"] = df["start"].apply(pd.to_numeric)
    df["start"] -= 1
    write_to_bed(df, "chromosome", "start", "stop", "strand", "parent", out_file)
    return None


if __name__ == "__main__":
    main()
