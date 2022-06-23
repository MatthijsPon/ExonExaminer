#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Convert the exons of a pandas dataframe (pickle) to a bed file.
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
    args = parser.parse_args()
    return args.gff3, args.out_file


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
    columns = [chrom, start, stop, strand, parent]
    df[columns].to_csv(out_file, sep="\t", index=False, header=False)
    return None


def parse_cds_from_gff(file):
    """Parse all the CDS regions from a gff3 file.

    :param file: iterable, open file object
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
            data["start"].append(line[3])
            data["stop"].append(line[4])
            data["strand"].append(line[6])
            data["parent"].append(determine_gff(line[8], "Parent=transcript:"))
    df = pd.DataFrame(data)
    return df


def main():
    """Main function."""
    gff, out_file = parse_arguments()

    with open(gff) as file:
        df = parse_cds_from_gff(file)

    # BedTools sees the stop index as non-inclusive, but gff3 sees it as inclusive
    df["stop"] = df["stop"].apply(pd.to_numeric)
    df["stop"] += 1
    write_to_bed(df, "chromosome", "start", "stop", "strand", "parent", out_file)
    return None


if __name__ == "__main__":
    main()
