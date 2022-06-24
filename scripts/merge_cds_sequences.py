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
    parser.add_argument("bed", type=str,
                        help="bed file to parse")
    parser.add_argument("out_file", type=str,
                        help="output filename")
    args = parser.parse_args()
    return args.bed, args.out_file


def parse_cds_into_full_seq(bed):
    """Parse a bed file containing coding sequences.

    :param bed: iterable, open bed file
    :return: dict, protein_id & sequence, yields one fasta sequence at a time
    """
    data = {
        "target": "",
        "seq": []
    }

    for line in bed:
        line = line.strip().split("\t")
        if data["target"] != line[3]:
            if data["seq"]:
                data["seq"] = "".join(data["seq"])
                yield data
            data = {
                "target": line[3],
                "seq": [line[6]]
            }
        else:
            if line[5] == "-":
                data["seq"].insert(0, line[6])
            else:
                data["seq"].append(line[6])
    if data["seq"]:
        data["seq"] = "".join(data["seq"])
        yield data


def main():
    """Main function."""
    bed, out_file = parse_arguments()
    with open(bed) as file, open(out_file, "w+") as f_out:
        for fasta in parse_cds_into_full_seq(file):
            f_out.write(">{}\n{}\n\n".format(fasta["target"], fasta["seq"]))
    return None


if __name__ == "__main__":
    main()
