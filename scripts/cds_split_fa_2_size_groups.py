#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

import argparse as arg
import pandas as pd


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine codon usage of ")
    parser.add_argument("fasta", type=str,
                        help="fasta file with CDS sequences")
    parser.add_argument("pickle", type=str,
                        help="pandas pickle file, parsed gff3")
    parser.add_argument("output_files", type=str, nargs="+",
                        help="output filenames (incl. group tuple, eg: (0,100) )")
    args = parser.parse_args()

    # Parse filenames into (filename, (start, stop))
    files = []
    for filename in args.output_files:
        if len(filename.split("_")[-1].split(".")[0].split(",")) != 2:
            raise ValueError("Please make sure the size limits are at the end of the filename and separated by a comma."
                             "\nFilename causing error: {}".format(filename))
        start, stop = filename.split("_")[-1].split(".")[0].split(",")
        files.append((filename, (int(start), int(stop))))

    return args.fasta, args.pickle, files


def parse_fasta(file, df):
    """"""
    fa_id = ""
    seq = []

    for idx, line in enumerate(file):
        line = line.strip()
        if line.startswith("#"):
            continue
        elif not line:
            continue
        elif line.startswith(">"):
            if fa_id and seq:
                yield fa_id, "".join(seq)
            line = line[1:].split(":")
            start, stop = line[3][:-3].split("-")
            start = str(int(start) + 1)
            try:
                fa_id = df.loc[(df["chromosome"] == line[2]) & (df["strand"] == line[3][-2]) &
                               (df["start"] <= start) & (df["stop"] >= stop), "id"].iloc[0]
            except IndexError:
                fa_id = ""
            if type(fa_id) != str:
                raise ValueError("Multiple exons map to the same location.")
            seq = []
        else:
            seq.append(line)

    if fa_id and seq:
        yield fa_id, "".join(seq)


def main():
    """Main function."""
    fa_in, pickle, out_files = parse_arguments()
    df = pd.read_pickle(pickle)
    # Only hand internal exons to parse_fasta
    df = df.loc[df["type"] == "exon"]
    df = df.loc[(df["three_prime"] == False) & (df["five_prime"] == False)]

    out = {fname: [] for fname, ranges in out_files}

    with open(fa_in) as file:
        count = 0
        for name, seq in parse_fasta(file, df):
            if count > 20:
                break
            count += 1
            # Determine size of exon
            size = int(df.loc[df["id"] == name, "size"])
            # Add to file of interest
            for fname, ranges in out_files:
                if ranges[0] <= size <= ranges[1]:
                    out[fname].append(">{}\n{}\n".format(name, seq))

    # Output everything to files
    for fname, ranges in out_files:
        with open(fname, "w+") as file:
            file.write("".join(out[fname]))

    return None


if __name__ == "__main__":
    main()