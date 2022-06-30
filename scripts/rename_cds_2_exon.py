#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import pandas as pd
import argparse as arg


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine codon usage of ")
    parser.add_argument("fasta", type=str,
                        help="fasta file with CDS sequences")
    parser.add_argument("pickle", type=str,
                        help="pandas pickle file, parsed gff3")
    parser.add_argument("output", type=str,
                        help="output filename")
    args = parser.parse_args()

    return args.fasta, args.pickle, args.output


def parse_fasta(file, df):
    """"""
    fa_id = ""
    seq = []

    for idx, line in enumerate(file):
        if idx % 100 == 0:
            print("Line {} ({:.3f}%)".format(idx, idx / 1702848 * 100))
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
    fa_in, pickle, fa_out = parse_arguments()

    df = pd.read_pickle(pickle)
    # Only hand internal exons to parse_fasta
    df = df.loc[df["type"] == "exon"]
    df = df.loc[(df["three_prime"] == False) & (df["five_prime"] == False)]

    with open(fa_in) as file, open(fa_out, "w+") as out_f:
        for exon_id, seq in parse_fasta(file, df):
            out_f.write(">{}\n{}\n".format(exon_id, seq))

    return None


if __name__ == "__main__":
    main()
