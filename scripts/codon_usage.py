#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import argparse as arg
import pandas as pd
from matplotlib import pyplot as plt


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine codon usage of ")
    parser.add_argument("fasta", type=str,
                        help="fasta file with CDS sequences")
    parser.add_argument("pickle", type=str,
                        help="parsed gff pickle file")
    parser.add_argument("out_dir", type=str,
                        help="output directory (incl. trailing \\")
    args = parser.parse_args()
    if not args.out_dir.endswith("/") or not args.out_dir.endswith("\\"):
        raise ValueError("Please end the out_dir with a trailing \\ or /.")
    return args.fasta, args.pickle, args.out_dir


def parse_fasta(file, df):
    """"""
    seqs = {}
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
                seqs[fa_id] = {
                    "id": fa_id,
                    "seq": "".join(seq)
                }
            line = line[1:].split(":")
            start, stop = line[3][:-3].split("-")
            start = str(int(start) + 1)
            try:
                fa_id = df.loc[(df["chromosome"] == line[2]) & (df["strand"] == line[3][-2]) &
                               (df["start"] <= start) & (df["stop"] >= stop), "id"].iloc[0]
            except IndexError:
                fa_id = ""
            if type(fa_id) != str:
                print(fa_id)
                raise ValueError("Multiple exons map to the same location.")
            seq = []
        else:
            seq.append(line)

    if fa_id and seq:
        seqs[fa_id] = {
            "id": fa_id,
            "seq": "".join(seq)
        }
    return seqs


def codon_aa_table():
    """"""
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codon_table = {
        "codon": codons,
        "amino_acid": [aa for aa in amino_acids]
    }
    return codon_table


def determine_codon_usage(sequences):
    """"""
    codons = {}

    for fa_id, fasta in sequences.items():
        len_seq = len(fasta["seq"])
        for i in range(0, len_seq, 3):
            if i + 3 > len_seq:
                # There is no full codon left in the sequence
                break
            codon = fasta["seq"][i:i+3]
            if codons.get(codon) is None:
                codons[codon] = 1
            else:
                codons[codon] += 1

    return codons


def codon_usage_df(sequences):
    """"""
    codon_usage = determine_codon_usage(sequences)
    total = sum(codon_usage.values())

    df = pd.DataFrame(codon_aa_table())
    df2 = pd.DataFrame(codon_usage, index=["count"])
    df2 = df2.transpose(copy=False)
    df2 = df2.reset_index()
    df2 = df2.rename(columns={"index": "codon"})
    df = df.merge(df2, how="left")
    df = df.fillna(0.0)
    df["fraction"] = df["count"] / total

    return df


def hbar_codon_usage(df, filename):
    """"""
    plt.rcParams.update({"font.size": 8, "figure.figsize": (19.2, 10.8)})
    df.plot.barh(x="codon", y="fraction")
    plt.savefig(filename)
    plt.close()


def main():
    """"""
    fasta, pickle, out_dir = parse_arguments()
    df = pd.read_pickle(pickle)
    # parse sequences
    with open(fasta) as file:
        seqs = parse_fasta(file, df)

    # divide into groups
    groups = {"size_<50": {key: seqs[key] for key in seqs.keys() & list(df.loc[df["size"] < 50, "id"])},
              "size_50-300": {key: seqs[key] for key in seqs.keys() &
                              list(df.loc[(df["size"] >= 50) & (df["size"] <= 300), "id"])},
              "size_300>": {key: seqs[key] for key in seqs.keys() & list(df.loc[df["size"] > 300, "id"])}}

    # determine codon usage per group
    data = {}
    for idx, group in groups.items():
        codon_usage = codon_usage_df(group)
        hbar_codon_usage(codon_usage, "{}hbar_graph_{}.png".format(out_dir, idx))

    return None


if __name__ == "__main__":
    main()
