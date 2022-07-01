#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import argparse as arg
import pandas as pd
from matplotlib import pyplot as plt
from cds_split_fa_2_size_groups import yield_fasta_record


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine codon usage of ")
    parser.add_argument("fasta", type=str,
                        help="fasta file with CDS sequences")
    parser.add_argument("out_dir", type=str,
                        help="output directory (incl. trailing \\")
    parser.add_argument("group_name", type=str,
                        help="group name to use in filenames")
    args = parser.parse_args()
    if not args.out_dir.endswith("/") and not args.out_dir.endswith("\\"):
        raise ValueError("Please end the out_dir with a trailing \\ or /.")
    return args.fasta, args.out_dir, args.group_name


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


def determine_codon_usage(file):
    """"""
    codons = {}

    for fa_id, fasta in yield_fasta_record(file):
        len_seq = len(fasta)
        for i in range(0, len_seq, 3):
            if i + 3 > len_seq:
                # There is no full codon left in the sequence
                break
            codon = fasta[i:i+3]
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
    fasta, out_dir, group = parse_arguments()

    with open(fasta) as file:
        cu_dict = determine_codon_usage(file)
    cu_df = codon_usage_df(cu_dict)
    hbar_codon_usage(cu_df, "{}hbar_graph_{}.png".format(out_dir, group))

    return None


if __name__ == "__main__":
    main()
