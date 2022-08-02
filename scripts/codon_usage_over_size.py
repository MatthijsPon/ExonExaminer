#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Analyse the codon usage over size
"""
import argparse as arg
import pandas as pd
from matplotlib import pyplot as plt
from cds_split_fa_2_size_groups import yield_fasta_record
import seaborn as sns


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine the difference in aa and codon usage of two fasta files.")
    parser.add_argument("fasta", type=str,
                        help="fasta file with CDS sequences")
    parser.add_argument("pickle", type=str,
                        help="pandas pickle file containing exon information")
    parser.add_argument("out_dir", type=str,
                        help="output directory (incl. trailing \\")
    args = parser.parse_args()
    if not args.out_dir.endswith("/") and not args.out_dir.endswith("\\"):
        raise ValueError("Please end the out_dir with a trailing \\ or /.")
    return args.fasta, args.pickle, args.out_dir


def codon_aa_table():
    """Create an amino acid dataframe, containing codons, amino acids and aa type.

    :return: pandas dataframe object, AA dataframe
    """
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    aa_type = {
        "*": "Stop",
        "R": "Positively charged",
        "H": "Positively charged",
        "K": "Positively charged",
        "D": "Negatively charged",
        "E": "Negatively charged",
        "S": "Polar uncharged",
        "T": "Polar uncharged",
        "N": "Polar uncharged",
        "Q": "Polar uncharged",
        "P": "Polar uncharged",
        "C": "Polar uncharged",
        "G": "Non-polar, aliphatic",
        "A": "Non-polar, aliphatic",
        "V": "Non-polar, aliphatic",
        "I": "Non-polar, aliphatic",
        "L": "Non-polar, aliphatic",
        "M": "Non-polar, aliphatic",
        "F": "Aromatic",
        "Y": "Aromatic",
        "W": "Aromatic"
    }
    codon_table = {
        "codon": codons,
        "amino_acid": [aa for aa in amino_acids],
        "aa_type": [aa_type[aa] for aa in amino_acids]
    }
    return pd.DataFrame(codon_table)


def fasta_2_dict(file):
    """Convert a fasta file to a dictionary with ID as key and seq as value.

    :param file: iterable, open file object
    :return: dict, fasta dictionary
    """
    fasta_dict = {}
    for fa_id, fasta in yield_fasta_record(file):
        fasta_dict[fa_id] = fasta

    return fasta_dict


def determine_codon_usage(fasta_dict):
    """Determine the usage of codons from all sequences in a fasta file

    :param file: iterable, open fasta file object
    :return: dict, str: int, codon: count
    """
    codons = {}

    for fa_id, fasta in fasta_dict.items():
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


def cu_2_fraction(dict):
    total = sum(dict.values())
    for key, value in dict.items():
        dict[key] = value / total
    return dict


def amino_acid_usage(dict):
    """"""
    df = codon_aa_table()
    df2 = pd.DataFrame(dict, index=["count"])
    df2 = df2.transpose(copy=False)
    df2 = df2.reset_index()
    df2 = df2.rename(columns={"index": "codon"})
    df = df.merge(df2, how="left")
    df = df.groupby(["amino_acid", "aa_type"]).sum()
    df = df.reset_index()
    df = df.rename(columns={"count": "fraction"})
    df = df.set_index(["amino_acid", "aa_type"])
    return df


def main():
    """Main function."""
    fasta, pickle, out_dir = parse_arguments()

    exon_df = pd.read_pickle(pickle)
    # Only take along internal exons
    exon_df = exon_df.loc[(exon_df["three_prime"] == False) &
                          (exon_df["five_prime"] == False)]
    with open(fasta) as file:
        fasta_dict = fasta_2_dict(file)

    step = 25
    usage_dict = {}
    aa_dfs = []
    for i in range(0, 2000, step):
        yf_exons = list(exon_df.loc[(exon_df["size"] >= i) & (exon_df["size"] < i + step), "id"])
        yf_fasta = {k: fasta_dict[k] for k in fasta_dict.keys() & yf_exons}
        yf_cu = determine_codon_usage(yf_fasta)
        usage_dict[i] = cu_2_fraction(yf_cu)
        aa_df_temp = amino_acid_usage(yf_cu)
        aa_df_temp = aa_df_temp.rename(columns={"fraction": i})
        aa_dfs.append(aa_df_temp)

    # Ignore codon usage, as it is too cluttered.
    # usage_df = pd.DataFrame(usage_dict)
    # usage_df = usage_df.transpose()
    # usage_df = usage_df.reset_index()
    #
    # plt.rcParams.update({"figure.figsize": (19.2, 16.8), "font.size": 30})
    # sns.lineplot(x="index", y="value", data=pd.melt(usage_df, ["index"]), hue="variable")
    # plt.setp(plt.legend().get_texts(), fontsize='20')  # Set legend font size
    # plt.title("All codons")
    # plt.xlabel("Exon size ({} nt wide buckets)".format(step))
    # plt.ylabel("Codon usage (fraction)")
    # plt.savefig("{}/codons_over_size_all_types.png".format(out_dir))
    # plt.close()

    aa_df = pd.concat(aa_dfs, axis=1)
    aa_df = aa_df.drop("*")
    # Plot each amino acid type separately
    for idx, group in aa_df.groupby("aa_type"):
        group = group.transpose()
        group = group.reset_index()
        plt.rcParams.update({"figure.figsize": (19.2, 16.8), "font.size": 30})
        sns.lineplot(x="index", y="value", data=pd.melt(group, ["index"]), hue="amino_acid")
        plt.title(idx)
        plt.setp(plt.legend().get_texts(), fontsize='20')  # Set legend font size

        plt.xlabel("Exon length ({} nt wide buckets)".format(step))
        plt.ylabel("Amino acid usage (fraction)")
        plt.savefig("{}/aa_over_size_{}.png".format(out_dir, idx))
        plt.close()

    aa_df = aa_df.transpose()
    aa_df = aa_df.reset_index()
    # Plot all amino acid types together
    plt.rcParams.update({"figure.figsize": (19.2, 16.8), "font.size": 30})
    sns.lineplot(x="index", y="value", data=pd.melt(aa_df, ["index"]), hue="amino_acid")
    plt.setp(plt.legend().get_texts(), fontsize='20')  # Set legend font size

    plt.title("All types")
    plt.xlabel("Exon length ({} nt wide buckets)".format(step))
    plt.ylabel("Amino acid usage (fraction)")
    plt.savefig("{}/aa_over_size_all_types.png".format(out_dir))
    plt.close()
    return None


if __name__ == "__main__":
    main()
