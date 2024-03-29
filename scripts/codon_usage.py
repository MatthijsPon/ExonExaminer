#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Analyse the codon usage difference between two fasta files
"""
import argparse as arg
import pandas as pd
from matplotlib import pyplot as plt
from cds_split_fa_2_size_groups import yield_fasta_record
import seaborn as sns
import statistics


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine the difference in aa and codon usage of two fasta files.")
    parser.add_argument("fasta1", type=str,
                        help="fasta file with CDS sequences")
    parser.add_argument("fasta2", type=str,
                        help="fasta file with CDS sequences to compare to fasta1")
    parser.add_argument("out_dir", type=str,
                        help="output directory (incl. trailing \\")
    parser.add_argument("group_name", type=str,
                        help="group name to use in filenames")
    args = parser.parse_args()
    if not args.out_dir.endswith("/") and not args.out_dir.endswith("\\"):
        raise ValueError("Please end the out_dir with a trailing \\ or /.")
    return args.fasta1, args.fasta2, args.out_dir, args.group_name


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


def determine_codon_usage(file):
    """Determine the usage of codons from all sequences in a fasta file

    :param file: iterable, open fasta file object
    :return: dict, str: int, codon: count
    """
    codons_frac = {}
    bases = "TCAG"
    codon_options = [a + b + c for a in bases for b in bases for c in bases]

    for fa_id, fasta in yield_fasta_record(file):
        codons = {}
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
        total = sum(codons.values())
        for codon in codon_options:
            if codons.get(codon) is None:
                if codons_frac.get(codon) is None:
                    codons_frac[codon] = [0.0]
                else:
                    codons_frac[codon].append(0.0)
            else:
                if codons_frac.get(codon) is None:
                    codons_frac[codon] = [float(codons.get(codon) / total)]
                else:
                    codons_frac[codon].append(float(codons.get(codon) / total))

    for codon, frac_list in codons_frac.items():
        codons_frac[codon] = statistics.fmean(frac_list)
    return codons_frac


def codon_usage_df(codon_usage):
    """Create a codon usage dataframe, combining codon usage and codon_aa_table()

    :param codon_usage: dict, codon usage counts
    :return: pandas dataframe object, codon usage dataframe including fractions
    """
    total = sum(codon_usage.values())

    df = codon_aa_table()
    df2 = pd.DataFrame(codon_usage, index=["count"])
    df2 = df2.transpose(copy=False)
    df2 = df2.reset_index()
    df2 = df2.rename(columns={"index": "codon"})
    df = df.merge(df2, how="left")
    df = df.fillna(0.0)
    df["fraction"] = df["count"] / total

    return df


def amino_acid_usage(cu_dict):
    """Calculate amino acid usage.

    :param cu_dict: dict, codon usage dictionary
    :return: amino acids usage dataframe
    """

    df = codon_aa_table()
    df2 = pd.DataFrame(cu_dict, index=["fraction"])
    df2 = df2.transpose(copy=False)
    df2 = df2.reset_index()
    df2 = df2.rename(columns={"index": "codon"})
    df = df.merge(df2, how="left")
    df = df.groupby(["amino_acid", "aa_type"]).sum()
    df = df.reset_index()

    return df


def hbar_aa_usage(df, filename):
    """"""
    df = df[df.amino_acid != "*"]
    df = df.sort_values(["amino_acid"])
    aa_order = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "P", "G",
                "A", "V", "I", "L", "M", "F", "Y", "W"]
    df.amino_acid = df.amino_acid.astype("category")
    df.amino_acid = df.amino_acid.cat.set_categories(aa_order)
    df = df.sort_values(["amino_acid"])
    plt.rcParams.update({"font.size": 22, "figure.figsize": (19.2, 10.8), "legend.loc": "lower right"})
    ax = sns.barplot(x="diff_perc", y="amino_acid", data=df, hue="aa_type", dodge=False)
    plt.xlim(-50, 50)
    plt.xlabel("Difference in amino acid usage (%)")
    plt.ylabel("Amino acid")
    plt.legend(title="Amino acid side chain type")
    plt.savefig(filename)
    plt.close()
    return None


def hbar_codon_usage(df, filename):
    """Create a horizontal codon usage bar graph, coloured on AA type

    :param df: pandas dataframe object, containing data to graph
    :param filename: str, output filename
    :return: None, figures are created.
    """
    aa_order = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "P", "G",
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
    df.amino_acid = df.amino_acid.astype("category")
    df.amino_acid = df.amino_acid.cat.set_categories(aa_order)
    df = df.sort_values(["amino_acid"])
    plt.rcParams.update({"font.size": 16, "figure.figsize": (19.2, 10.8), "legend.loc": "lower left"})
    ax = sns.barplot(x="diff_frac", y="codon", data=df, hue="aa_type", dodge=False)
    plt.savefig(filename)
    plt.close()
    return None


def main():
    """Main function."""
    fasta1, fasta2, out_dir, group = parse_arguments()

    with open(fasta1) as file:
        cu_dict1 = determine_codon_usage(file)
    with open(fasta2) as file:
        cu_dict2 = determine_codon_usage(file)

    # Amino acid usage
    au_df1 = amino_acid_usage(cu_dict1)
    au_df2 = amino_acid_usage(cu_dict2)
    au_df = pd.merge(au_df1, au_df2, on=["amino_acid", "aa_type"], how="inner")
    au_df["diff_frac"] = au_df["fraction_y"] - au_df["fraction_x"]
    au_df["diff_perc"] = (au_df["fraction_y"] - au_df["fraction_x"]) / au_df["fraction_x"] * 100
    hbar_aa_usage(au_df, "{}aa_diff_hbar_{}.png".format(out_dir, group))
    # Codon usage
    cu_df1 = codon_usage_df(cu_dict1)
    cu_df2 = codon_usage_df(cu_dict2)
    cu_df = pd.merge(cu_df1, cu_df2, on=["codon", "amino_acid", "aa_type"], how="inner")
    cu_df["diff_frac"] = cu_df["fraction_y"] - cu_df["fraction_x"]

    hbar_codon_usage(cu_df, "{}codon_usage_hbar_{}.png".format(out_dir, group))

    return None


if __name__ == "__main__":
    main()
