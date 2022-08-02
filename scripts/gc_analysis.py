#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import argparse as arg
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """

    parser = arg.ArgumentParser(description="Analyse the gc contents of exons.")
    parser.add_argument("pickle", type=str,
                        help="exon incorporation pickle file")
    parser.add_argument("bed_file", type=str,
                        help="bed file contraining fasta sequences of exons")
    parser.add_argument("quartile_file", type=str,
                        help="text file containing species 5%/95% quartile information")
    parser.add_argument("inc_border_file", type=str,
                        help="text file containing the borders of the positive incorporation region")
    parser.add_argument("out_dir", type=str,
                        help="output directory")

    args = parser.parse_args()
    return args.pickle, args.bed_file, args.quartile_file, \
           args.inc_border_file, args.out_dir


def average_values(df, x, y, minimal=None, maximal=None):
    """Calculate the average values of a dataframe column

    :param df: pandas dataframe object
    :param x: column index, column to apply min and max values to
    :param y: column index, column to average
    :param minimal: int, minimal value of x (default: None)
    :param maximal: int, maximal value of x (default: None)
    :return: float, mean of column y
    """
    if minimal is not None and maximal is not None:
        mean = df.loc[(df[x] >= minimal) & (df[x] <= maximal), y].mean()
    elif minimal is not None:
        mean = df.loc[df[x] >= minimal, y].mean()
    elif maximal is not None:
        mean = df.loc[df[x] <= maximal, y].mean()
    else:
        mean = df[y].mean()
    return mean


def scatterplot_roll_avg(data, x, y, window, out_file,
                         Q5=49, Q95=288, inc_bottom=69, inc_top=199,
                         x_min=None, x_max=None, hue=None):
    """Create a scatter plot with a moving average.

    :param data: 2d object (e.g. pandas dataframe) with data to plot
    :param x: str/index, index to access x values of data
    :param y: str/index, index to access y values of data
    :param window: int, no. of data points to use in the rolling average
    :param out_file: str, output filename
    :param Q5: int, 5% quartile size
    :param Q95: int, 95% quartile size
    :param inc_bottom: int, size corresponding to the bottom of the cumulative bar plot
    :param inc_top: int, size corresponding to the top of the cumulative bar plot
    :param x_min: int, limit of the left side of the x scale
    :param x_max: int, limit of the right side of the x scale
    :param hue: str/index, index to access hue values of data
    :return: None, files are created
    """
    plt.rcParams.update({"font.size": 16, "figure.figsize": (19.2, 10.8)})
    if x_max:
        plt.xlim(right=x_max)
    elif x_min:
        plt.xlim(left=x_min)
    elif x_min and x_max:
        plt.xlim(x_min, x_max)
    if hue:
        # Linegraph with continuous colour palette
        cmap = sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)
        points = plt.scatter(x=data[x], y=data[y], c=data[hue], cmap=cmap)
        cbar = plt.colorbar(points)
        cbar.set_label("Exon incorporation ratio", rotation=270)
    else:
        # Linegraph with continuous colour palette
        sns.set_style("whitegrid")
        sns.scatterplot(x=x, y=y, data=data, color="b")

    # Quartile lines
    plt.axvline(x=Q5, color="r", linestyle=":", label="Q05 ({} nt)".format(Q5))
    plt.axvline(x=Q95, color="g", linestyle=":", label="Q95 ({} nt)".format(Q95))
    plt.axvline(x=inc_bottom, color="y", linestyle=":", label="Positive incorporation start ({} nt)".format(inc_bottom))
    plt.axvline(x=inc_top, color="purple", linestyle=":", label="Positive incorporation end ({} nt)".format(inc_top))

    # rolling average line
    data['MA_group'] = data[y].rolling(window=window).mean()
    plt.plot(data[x], data["MA_group"], color="r", label="Rolling average of {} data points".format(window))

    plt.xlabel("Exon length (nt)")
    plt.ylabel("Average GC content (%)")
    plt.legend()
    plt.savefig(out_file)
    plt.close()
    return None


def analyse_gc(gc_df, inc_df, Q5, Q95, inc_border_low, inc_border_high, outdir):
    """"""

    # Get mean GC content per size
    gc_mean = gc_df.groupby("seq_len").mean()
    gc_mean = gc_mean.reset_index()
    gc_mean["gc_perc"] = gc_mean["gc_perc"] * 100

    scatterplot_roll_avg(gc_mean, "seq_len", "gc_perc", 20, "{}scatter_roll_avg_20_by_size.png".format(outdir),
                         Q5, Q95, inc_border_low, inc_border_high, x_max=2000)

    # Get mean exon incorporation per size
    inc_mean = inc_df.groupby("exon_size").mean()
    inc_mean = inc_mean.reset_index()
    gc_inc_merge = pd.merge(gc_mean, inc_mean, left_on="seq_len", right_on="exon_size", how="inner")

    scatterplot_roll_avg(gc_inc_merge, "seq_len", "gc_perc", 20,
                         "{}scatter_roll_avg_20_by_size_colour_exon_inc.png".format(outdir),
                         Q5, Q95, inc_border_low, inc_border_high, x_max=2000, hue="ratio")

    return None


def parse_bed(filename):
    """Parse bedtools nuc .bed file into a dataframe

    :param filename: str, filename of bed file
    :return: pandas dataframe object, parsed bed file with renamed columns
    """
    columns = {
        "#1_usercol": "chrom",
        "2_usercol": "start",
        "3_usercol": "stop",
        "4_usercol": "id",
        "5_pct_at": "at_perc",
        "6_pct_gc": "gc_perc",
        "7_num_A": "num_a",
        "8_num_C": "num_c",
        "9_num_G": "num_g",
        "10_num_T": "num_t",
        "11_num_N": "num_N",
        "12_num_oth": "num_oth",
        "13_seq_len": "seq_len",
    }
    df = pd.read_csv(filename, sep="\t", header=0)
    return df.rename(columns=columns)


def main():
    """Main function."""
    pickle, bed, quartile_file, inc_border_file, out_dir = parse_arguments()

    inc_df = pd.read_pickle(pickle)
    gc_df = parse_bed(bed)

    # Read in quartiles
    with open(quartile_file) as file:
        Q5, Q95 = file.readline().strip().split("\t")
    Q5, Q95 = int(Q5), int(Q95)

    # Read in positive exon incorporation borders
    with open(inc_border_file) as file:
        inc_border_low, inc_border_high = file.readline().strip().split("\t")
    inc_border_low, inc_border_high = int(inc_border_low), int(inc_border_high)

    analyse_gc(gc_df, inc_df, Q5, Q95, inc_border_low, inc_border_high, out_dir)

    # comparative_analysis(filenames)
    return None


if __name__ == "__main__":
    main()
