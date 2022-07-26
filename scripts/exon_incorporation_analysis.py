#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: This script parses an ENSEMBL gff3 file and determines the average
             inclusion of exons compared to their length.
"""
import math

from helpers import get_internal_exons, yield_single_gene, determine_target_exon
from matplotlib import pyplot as plt

import os
import argparse as arg
import pandas as pd
import scipy.stats as stats
import statistics
import seaborn as sns


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine the average inclusion of exons compared to their length.")
    parser.add_argument("pickle_file", type=str,
                        help="pandas pickle file created by exon_incorporation_pickle.py")
    parser.add_argument("quartile_file", type=str,
                        help="text file containing the Q5 and Q95, created by gff_statistics.py")
    parser.add_argument("out_dir", type=str,
                        help="output directory (incl. trailing \"/\")")
    args = parser.parse_args()

    if not os.path.exists(args.pickle_file):
        raise FileNotFoundError("Pickle file does not exist.")

    return args.pickle_file, args.quartile_file, args.out_dir


def cumulative_barplot(data, x, y, filename):
    """Create a cumulative bar plot.

    :param data: pandas dataframe object, with data to plot
    :param x: str/index, index to access x values of data
    :param y: str/index, index to access y values of data
    :param filename: str, output filename for figure
    :return: None, files are created
    """
    # Get highest and lowest value of cumulative bar plot
    low = int(data.loc[data[y] == data[y].min(), "size"])
    high = int(data.loc[data[y] == data[y].max(), "size"])

    plt.rcParams.update({"font.size": 16, "figure.figsize": (19.2, 10.8)})
    # Update axes
    plt.xticks([0, 50, 150, 250, 350, 500, 750, 1000, 1250, 1500, 1750, 2000])
    plt.xlabel("Size (nt)")
    plt.ylabel("Cumulative ratio (all genes)")
    plt.yticks([0])
    plt.axvline(x=low, color="y", linestyle=":", label="Positive incorporation start ({} nt)".format(low))
    plt.axvline(x=high, color="purple", linestyle=":", label="Positive incorporation end ({} nt)".format(high))
    plt.bar(x=data[x], height=data[y])
    plt.legend()
    plt.savefig(filename)
    plt.close()
    return low, high


def scatterplot_roll_avg(data, x, y, hue, window, out_file,
                         Q5=49, Q95=288, inc_bottom=69, inc_top=199,
                         x_min=None, x_max=None):
    """Create a scatter plot with a moving average.

    :param data: 2d object (e.g. pandas dataframe) with data to plot
    :param x: str/index, index to access x values of data
    :param y: str/index, index to access y values of data
    :param hue: str/index, index to access hue values of data
    :param window: int, no. of data points to use in the rolling average
    :param out_file: str, output filename
    :param Q5: int, 5% quartile size
    :param Q95: int, 95% quartile size
    :param inc_bottom: int, size corresponding to the bottom of the cumulative bar plot
    :param inc_top: int, size corresponding to the top of the cumulative bar plot
    :param x_min: int, limit of the left side of the x scale
    :param x_max: int, limit of the right side of the x scale
    :return: None, files are created
    """
    plt.rcParams.update({"font.size": 20, "figure.figsize": (20, 7), "legend.loc": "upper left"})
    plt.grid()
    if x_max:
        plt.xlim(right=x_max)
    elif x_min:
        plt.xlim(left=x_min)
    elif x_min and x_max:
        plt.xlim(x_min, x_max)

    # Linegraph with continuous colour palette
    cmap = sns.color_palette("crest", as_cmap=True)
    points = plt.scatter(x=data[x], y=data[y], c=data[hue], cmap=cmap)
    cbar = plt.colorbar(points)
    cbar.set_label("No. of transcripts in gene (max. 10)", rotation=270)

    # Quartile lines
    plt.axvline(x=Q5, color="r", linestyle=":", label="Q05 ({} nt)".format(Q5))
    plt.axvline(x=Q95, color="g", linestyle=":", label="Q95 ({} nt)".format(Q95))
    plt.axvline(x=inc_bottom, color="y", linestyle=":", label="Positive incorporation start ({} nt)".format(inc_bottom))
    plt.axvline(x=inc_top, color="purple", linestyle=":", label="Positive incorporation end ({} nt)".format(inc_top))

    # rolling average line
    data['MA_group'] = data[y].rolling(window=window).mean()
    plt.plot(data[x], data["MA_group"], color="r", label="Rolling average of {} data points".format(window))

    plt.xlabel("Exon size (nt)")
    plt.ylabel("Average exon incorporation")
    plt.legend()
    plt.savefig(out_file)
    plt.close()
    return None


def main():
    """Main function."""
    pickle_file, quartile_file, out_dir = parse_arguments()

    results = pd.read_pickle(pickle_file)

    # Output genes which have a higher and lower incorporation ratio than 0
    # for further analysis in the expression script
    ratio_higher = results.loc[results["ratio"] >= 0, "id"]
    ratio_lower = results.loc[results["ratio"] < 0, "id"]
    with open("{}/genes_ratio_biggerequal_0.csv".format(out_dir), "w+") as outfile:
        outfile.write(",".join(ratio_higher))
    with open("{}/genes_ratio_smaller_0.csv".format(out_dir), "w+") as outfile:
        outfile.write(",".join(ratio_lower))


    data = {
        "size": [],
        "sum_smaller": [],
        "mean_group": [],
        "mean_trans": [],
        "count": []
    }
    step = 1
    sizes_of_interest = [i for i in range(0, 2000, step)]
    for size in sizes_of_interest:
        subset = results.loc[results["exon_size"] >= size]
        subset_small = results.loc[results["exon_size"] <= size]

        if len(list(subset.loc[subset["exon_size"] < (size + step), "ratio"])) > 0:
            data["size"].append(size)
            data["sum_smaller"].append(sum(list(subset_small["ratio"])))
            data["mean_group"].append(statistics.fmean(list(subset.loc[subset["exon_size"] < (size + step), "ratio"])))
            # Take the number of transcripts, with 10 as the max
            # (otherwise the colour bar gets saturated by extreme values)
            data["mean_trans"].append(
                min(statistics.fmean(list(subset.loc[subset["exon_size"] < (size + step), "no_trans"])), 10))
            data["count"].append(len(subset.index))
    out = pd.DataFrame(data)

    # Cumulative bar plot
    inc_border_low, inc_border_high = \
        cumulative_barplot(out, "size", "sum_smaller",
                           "{}/cumulative_barplot.png".format(out_dir))

    # Get quartile values and convert to integer
    with open(quartile_file) as file:
        Q5, Q95 = file.readline().strip().split("\t")

    Q5, Q95 = int(Q5), int(Q95)

    # Scatterplots with different rolling average windows
    scatterplot_roll_avg(out, "size", "mean_group", "mean_trans", 15,
                         "{}/scatterplot_rolling_avg_window15.png".format(out_dir),
                         Q5, Q95, inc_border_low, inc_border_high)
    scatterplot_roll_avg(out, "size", "mean_group", "mean_trans", 20,
                         "{}/scatterplot_rolling_avg_window20.png".format(out_dir),
                         Q5, Q95, inc_border_low, inc_border_high)
    scatterplot_roll_avg(out, "size", "mean_group", "mean_trans", 25,
                         "{}/scatterplot_rolling_avg_window25.png".format(out_dir),
                         Q5, Q95, inc_border_low, inc_border_high)

    # Output lowest and highest value, to determine peaks
    with open("{}/exon_incorporation_positive_borders.txt".format(out_dir), "w+") as file:
        file.write("{}\t{}".format(inc_border_low, inc_border_high))

    return None


if __name__ == '__main__':
    main()


