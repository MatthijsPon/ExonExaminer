#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: This script parses an ENSEMBL gff3 file and determines the average
             inclusion of exons compared to their length.
"""
from helpers import get_internal_exons, yield_single_gene, determine_target_exon
from matplotlib import pyplot as plt

import os
import argparse as arg
import pandas as pd
import scipy.stats as stats
import statistics


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Determine the average inclusion of exons compared to their length.")
    parser.add_argument("pickle_file", type=str,
                        help="pandas pickle file create by exon_incorporation_pickle.py")
    parser.add_argument("out_dir", type=str,
                        help="output directory (incl. trailing \"/\"), default: ./data/out/")
    parser.add_argument("sizes_of_interest", type=int, nargs="+",
                        help="size cutoff(s) for selecting exon sizes")
    args = parser.parse_args()
    return args.pickle_file, args.out_dir, args.sizes_of_interest


def investigate_normality_plots(data, column_of_interest, out_dir, filename):
    """Create a histogram and QQplot vs. normal distribution for investigating normality on the data

    :param data: pandas dataframe, dataframe containing data
    :param column_of_interest: str, column of interest in the dataframe
    :param out_dir: str, output directory
    :param filename: str, start of filename for images to be saved to
    :return: None, images are created
    """
    # Histogram
    plt.rcParams.update({"font.size": 20, "figure.figsize": (19.2, 10.8)})
    data[column_of_interest].hist()  # Histogram
    plt.savefig("{}/{}_histogram.png".format(out_dir, filename))
    plt.close()
    # QQ plot
    plt.rcParams.update({"font.size": 20, "figure.figsize": (19.2, 10.8)})
    stats.probplot(data[column_of_interest], dist="norm", plot=plt)  # QQ plot
    plt.savefig("{}/{}_QQ_plot_vs_normal_distribution.png".format(out_dir, filename))
    plt.close()
    return None


def statistical_information(result_df, out_file, out_dir, size):
    """Calculate Welsch T-test if possible and output some statistical information

    :param result_df: pandas dataframe, ratios and sizes
    :param out_file: str, file to output information to
    :param out_dir: str, output directory path
    :param size: int, size value to split the data on
    :return: None, all information is outputted into outfile and figure
    """
    total = None
    out_file.write("----- Exons of interest: {}+ bp -----\n".format(size))
    subset = result_df.loc[result_df["exon_size"] >= size]
    if size > 0:
        total = result_df[~result_df["exon_size"].isin(subset["exon_size"])]

    out_file.write("Total average expression difference {}: {:.2f}\n"
                   "".format(size, statistics.fmean(list(subset["ratio"]))))
    out_file.write("No. of internal exons examined: {}\n".format(len(subset.index)))

    # Get mean, median, mode, quartiles, etc. for target group
    info_sub = subset.describe().transpose()
    out_file.write("\nInformation target group:\n")
    out_file.write(info_sub.to_string())
    out_file.write("\n")

    if total is not None:
        # Get mean, median, mode, quartiles, etc. for leftover group
        info_total = total.describe().transpose()
        out_file.write("\nInformation rest group:\n")
        out_file.write(info_total.to_string())
        out_file.write("\n")

        # Perform Welch T-test
        welch_t, p_val = stats.ttest_ind(total["ratio"], subset["ratio"], equal_var=False)
        out_file.write("Welch t-statistic: {}, p-value: {}\n".format(welch_t, p_val))

    out_file.write("\n\n")
    plt.rcParams.update({"font.size": 20, "figure.figsize": (19.2, 10.8)})
    subset.plot(x="exon_size", y="ratio", style="o")
    plt.savefig("{}/internal_exon_ratio_{}+.png".format(out_dir, size))
    plt.close()
    return None


def scatterplot_roll_avg(data, x, y, window, out_dir, filename):
    """Create a scatter plot with a moving average.

    :param data: 2d object (e.g. pandas dataframe) with data to plot
    :param x: str/index, index to access x values of data
    :param y: str/index, index to access y values of data
    :param window: int, no. of data points to use in the rolling average
    :param out_dir: str, output directory
    :param filename: str, output filename
    :return: None, files are created
    """
    plt.rcParams.update({"font.size": 20, "figure.figsize": (19.2, 10.8)})
    plt.grid()
    # Line graph per group
    plt.scatter(x=data[x], y=data[y])
    # Quartile lines
    plt.axvline(x=49, color="r", linestyle=":", label="Q05 (size: 49)")
    plt.axvline(x=122, color="y", linestyle=":", label="Q50 (size: 122)")
    plt.axvline(x=288, color="g", linestyle=":", label="Q95 (size: 288)")
    # rolling average line
    data['MA_group'] = data[y].rolling(window=window).mean()
    plt.plot(data[x], data["MA_group"], color="r", label="Rolling average of {} data points".format(window))
    plt.xlabel("Exon size (nt)")
    plt.ylabel("Average usage")
    plt.legend()
    plt.savefig("{}/scatterplot_ratio_{}.png".format(out_dir, filename))
    plt.close()
    return None


def cumulative_barplot(data, x, y, out_dir, filename):
    """Create a comolative barplot.

    :param data: pandas dataframe object, with data to plot
    :param x: str/index, index to access x values of data
    :param y: str/index, index to access y values of data
    :param out_dir: str, output directory
    :param filename: str, output filename
    :return: None, files are created
    """
    plt.rcParams.update({"font.size": 16, "figure.figsize": (19.2, 10.8)})
    # Update axes
    plt.xticks([0, 50, 150, 250, 350, 500, 750, 1000, 1250, 1500, 1750, 2000])
    plt.xlabel("Size (nt)")
    plt.ylabel("Cumulative ratio (all genes)")
    plt.yticks([0])
    plt.bar(x=data[x], height=data[y])
    # Human genome quartiles (.05, .50, .95)
    plt.axvline(x=49, color="r", linestyle=":", label="Q05 (size: 49)")
    plt.axvline(x=122, color="y", linestyle=":", label="Q50 (size: 122)")
    plt.axvline(x=288, color="g", linestyle=":", label="Q95 (size: 288)")
    plt.legend()
    plt.savefig("{}/{}.png".format(out_dir, filename))
    plt.close()

    # Output lowest and highest value, to determine peaks
    print(data.loc[data[y] == data[y].min()])
    print(data.loc[data[y] == data[y].max()])
    return None


def main():
    """Main function."""
    pickle_file, out_dir, sizes_of_interest = parse_arguments()

    if not os.path.exists(pickle_file):
        raise FileNotFoundError("Pickle file does not exist.")
    results = pd.read_pickle(pickle_file)

    ratio_higher = results.loc[results["ratio"] >= 0, "id"]
    ratio_lower = results.loc[results["ratio"] < 0, "id"]

    with open("{}/genes_ratio_biggerequal_0.csv".format(out_dir), "w+") as outfile:
        outfile.write(",".join(ratio_higher))
    with open("{}/genes_ratio_smaller_0.csv".format(out_dir), "w+") as outfile:
        outfile.write(",".join(ratio_lower))

    investigate_normality_plots(results, "ratio", out_dir, "ratio")

    # Output Welch's T-test and dataframe distributions to file
    with open("{}/statistical_information.txt".format(out_dir), "w+") as out_file:
        for size in sizes_of_interest:
            statistical_information(results, out_file, out_dir, size)

    data = {
        "size": [],
        "sum_smaller": [],
        "mean_group": [],
        "count": []
    }
    # TEMP
    # TODO remove single transcript genes?
    step = 1
    sizes_of_interest = [i for i in range(0, 2000, step)]
    for size in sizes_of_interest:
        subset = results.loc[results["exon_size"] >= size]
        subset_small = results.loc[results["exon_size"] <= size]

        if len(list(subset.loc[subset["exon_size"] < (size + step), "ratio"])) > 0:
            data["size"].append(size)
            data["sum_smaller"].append(sum(list(subset_small["ratio"])))
            data["mean_group"].append(statistics.fmean(list(subset.loc[subset["exon_size"] < (size + step), "ratio"])))
            data["count"].append(len(subset.index))
    out = pd.DataFrame(data)

    scatterplot_roll_avg(out, "size", "mean_group", 15, out_dir, "step1_window15")
    scatterplot_roll_avg(out, "size", "mean_group", 20, out_dir, "step1_window20")
    scatterplot_roll_avg(out, "size", "mean_group", 25, out_dir, "step1_window25")
    cumulative_barplot(out, "size", "sum_smaller", out_dir, "barplot_cumulative")


if __name__ == '__main__':
    main()
