#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Create figures and plots
"""

import matplotlib.pyplot as plt
import seaborn as sns


def bucket_barplot(plot_dict, fig_name):
    """Make four bar plots of the bucketed exons

    :param plot_dict: a dictionary containing the bucketed dataframes
    :param fig_name: str, name of output file (excl. ".png"), incl leading folder location
    :return: None, function creates .png file
    """

    name_dict = {
        (False, False): ((0, 0), "Internal exons"),
        (False, True): ((0, 1), "3' exons"),
        (True, False): ((1, 0), "5' exons"),
        (True, True): ((1, 1), "both-prime exons"),
    }
    # Make figure containing 2x2 subplots
    fig, axes = plt.subplots(2, 2, sharey="all", figsize=(19.2, 10.8))

    # Make each subplot
    for idx, value in plot_dict.items():
        # Create the plot and set name and axis titles
        value.plot.bar(ax=axes[name_dict[idx][0]], title=name_dict[idx][1],
                       xlabel="bin borders (bp)", ylabel="relative amount of exons in bin")
        # Change the x-axis values
        axes[name_dict[idx][0]].set_xticklabels(axes[name_dict[idx][0]].xaxis.get_majorticklabels(),
                                                rotation=0, fontsize=8)

    fig.savefig("{}.png".format(fig_name))
    # Close figure, otherwise object keeps existing and bins are added cumulatively
    plt.close()


def exon_len_histograms(exon_prime_df, fig_name, bins, ylog=False):
    """Make four histograms of the exon lengths

    :param exon_prime_df: exon dataframes, grouped on primes
    :param fig_name: str, name of output file (excl. ".png"), incl leading folder location
    :param bins: list of ints, bin borders
    :param ylog: bool, log transforms y-axis if True
    :return: None, function creates .png file
    """

    name_dict = {
        (False, False): ((0, 0), "Internal exons"),
        (False, True): ((0, 1), "3' exons"),
        (True, False): ((1, 0), "5' exons"),
        (True, True): ((1, 1), "both-prime exons"),
    }
    # Make figure containing 2x2 subplots
    if ylog:
        fig, axes = plt.subplots(2, 2, sharex="all", sharey="all", figsize=(19.2, 10.8))
    else:
        fig, axes = plt.subplots(2, 2, sharex="all", figsize=(19.2, 10.8))

    # Make each subplot
    for idx, value in exon_prime_df.items():
        # Create the plot and set name and axis titles
        value["size"].plot.hist(ax=axes[name_dict[idx][0]], bins=bins, title=name_dict[idx][1])
        axes[name_dict[idx][0]].set_xlabel("exon size (500bp buckets)")

        # log transform y axes
        if ylog:
            axes[name_dict[idx][0]].set_yscale('log')
            axes[name_dict[idx][0]].set_ylabel('log(Frequency)')

    fig.savefig("{}.png".format(fig_name))
    # Close figure, otherwise object keeps existing and bins are added cumulatively
    plt.close()


def expression_heatmap(expression_df, filename):
    """"""
    # Drop the name column, which is ENSEMBL gene codes
    expression_df = expression_df.drop(columns="Name")
    # Set the actual gene names as index for the figure
    expression_df = expression_df.set_index("Description")
    # Set figure size before creating figure, otherwise it won't work
    # plt.figure(figsize=(19.2, 10.8))
    plt.rcParams.update({"font.size": 8, "figure.figsize": (30, 20)})
    sns.heatmap(expression_df)
    plt.savefig("{}.png".format(filename))
    plt.close()

    # Normalize heatmap
    df_norm_row = expression_df.apply(lambda x: x / x.max(), axis=1)

    # Check what values are causing NA during normalization:
    # with pd.option_context("display.max_rows", None, "display.max_columns", None):
    #     print(expression_df[df_norm_row.isna().any(axis=1)])

    # Values with 0 in the whole row are converted to NaN (as 0 / 0 occurs), so replace those with 0
    df_norm_row.fillna(0, inplace=True)
    plt.rcParams.update({"font.size": 8, "figure.figsize": (30, 20)})
    sns.heatmap(df_norm_row)
    plt.savefig("{}_normalized.png".format(filename))
    plt.close()
