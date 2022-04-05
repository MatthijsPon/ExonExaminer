#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Create figures and plots
"""

import matplotlib.pyplot as plt


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
        value["size"].plot.hist(ax=axes[name_dict[idx][0]], bins=bins, title=name_dict[idx][1],
                                xlabel="size (bp)", ylabel="relative amount of exons in bin")
        # log transform y axes
        if ylog:
            axes[name_dict[idx][0]].set_yscale('log')

    fig.savefig("{}.png".format(fig_name))
    # Close figure, otherwise object keeps existing and bins are added cumulatively
    plt.close()
