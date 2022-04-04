#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Create figures and plots
"""

import matplotlib.pyplot as plt


def bucket_histogram(plot_dict, fig_name):
    """Make a histogram of the bucketed exons

    :param plot_dict: a dictionary containing the bucketed dataframes
    :param fig_name: str, name of output file (excl. ".png"), incl leading folder location
    :return: None, function creates .png file
    """

    name_dict = {
        (False, False): ((0, 0), "Non-prime exons"),
        (False, True): ((0, 1), "Three-prime exons"),
        (True, False): ((1, 0), "Five-prime exons"),
        (True, True): ((1, 1), "Five- and three-prime exons"),
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
