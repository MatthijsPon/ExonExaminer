#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

from transcript_expression import determine_target_exon, yield_single_gene
from matplotlib import pyplot as plt
import statistics
import pandas as pd
import scipy.stats as stats
import os


def get_internal_exons(gene_df, remove=None):
    """Create a list of all non-prime exons in the gene

    :param gene_df: pandas dataframe, gene dataframe
    :param remove: str, exon ID (default: None)
    :return: set of str, exon IDs of non-prime exons
    """
    non_prime = set()
    only_exons = gene_df[gene_df["type"] == "exon"]
    # Select exons per transcript and apply True only if they are prime
    for idx, group in only_exons.groupby("parent"):
        exons_transcript = list(group["id"])
        exons_transcript.pop(0)  # Remove the first exon
        if len(exons_transcript) >= 1:  # Cannot pop twice if there is only one exon
            exons_transcript.pop(-1)  # Remove the last exon
        # All exons that have a non-prime function at least once are possible targets
        for exon in exons_transcript:
            non_prime.add(exon)
    if len(non_prime) >= 1:  # If there are internal exons return the largest, otherwise return none.
        if remove:
            non_prime.remove(remove)
            if len(non_prime) == 0:
                return None
        return non_prime
    return None


def exon_transcript_presence(gene_df, target_exon):
    """Determine in which transcripts an exon is present

    :param gene_df: pandas dataframe, gene dataframe
    :param target_exon: str, exon ID
    :return: numpy list, IDs of transcripts containing target exon
    """
    return gene_df.loc[gene_df["id"] == target_exon, "parent"].nunique()


def ratio_transcripts_exon(gene_df, target_exon=None):
    """Calculate the ratio of transcripts containing exon vs all transcripts

    :param gene_df: pandas dataframe, gene dataframe
    :param target_exon: str, exon ID (default: None)
    :return: float, ratio of transcripts containing exon
    """
    if not target_exon:
        target_exon = determine_target_exon(gene_df, target="largest_internal")
        if not target_exon:
            return None
    # Determine how many transcripts are present
    all_trans = gene_df.loc[gene_df["type"] == "transcript", "id"].nunique()
    # Determine how many transcripts contain exon
    target_trans = gene_df.loc[gene_df["id"] == target_exon, "parent"].nunique()
    target_ratio = float(target_trans) / all_trans
    # Get all non-target internal exons
    internal_exons = get_internal_exons(gene_df)  # , remove=target_exon)
    if not internal_exons:
        return None
    internal_ratios = []
    for internal_exon in internal_exons:
        int_ex_pres = exon_transcript_presence(gene_df, internal_exon)
        internal_ratios.append(float(int_ex_pres) / all_trans)
    internal_ratio = statistics.fmean(internal_ratios)
    return target_ratio - internal_ratio


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
    # Sample size too big for accurate normality testing.
    # stat, p = stats.shapiro(info["ratio"])
    # print("Shapiro Wilkins test statistic: {:.3f}, p-val {:.3f}".format(stat, p))


def main():
    """Main function."""

    # This is a temporary location for command line arguments, replace this later with a command line parser
    GFF_FILE = "./data/Homo_sapiens.GRCh38.106.chr.gff3"
    OUT_DIR = "./data/out/transcript_ratio"
    TEMP_DIR = "./data/temp"
    SIZES_OF_INTEREST = [0, 100, 200, 250, 400, 500, 1000, 2500, (1000, 10000)]

    if not os.path.exists("{}/longest_exon_transcript_ratio.pickle".format(TEMP_DIR)):
        exon_sizes = []
        ratios = []
        with open(GFF_FILE) as file:
            count = 0
            for gene in yield_single_gene(file):
                count += 1
                if count % 1000 == 0:
                    print("Gene {}".format(count))
                internal_exons = get_internal_exons(gene)
                if internal_exons is None:
                    continue  # If there are no internal exons, skip the ratios
                for exon in internal_exons:
                    ratio = ratio_transcripts_exon(gene, target_exon=exon)
                    if ratio:
                        ratios.append(ratio)
                        size = gene.loc[gene["id"] == exon, "size"].iloc[0]
                        exon_sizes.append(int(size))
        # Create DF from exon size and ratio
        info = pd.DataFrame({"exon_size": exon_sizes, "ratio": ratios})
        info.to_pickle("{}/longest_exon_transcript_ratio.pickle".format(TEMP_DIR))  # Save DF
    else:
        info = pd.read_pickle("{}/longest_exon_transcript_ratio.pickle".format(TEMP_DIR))  # Load DF if exists

    investigate_normality_plots(info, "ratio", OUT_DIR, "ratio")

    # Output Welch's T-test and dataframe distributions to file
    with open("{}/statistical_information.txt".format(OUT_DIR), "w+") as out_file:
        for size in SIZES_OF_INTEREST:
            total = None
            if type(size) == int:
                out_file.write("----- Exons of interest: {}+ bp -----\n".format(size))
                subset = info.loc[info["exon_size"] >= size]
                if size > 0:
                    total = info[~info["exon_size"].isin(subset["exon_size"])]
            else:
                out_file.write("----- Exons of interest: between {} and {} bp -----\n"
                               "No Welch T-test will be performed on this region.\n".format(size[0], size[1]))
                subset = info.loc[(info["exon_size"] >= size[0]) & (info["exon_size"] <= size[1])]

            out_file.write("Total average expression difference {}: {:.2f}\n"
                           "".format(size, statistics.fmean(list(subset["ratio"]))))
            out_file.write("No. of genes containing internal exons: {}\n".format(len(subset.index)))

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
            plt.savefig("{}/internal_exon_ratio_{}.png".format(OUT_DIR, size))
            plt.close()
    print("Finished")


if __name__ == '__main__':
    main()
