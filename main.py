#!/usr/env/bin python3
"""
Author: Matthijs Pon
Description: Script to parse an GFF file and investigate different characteristics
             of exons.
"""
import os
import time
from calculate import *
from plots import *
from parse import *


def main(filename, outdir="./out/", clocking=False):
    start_time = None

    if clocking:
        start_time = time.perf_counter()

    # Either parse the file or load existing dataframe to save time
    if not os.path.exists("./data/exons_df.pickle"):
        # Parse file into df
        with open(filename) as file:
            exon_df = parse_gff_file(file, start_time)
        # Save df
        exon_df.to_pickle("./data/exons_df.pickle")
    else:
        # Load parsed dataframe
        exon_df = pd.read_pickle("./data/exons_df.pickle")

    # Calculate how many exons, transcripts and genes are present
    exons, transcripts, genes = exon_trans_gene_amount(exon_df)
    # Output information
    with open("{}general_statistics.txt".format(outdir), "w+") as f_out:
        f_out.write("No. of genes: \t\t{}\nNo. of transcripts: {}\n"
                    "No. of exons: \t\t{}\n".format(genes, transcripts, exons))
        f_out.write("Avg. no. of transcripts per gene: {:.2f}\nAvg. no. of exons per gene: {:.2f}\n"
                    "".format((transcripts / genes), (exons / genes)))

    # Add primes to dict
    prime_df_dict = prime_dataframes(exon_df)
    if start_time:
        print("Time to determine primes and add to dict: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    # Create length histograms
    bins = [i for i in range(0, 36000, 500)]
    exon_len_histograms(prime_df_dict, "./images/exon_length_histogram", bins)
    exon_len_histograms(prime_df_dict, "./images/exon_length_histogram_xlog_transformed", bins, ylog=True)

    if start_time:
        print("Time to create length histograms: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    # Bucket exons and create bar plots
    exon_buckets = [0, 100, 200, 300, 400, 500, 1000, 2000, 5000, 50000]
    bucket_dict = bucket_exons_and_plot(prime_df_dict, exon_buckets)
    name_dict = {
        (False, False): "internal exons",
        (False, True): "3' exons",
        (True, False): "5' exons",
        (True, True): "double-prime exons",
    }
    # Output information
    with open("{}bucket_information.txt".format(outdir), "w+") as f_out:
        for idx, bucket_series in bucket_dict.items():
            f_out.write("Bucket information for {}:\n{}\n".format(name_dict[idx], bucket_series))
            f_out.write("Total exons in all buckets: {}\n\n".format(sum(bucket_series)))

    if start_time:
        print("Time to bucket exons and plot: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    # Variable determining how many exons are taken along
    exon_quantity = 500
    # Retrieve associated genes and plot sizes of the largest exons
    gene_list = largest_exon_genes(exon_df, exon_quantity, "./images/top_{}_barplot_size".format(exon_quantity))
    # Create heatmaps of the genes containing the top N largest exons

    duplicates, excluded = heatmap_expression_genes(gene_list, "./data/GTEx_gene_median_tpm.gct.gz",
                                                    "./images/heatmaps/top_{}_heatmap_gene_expression".format(exon_quantity))
    # Do the same for a random subset of the same size
    seeds = [i for i in range(1, 11)]
    for seed in seeds:
        random_gene_list = exon_df.sample(n=exon_quantity, random_state=seed)
        heatmap_expression_genes(random_gene_list["geneID"].tolist(), "./data/GTEx_gene_median_tpm.gct.gz",
                                 "./images/heatmaps/random_{}_{}_heatmap_gene_expression".format(seed, exon_quantity))

    # Output information
    with open("{}expression_information.txt".format(outdir), "a+") as f_out:
        f_out.write("Out of the top {} genes, the following genes are present multiple times:\n".format(exon_quantity))
        duplicates.sort(reverse=True, key=lambda x: x[1])
        for item, count in duplicates:
            f_out.write("{}\t{}\n".format(item, count))
        f_out.write("\n\nThe following genes were not found in the GTEx expression data:\n")
        for item in excluded:
            f_out.write("-{}\n".format(item))
        f_out.write("\n")

    if start_time:
        print("Time to create expression heatmaps: {:.2f} seconds\n".format(time.perf_counter() - start_time))


if __name__ == "__main__":
    # If you change the input file, remove the pickle file!
    main("./data/GRCh38_latest_genomic.gff", clocking=True)
