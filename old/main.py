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

    exon_df = parse_or_load_dataframe(filename, "./data/exons_df.pickle", start_time=start_time)

    # Calculate how many exons, transcripts and genes are present
    output_type_counts(exon_df, filename="{}general_statistics.txt".format(outdir))

    # Split dataframe in four different dataframes based on prime status
    prime_df_dict = prime_dataframes(exon_df)

    if start_time:
        print("Time to determine primes and add to dict: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    # Create length histograms
    bins = [i for i in range(0, 36000, 50)]
    exon_len_histograms(prime_df_dict, "./images/exon_length_histogram", bins)
    exon_len_histograms(prime_df_dict, "./images/exon_length_histogram_y_limited", bins, max_height=3000)

    if start_time:
        print("Time to create length histograms: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    # Bucket exons and create bar plots
    exon_buckets = [0, 100, 200, 300, 400, 500, 1000, 2000, 5000, 50000]
    bucket_dict = bucket_exons_and_plot(prime_df_dict, exon_buckets)
    name_dict = {
        (False, False): "Internal_exons",
        (False, True): "3'_exons",
        (True, False): "5'_exons",
        (True, True): "Single_exon_genes",
    }
    # Output information
    with open("{}bucket_information.txt".format(outdir), "w+") as f_out:
        for idx, bucket_series in bucket_dict.items():
            f_out.write("Bucket information for {}:\n{}\n".format(name_dict[idx], bucket_series))
            f_out.write("Total exons in all buckets: {}\n\n".format(sum(bucket_series)))

    if start_time:
        print("Time to bucket exons and plot: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    # Variable determining how many exons are taken along
    exon_quantity = 100

    for idx, prime_dict in prime_df_dict.items():
        # Retrieve associated genes and plot sizes of the largest exons
        gene_list = largest_exon_genes(prime_dict, exon_quantity, "./images/top_{}_{}_barplot_size".format(exon_quantity,
                                                                                                        name_dict[idx]))
        # Output gene list to file
        with open("./out/latest_{}_gene_list_usedforheatmap.txt".format(name_dict[idx]), "w+") as file:
            file.write(", ".join(gene_list))

        # Create heatmaps of the genes containing the top N largest exons
        duplicates, excluded = heatmap_expression_genes(gene_list, "./data/GTEx_gene_median_tpm.gct.gz",
                                                        "./images/heatmaps/top_{}_{}_heatmap_gene_expression"
                                                        "".format(name_dict[idx], exon_quantity))
        # Do the same for a random subset of the same size
        seeds = [i for i in range(1, 11)]
        for seed in seeds:
            random_gene_list = prime_dict.sample(n=exon_quantity, random_state=seed)
            heatmap_expression_genes(random_gene_list["geneID"].tolist(), "./data/GTEx_gene_median_tpm.gct.gz",
                                     "./images/heatmaps/random_{}_{}_heatmap_gene_expression".format(seed, exon_quantity))

        # Output information
        with open("{}expression_information_{}.txt".format(outdir, name_dict[idx]), "a+") as f_out:
            f_out.write("Out of the top {} genes, the following genes are present multiple times:\n"
                        "".format(exon_quantity))
            duplicates.sort(reverse=True, key=lambda x: x[1])
            for item, count in duplicates:
                f_out.write("{}\t{}\n".format(item, count))
            f_out.write("\n\nThe following genes were not found in the GTEx expression data:\n")
            for item in excluded:
                f_out.write("-{}\n".format(item))
            f_out.write("\n")

        if start_time:
            print("Time to create expression heatmaps {}: {:.2f} seconds\n".format(name_dict[idx],
                                                                                   time.perf_counter() - start_time))

    # Get genes per bucket:
    gene_list = genes_between_size(exon_df, 5000)
    print(len(gene_list))
    with open("./out/latest_gene_list_usedforheatmap.txt", "w+") as file:
        file.write(", ".join(gene_list))

    # TODO: think of a way to statistically calculate if the genes are differentially expressed in certain cell types
    # TODO: alternate names of genes have to be included for the expression data
    #   TODO: gene_synonym should be added somehow. Perhaps in a separate dict.
    #         This is located in gff[8] under gene_synonym=
    # TODO: split the expression on 5' and 3'? Or nah?
    # TODO: GO enrichment analysis on the top X genes?


if __name__ == "__main__":
    # If you change the input file, remove the pickle file!
    main("./data/GRCh38_latest_genomic.gff", clocking=True)
