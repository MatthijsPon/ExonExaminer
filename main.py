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


def main(filename, clocking=False):
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
    # TODO output stats to file instead of stdout
    print("No. of genes: \t\t{}\nNo. of transcripts: {}\n"
          "No. of exons: \t\t{}\n".format(genes, transcripts, exons))
    print("Avg. no. of transcripts per gene: {:.2f}\nAvg. no. of exons per transcript: {:.2f}\n"
          "".format((transcripts / genes), (exons / transcripts)))

    # Add primes to dict
    prime_df_dict = prime_dataframes(exon_df)
    if start_time:
        print("Time to determine primes and add to dict: {:.2f} seconds".format(time.perf_counter() - start_time))

    # Bucket exons and create histograms
    exon_buckets = [0, 100, 200, 300, 400, 500, 1000, 2000, 5000, 50000]
    bucket_exons_and_plot(prime_df_dict, exon_buckets)
    if start_time:
        print("Time to bucket exons and plot: {:.2f} seconds".format(time.perf_counter() - start_time))


if __name__ == "__main__":
    # If you change the input file, remove the pickle file!
    main("./data/GRCh38_latest_genomic.gff", clocking=True)
