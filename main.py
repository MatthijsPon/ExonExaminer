#!/usr/env/bin python3
"""
Author: Matthijs Pon
Description: Script to parse an GFF file and investigate different characteristics
             of exons.
"""
import os
from calculate import *
from plots import *
from parse import *


def main(filename):
    # Either parse the file or load existing dataframe to save time
    if not os.path.exists("./data/exons_df.pickle"):
        # Parse file into df
        with open(filename) as file:
            exon_df = parse_gff_file(file)
        # Save df
        exon_df.to_pickle("./data/exons_df.pickle")
    else:
        # Load parsed dataframe
        exon_df = pd.read_pickle("./data/exons_df.pickle")

    # Calculate how many exons, transcripts and genes are present
    exons, transcripts, genes = exon_trans_gene_amount(exon_df)
    # Output information
    print(exons, transcripts, genes)

    # Add primes to dict
    prime_df_dict = prime_dataframes(exon_df)

    # Bucket exons and create histograms
    exon_buckets = [0, 100, 200, 300, 400, 500, 1000, 2000, 5000, 50000]
    bucket_exons_and_plot(prime_df_dict, exon_buckets)


if __name__ == "__main__":
    # If you change the input file, remove the pickle file!
    main("./data/GRCh38_latest_genomic.gff")
