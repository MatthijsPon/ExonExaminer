#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description: Parse the exons of a gff file into a pandas dataframe
"""
import time

import pandas as pd


def yield_gff_line(file_object):
    """Iterator, Yields a single non-comment line of a tsv gff file

    :param file_object: iterable, open file object or list of lines
    :return: yields a list, a gff line split on \t
    """
    for line in file_object:
        if line.startswith("#"):
            # Skip the comment line, which is not of interest
            continue
        line = line.strip()
        yield line.split("\t")


def present_gff_elements(file_object):
    """Parse all unique gff elements (gff[2]).
    
    :param file_object: open file object, gff file
    :return: set containing strings, gff elements
    """
    elements = set()
    for line in yield_gff_line(file_object):
        elements.add(line[2])
    return elements


def parse_exon_gff_to_dict(gff_list):
    """Parse a single exon gff line to a dict.

    :param gff_list: gff line split on '\t'
    :return: dict, column name: data
    """
    split_gff_8 = gff_list[8].split(";")
    gene_id = ""

    if split_gff_8[4].startswith("gene="):
        gene_id = split_gff_8[4].split("=")[1]
    else:
        # Sometimes, the order of the gff[8] is different, this is for these edge cases
        for item in split_gff_8:
            if item.startswith("gene="):
                gene_id = item.split("=")[1]

    exon_dict = {
        "chromosomeID": gff_list[0],
        "exonID": ["-".join(split_gff_8[0].split("-")[1:])],
        "transcriptID": "-".join(split_gff_8[1].split("-")[1:]),
        "geneID": gene_id,
        "start": int(gff_list[3]),
        "stop": int(gff_list[4]),
        "strand": gff_list[6],
        "size": int(gff_list[4]) - int(gff_list[3])
    }
    return exon_dict


def determine_prime_exons(exon_df):
    """Determine 5' and 3' exons by grouping exon based on transcripts and selecting the first and last exon.

    :param exon_df: pandas dataframe object
    :return: exon dataframe with the primes included
    """
    five_prime = []
    three_prime = []

    for idx, group in exon_df.groupby("transcriptID"):
        first = group.head(1)
        last = group.tail(1)
        # Collect the indices of 5' and 3' exons
        if first.iloc[0]["strand"] == "+":
            five_prime.append(first.index.values.astype(int)[0])
            three_prime.append(last.index.values.astype(int)[0])
        else:
            five_prime.append(last.index.values.astype(int)[0])
            three_prime.append(first.index.values.astype(int)[0])

    add_five_prime = [False for i in range(len(exon_df.index))]
    add_three_prime = [False for i in range(len(exon_df.index))]

    # There has to be a more effective way with list comprehension, but I do not know
    for i in five_prime:
        add_five_prime[i] = True
    for i in three_prime:
        add_three_prime[i] = True

    # Since the dataframes are 0-count as are all python objects, the lists can simply be added
    exon_df["five_prime"] = add_five_prime
    exon_df["three_prime"] = add_three_prime

    return exon_df


def merge_duplicate_exon_rows(group):
    """Merge duplicate exons based on start and stop location.

    :param group: dataframe containing only duplicate rows
    :return: one single dict, which can be turned into a dataframe
    """
    merged_dict = {
        "chromosomeID": "",
        # Sets are unique and will remove any accidental double exonID/transcriptIDs
        "exonID": set(),
        "transcriptID": set(),
        "geneID": "",
        "start": 0,
        "stop": 0,
        "strand": "",
        "size": 0,
        "five_prime": False,
        "three_prime": False
    }

    for index, row in group.iterrows():
        # Sanity checks for actual duplication events:
        if merged_dict["start"] != 0 and row["start"] != merged_dict["start"]:
            raise ValueError("For some reason, the start values are different")

        if merged_dict["stop"] != 0 and row["stop"] != merged_dict["stop"]:
            raise ValueError("For some reason, the stop values are different")

        if merged_dict["strand"] and merged_dict["strand"] != row["strand"]:
            print("start: {}\nstop: {}".format(row["start"], row["stop"]))
            raise ValueError("Exons are on different strands but share the same locations")

        merged_dict["chromosomeID"] = row["chromosomeID"]
        merged_dict["exonID"].add(row["exonID"][0])
        merged_dict["transcriptID"].add(row["transcriptID"])
        merged_dict["geneID"] = row["geneID"]
        merged_dict["start"] = row["start"]
        merged_dict["stop"] = row["stop"]
        merged_dict["size"] = row["size"]
        merged_dict["strand"] = row["strand"]

        # If either exon acts as a five or three prime, set value to true
        if row["five_prime"]:
            merged_dict["five_prime"] = True
        if row["three_prime"]:
            merged_dict["three_prime"] = True

    # Convert exonID and transcriptID back to list
    merged_dict["exonID"] = list(merged_dict["exonID"])
    merged_dict["transcriptID"] = list(merged_dict["transcriptID"])
    return merged_dict


def merge_duplicate_exon_df(exons_df):
    """Merge duplicate exon rows within a dataframe.

    :param exons_df: pandas dataframe object
    :return: exons dataframe with unique exons (based on location)
    """
    # Possible improvement: determine which rows are duplicates, drop them
    # based on index and add back the merged rows to save memory (maybe time?)

    unique_list = []
    for idx, group in exons_df.groupby(["chromosomeID", "start", "stop", "strand"]):
        if len(group) > 1:
            add = merge_duplicate_exon_rows(group)
            unique_list.append(add)
        else:
            add = group.to_dict(orient="records")[0]
            # Transform the transcriptID into a list which also happens in merge_duplicate_exon_rows()
            add["transcriptID"] = [add["transcriptID"]]
            unique_list.append(add)

    return pd.DataFrame(unique_list)


def parse_gff_file(file_object, start_time=None):
    """Iterate through gff lines and parse exons into a pandas dataframe containing relevant information.

    :param file_object: iterable, open file object or list of lines
    :param start_time: time object, only present if clocking is on in main function
    :return: exon_df: dataframe containing deduplicated exons
    """
    genes = {}
    transcripts = {}
    exons = []
    for gff in yield_gff_line(file_object):

        if gff[2] in ["gene", "pseudogene"]:
            # Add genes to the gene dict
            gene_id = gff[8].split(";")[0][3:]
            # If more information is needed, add it here!
            # Init empty list for the insertions of transcripts
            genes[gene_id] = []

        elif gff[2] in ["transcript", "mRNA", "primary_transcript"]:
            parent_gene = gff[8].split(";")[1][7:]
            # Only transcripts that belong to selected genes are interesting
            if genes.get(parent_gene) is not None:
                trans_id = gff[8].split(";")[0][3:]
                # Add transcript to gene information
                trans_list = genes.get(parent_gene)
                trans_list.append(trans_id)
                genes[parent_gene] = trans_list
                # Init empty list for the insertion of exons
                transcripts[trans_id] = []

        elif gff[2] == "exon":
            parent_trans = gff[8].split(";")[1][7:]
            # Only exons that belong to selected transcripts are interesting
            if transcripts.get(parent_trans) is not None:
                exons.append(parse_exon_gff_to_dict(gff))
    if start_time:
        print("Time to read in file: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    exon_df = pd.DataFrame(exons)
    if start_time:
        print("Time to convert to dataframe: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    exon_df = determine_prime_exons(exon_df)
    if start_time:
        print("Time to determine prime exons: {:.2f} seconds\n".format(time.perf_counter() - start_time))

    exon_df = merge_duplicate_exon_df(exon_df)
    if start_time:
        print("Time to merge duplicates: {:.2f} seconds\n".format(time.perf_counter() - start_time))
    return exon_df


# Testing on a small dataset
def testing():
    with open("./data/GRCh38_latest_genomic.gff") as file:
        exon = parse_gff_file(file, start_time=time.perf_counter())


if __name__ == '__main__':
    testing()
