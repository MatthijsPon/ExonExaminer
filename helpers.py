#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
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


def determine_gff(gff_8, identifier, further_split=None):
    """Select a specific object from the information column (column 8 in 0-count)

    :param gff_8: list of str, information column, split
    :param identifier: str, identifier of column of interest (e.g. ID=). Everything
                       identifier will be removed from return string
    :param further_split: str, seperator on which the column should be split on
                          after removal of the identifier (default: None)
    :return: str, information following the identifier, otherwise None
    """
    for item in gff_8.split(";"):
        if item.startswith(identifier):
            if further_split:
                return item[len(identifier):].split(further_split)
            return item[len(identifier):]
    return None


def yield_single_gene(file_object):
    """Yield a single gene from a gff3 file, including related protein coding transcripts and exons.

    :param file_object: open file object, gff3 file
    :return: pandas dataframe, gene information containing type, id, parent, size and alias columns
    """
    gene_dict = {
        "type": [],
        "id": [],
        "chromosome": [],
        "start": [],
        "stop": [],
        "strand": [],
        "parent": [],
        "size": [],
    }

    # A dict needs all three to be complete
    gene_id = ""
    transcript = []
    exon = []

    for line in yield_gff_line(file_object):
        if line[2] == "gene":
            if gene_id and transcript and exon:
                yield pd.DataFrame.from_dict(gene_dict)
            # Reset stats
            gene_dict = {
                "type": [],
                "id": [],
                "chromosome": [],
                "start": [],
                "stop": [],
                "strand": [],
                "parent": [],
                "size": [],
            }
            gene_id = determine_gff(line[8], "ID=gene:")
            transcript = []
            exon = []
            # Add data to dict
            gene_dict["type"].append("gene")
            gene_dict["id"].append(gene_id)
            gene_dict["chromosome"].append(line[0])
            gene_dict["start"].append(line[3])
            gene_dict["stop"].append(line[4])
            gene_dict["strand"].append(line[6])
            gene_dict["parent"].append(None)
            gene_dict["size"].append(None)
        elif line[2] == "mRNA":
            # only add if there is a parent gene and the mRNA is protein_coding
            if determine_gff(line[8], "Parent=gene:") == gene_id and \
                    determine_gff(line[8], "biotype=") == "protein_coding":
                trans_id = determine_gff(line[8], "ID=transcript:")
                transcript.append(trans_id)
                # Add data to dict
                gene_dict["type"].append("transcript")
                gene_dict["id"].append(trans_id)
                gene_dict["chromosome"].append(line[0])
                gene_dict["start"].append(line[3])
                gene_dict["stop"].append(line[4])
                gene_dict["strand"].append(line[6])
                gene_dict["parent"].append(gene_id)
                gene_dict["size"].append(None)
        elif line[2] == "exon":
            if determine_gff(line[8], "Parent=transcript:") in transcript:
                exon_id = determine_gff(line[8], "Name=")
                exon.append(exon_id)
                gene_dict["type"].append("exon")
                gene_dict["id"].append(exon_id)
                gene_dict["chromosome"].append(line[0])
                gene_dict["start"].append(line[3])
                gene_dict["stop"].append(line[4])
                gene_dict["strand"].append(line[6])
                gene_dict["parent"].append(determine_gff(line[8], "Parent=transcript:"))
                gene_dict["size"].append((int(line[4]) - int(line[3]) + 1))  # The size is incl. the last nucleotide

    if gene_id and transcript and exon:
        yield pd.DataFrame.from_dict(gene_dict)


def get_internal_exons(gene_df, remove=None):
    """Create a list of all internal exons in the gene dataframe

    :param gene_df: pandas dataframe, gene dataframe
    :param remove: str, exon ID (default: None)
    :return: set of str, exon IDs of internal exons (otherwise None)

    Select exons per transcript and determine which are internal. An internal exon has a non-prime position in at
    least one transcript. Returns all unique exon IDs.
    """
    non_prime = set()
    only_exons = gene_df[gene_df["type"] == "exon"]

    for idx, group in only_exons.groupby("parent"):
        exons_transcript = list(group["id"])
        exons_transcript.pop(0)  # Remove the first exon
        if len(exons_transcript) >= 1:  # Cannot pop twice if there is only one exon
            exons_transcript.pop(-1)  # Remove the last exon
        for exon in exons_transcript:
            non_prime.add(exon)
    # Only return a value when there are internal exons
    if len(non_prime) >= 1:
        if remove:
            non_prime.remove(remove)
            if len(non_prime) == 0:
                return None
        return non_prime
    return None


def determine_target_exon(gene_df, target):
    """Determine which exon to target

    :param gene_df: pandas dataframe, gene dataframe
    :param target: str, "largest" or "largest_internal" to decide target
    :return: str, geneID if found, otherwise None
    """
    if target == "largest":
        return gene_df[gene_df["size"] == gene_df["size"].max()]["id"].iloc[0]  # Select largest exon
    elif target == "largest_internal":
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
            subset_internal = gene_df[gene_df["id"].isin(non_prime)]  # Select only internal exons from gene_df
            # Select largest exon from the subset
            largest = subset_internal[subset_internal["size"] == subset_internal["size"].max()]["id"].iloc[0]
            return largest
        return None
    else:
        raise ValueError("Incorrect target given in determine_target_exon().")