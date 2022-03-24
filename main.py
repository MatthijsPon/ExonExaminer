#!/usr/env/bin python3
"""
Author: Matthijs Pon
Description: Script to parse an GFF file and investigate different characteristics
             of exons.
"""


def yield_gff_line(file):
    """Iterator, Yields a single non-comment line of a tsv gff file

    :param file: iterable, open file object or list of lines
    :return: yields a list, a gff line split on \t
    """
    for line in file:
        if line.startswith("#"):
            # Remove the comment line, which is unnecessary
            continue
        line = line.strip()
        yield line.split("\t")


def parse_gff_file(file):
    """Iterate through gff lines and extract genes, transcripts and exons

    :param file: iterable, open file opbject or list of lines
    :return: genes: dict, gene id coupled to list of child transcripts.
             transcripts: dict, transcript id coupled to list of child exons.
             exons: dict, exon id coupled to list of start location, end location and length
    """
    genes = {}
    transcripts = {}
    exons = {}
    for gff in yield_gff_line(file):
        if gff[2] == "gene":
            # Add genes to the gene dict
            gene_id = gff[8].split(";")[1].split(",")[0][7:]
            # If more information is needed, add it here!
            # Init empty list for the insertions of transcripts
            genes[gene_id] = []
        elif gff[2] == "transcript" or gff[2] == "mRNA" or gff[2] == "primary_transcript":
            parent_gene = gff[8].split(";")[1][7:]
            if genes.get(parent_gene) is not None:
                # Only transcripts that belong to coding genes are interesting
                trans_id = gff[8].split(";")[0][3:]
                # Add transcript to gene information
                trans_list = genes.get(parent_gene).append(trans_id)
                genes[parent_gene] = trans_list
                # Init empty list for the insertion of exons
                transcripts[trans_id] = []
        elif gff[2] == "exon":
            parent_trans = gff[8].split(";")[1][7:]
            if transcripts.get(parent_trans) is not None:
                # Only transcripts that belong to coding genes are interesting?
                exon_id = gff[8].split(";")[0][3:]
                # Add transcript to gene information
                exon_list = transcripts.get(parent_trans).append(exon_id)
                transcripts[parent_trans] = exon_list
                # Insert exon information
                ### More information can be added here if needed in analyses ###
                exons[exon_id] = (int(gff[3]), int(gff[4]))

    ### Double exons should be cleaned up in analyses ###
    return genes, transcripts, exons


def main(filename="./data/GRCh38_latest_genomic.gff"):
    with open(filename) as file:
        genes, transcripts, exons = parse_gff_file(file)
        print(len(genes))
        print(len(transcripts))
        print(len(exons))


if __name__ == "__main__":
    main()
