#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

import argparse as arg
import numpy as np
from helpers import yield_gff_line, determine_gff
import pandas as pd
from numpy.core.defchararray import zfill


def parse_arguments():
    """Parse command line arguments

    :return: parsed command line options
    """
    parser = arg.ArgumentParser(description="Rewrite an biomart ortholog gene tsv to ExOrthist required cluster file.")
    parser.add_argument("infile", type=str,
                        help="Input file")
    parser.add_argument("outfile", type=str,
                        help="Output file")
    parser.add_argument("genome_names", type=str, nargs="+",
                        help="Strings to put as genomes in the output file (in order they appear in the file).")
    args = parser.parse_args()
    return args.infile, args.outfile, args.genome_names


def parse_gff(file):
    """Parse gene ID and alternative names from gff3 files

    :param file: iterable, open gff3 file object
    :return: pandas dataframe object, gff3 ids and names
    """
    data = {
        "ids": [],
        "names": []
    }
    for line in yield_gff_line(file):
        if not line[2] == "gene":
            continue
        data["ids"].append(determine_gff(line[8], "ID=gene:"))
        name = determine_gff(line[8], "Name=")
        if name:
            data["names"].append(name)
        else:
            data["names"].append("")
    return pd.DataFrame(data)


def main():
    """Main function"""
    # Parse args
    infile, outfile, genome_names = parse_arguments()
    # Read in orthologs
    df = pd.read_csv(infile, sep="\t")
    df.dropna(axis=0, inplace=True)
    df.reset_index(drop=True, inplace=True)
    ### Get gene synonyms for both human and mouse ###
    # Normally you would include this in BioMart, but the service was down
    with open("./data/Homo_sapiens.GRCh38.106.chr.gff3") as file:
        gff_human = parse_gff(file)
    with open("./data/Mus_musculus.GRCm39.106.chr.gff3") as file:
        gff_mouse = parse_gff(file)

    gff_human.set_index("ids", inplace=True)
    gff_mouse.set_index("ids", inplace=True)

    df = df.join(gff_human, on="Gene stable ID")
    df.rename(columns={"Gene stable ID": "human ID", "names": "human alt name"}, inplace=True)
    df = df.join(gff_mouse, on="Mouse gene stable ID")
    df.rename(columns={"Mouse gene stable ID": "mouse ID", "names": "mouse alt name"}, inplace=True)
    df.replace('', np.nan, inplace=True)
    df.dropna(inplace=True)
    df.to_csv("./data/temp/hum_mou_no_NA.tsv", columns=["human ID", "human alt name", "mouse ID", "mouse alt name"],
              sep="\t", index=False)
    ### End of gene synonyms ###

    # Output to ExOrthist required file
    with open("./data/temp/hum_mou_no_NA.tsv") as inf, open(outfile, "w+") as outf:
        for idx, line in enumerate(inf):
            if idx == 0:
                continue
            base = "clst{}".format(zfill(str(idx), 6))
            items = line.strip().split("\t")
            if len(items) > 2:
                print(items)
                for idy in range(0, len(items), 2):
                    outf.write(base + "\t{}\t{}\t{}\n".format(genome_names[idy // 2], items[idy], items[idy + 1]))


if __name__ == "__main__":
    main()
