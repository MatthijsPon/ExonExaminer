#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""
import pandas as pd


def main():
    pickle_file = "./data/temp/exon_incorporation.pickle"
    homolog_file = "./data/GRCh38.106_GRCm39.106_biomart.txt"

    dataframe = pd.read_pickle(pickle_file)
    dataframe.sort_values("exon_size", ascending=False, inplace=True)
    no_genes = 100
    for i in range(no_genes, no_genes + 50):
        out = dataframe.head(i)
        if out["parent_gene"].nunique() == 100:
            break

    genes = list(out["parent_gene"].unique())
    with open("./data/gene_list.txt", "w+") as file:
        for gene in genes:
            file.write("{}\n".format(gene))
    # print(genes)
    # hum_mou = pd.read_csv(homolog_file, sep="\t")
    # hum_mou.rename(columns={"Gene stable ID": "Human gene ID", "Mouse gene stable ID": "Mouse gene ID"}, inplace=True)
    # print(hum_mou)
    # genes = hum_mou[hum_mou["Human gene ID"].isin(genes)]
    # genes.reset_index(inplace=True, drop=True)
    # genes.to_csv("./data/subset_100_hum_mou.tsv", sep="\t", index=False)
    # return None


if __name__ == "__main__":
    main()
