#!/usr/bin/env python3
"""
Author: Matthijs Pon
Description:
"""

from parse import yield_gff_line


def present_gff_elements(file_object):
    """"""
    elements = set()
    for line in yield_gff_line(file_object):
        elements.add(line[2])
    return elements


def get_children(role, role_connection_dict):
    children = role_connection_dict.get(role)
    options = set()
    if children is None:
        return [role]
    # for item in get_children(role, role_connection_dict):
    #     if item == role or item in :
    #         continue
    #     options.add([role] + get_children(item, role_connection_dict))


def present_gff_elements_hierarchically(file_object):
    """"""
    id_role_and_parent = {}
    full_parents = set()
    role_n_parent_role = {}

    for line in yield_gff_line(file_object):
        element = line[2]
        parent = None
        elem_id = None
        for sub in line[8].split(";"):
            if sub.startswith("ID="):
                elem_id = "".join(sub.split("=")[1:])
            elif sub.startswith("Parent="):
                parent = "".join(sub.split("=")[1:])
        id_role_and_parent[elem_id] = [element, parent]
        if parent is None:
            full_parents.add(element)
        else:
            parent_role = id_role_and_parent.get(parent)[0]
            if not role_n_parent_role.get(parent_role):
                role_n_parent_role[parent_role] = set()
            role_n_parent_role[parent_role].add(element)
    options = set()
    for parent in full_parents:
        for option in get_children(parent, role_n_parent_role):
            options.add(option)
    print(options)


def count_N_X_exon(file):
    exp = 0
    curated = 0
    for line in yield_gff_line(file):
        if line[2] == "exon":
            if line[8].split(";")[0].split("=")[1].startswith("exon-X"):
                exp += 1
            elif line[8].split(";")[0].split("=")[1].startswith("exon-N"):
                curated += 1
    print("exp: {}\ncur: {}".format(exp, curated))


def element_previous(file, string):
    """"""
    type_precursors = set()
    previous = None
    for line in yield_gff_line(file):

        if line[2] == string and previous != string:
            type_precursors.add(previous)

        previous = line[2]
    return type_precursors


def main(gff_file, out_dir="./out/"):
    #
    # # Parse all gff elements
    # with open(gff_file) as file:
    #     elements = present_gff_elements(file)
    # # Write to file
    # with open("{}gff_elements.txt".format(out_dir), "w+") as out_f:
    #     for item in elements:
    #         out_f.write("{}\n".format(item))
    #
    items = {'J_gene_segment', 'lnc_RNA', 'RNase_P_RNA', 'snRNA', 'pseudogene', 'rRNA', 'vault_RNA', 'primary_transcript', 'Y_RNA', 'tRNA', 'scRNA', 'snoRNA', 'transcript', 'V_gene_segment', 'ncRNA', 'antisense_RNA', 'D_gene_segment', 'telomerase_RNA', 'mRNA', 'C_gene_segment', 'miRNA', 'RNase_MRP_RNA'}
    for item in items:
        print(item)
        with open(gff_file) as file:
            print(element_previous(file, item))


if __name__ == "__main__":
    main("../data/GRCh38_latest_genomic.gff")
