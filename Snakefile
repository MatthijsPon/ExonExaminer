SPECIES = ["mouse", "human", "cat", "dog", "bonobo", "chimpanzee", "c.elegans", "zebrafish", "chicken", "arabidopsis_thaliana"]
SPECIES_SMALL = ["mouse", "human", "chicken"]
CDS_SIZES = ["0,50", "50,300", "300,10000", "0,48", "49,288", "289,10000", "50,175", "175,300", "300,1000", "1000,5000"]

# Rule to gather all output files from rules
# TODO cleanup this rule
rule all:
    input:
        expand("output/exon_incorporation/{species}/statistical_information.txt", species=SPECIES),
        expand("output/exon_gc/{species}_gc_exons.bed", species=SPECIES_SMALL),
        "output/exon_gc/full_gc_analysis.txt",
        expand("output/codon_usage/{species}/hbar_graph_{cds}.png", species=SPECIES_SMALL, cds=CDS_SIZES),


rule exon_incorporation_pickle:
    input:
        "input/ENSEMBL/{species}.gff3"
    output:
        "output/exon_incorporation/{species}_exon_usage.pickle"
    shell:
        "python3 scripts/exon_incorporation_pickle.py {input} {output}"


rule exon_incorporation_script:
    input:
        "output/exon_incorporation/{species}_exon_usage.pickle"
    output:
        "output/exon_incorporation/{species}/statistical_information.txt"
    params:
        out_dir="output/exon_incorporation/{species}/",
        sizes="0 100 250 500 1000 2500"
    shell:
        "python3 scripts/exon_incorporation_analysis.py {input} {params.out_dir} {params.sizes}"


rule parse_gff3:
    input:
        "input/ENSEMBL/{species}.gff3"
    output:
        "output/parsed_gff/{species}.pickle"
    shell:
        "python3 scripts/parse_gff3.py {input} {output}"


rule gff_statistics:
    input:
        "output/parsed_gff/{species}.pickle"
    output:
        "output/statistics/{species}_exon_statistics.txt"
    params:
        file_pre="{species}",
        out_dir="output/statistics"
    shell:
        "python3 scripts/exon_statistics {input} {params.out_dir} {params.file_pre}"


rule exons_2_bed:
    input:
        "output/parsed_gff/{species}.pickle"
    output:
        "output/exon_gc/{species}_exons.bed"
    shell:
        "python3 scripts/internal_exons_2_bed.py {input} {output}"


rule temp_gunzip_fasta:
    input:
        "input/ENSEMBL/fasta/{species}.fa.gz" 
    output:
        "input/ENSEMBL/fasta/{species}.fa"
    shell:
        "gunzip -c {input} > {output}"


rule calc_gc_exons:
    input:
        fa="input/ENSEMBL/fasta/{species}.fa",
        bed="output/exon_gc/{species}_exons.bed"
    output:
        "output/exon_gc/{species}_gc_exons.bed"
    shell:
        "bedtools nuc -fi {input.fa} -bed {input.bed} > {output} && rm {input.fa}"


# TODO finish this rule and script
rule gc_exons_analysis:
    input:
        gc=expand("output/exon_gc/{species}_gc_exons.bed", species=SPECIES_SMALL),
        usage=expand("output/exon_incorporation/{species}_exon_usage.pickle", species=SPECIES_SMALL)
    output:
        "output/exon_gc/full_gc_analysis.txt"
    shell:
        "python3 scripts/full_gc_analysis.py {output} {input.gc} {input.usage}"


rule cds_2_bed_phase_aware:
    input:
        "input/ENSEMBL/{species}.gff3"
    output:
        "output/codon_usage/{species}_cds_phase_aware.bed"
    shell:
        "python3 scripts/parse_cds_2_bed.py -p {input} {output}"

rule cds_get_fasta_phase_aware:
    input:
        bed="output/codon_usage/{species}_cds_phase_aware.bed",
        fa="input/ENSEMBL/fasta/{species}.fa"
    output:
        "output/codon_usage/{species}_cds_seq_phase_aware.fa"
    shell:
        "bedtools getfasta -s -name -fi {input.fa} -bed {input.bed} > {output}"

rule rename_cds_fasta:
    input:
        fa="output/codon_usage/{species}_cds_seq_phase_aware.fa",
        pickle="output/parsed_gff/{species}.pickle"
    output:
        "output/codon_usage/{species}_cds_seq_renamed.fa"
    shell:
        "python3 scripts/rename_cds_2_exon.py {input.fa} {input.pickle} {output}"

rule cds_divide_into_groups:
    input:
        fa="output/codon_usage/{species}_cds_seq_renamed.fa",
        pickle="output/parsed_gff/{species}.pickle"
    output:
        ["output/codon_usage/{species}/group_" + cds + ".fa" for cds in CDS_SIZES]
    shell:
        "python3 scripts/cds_split_fa_2_size_groups.py {input.fa} {input.pickle} {output}"


rule analyze_codon_usage:
    input:
        "output/codon_usage/{species}/group_{cds}.fa"
    output:
        "output/codon_usage/{species}/hbar_graph_{cds}.png"
    params:
        out_dir="output/codon_usage/{species}/",
        cds="{cds}"
    shell:
        "python3 scripts/codon_usage.py {input} {params.out_dir} {params.cds}"

