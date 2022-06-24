SPECIES = ["mouse", "human", "cat", "dog", "bonobo", "chimpanzee", "c.elegans", "zebrafish", "chicken", "arabidopsis_thaliana"]
SPECIES_TEMP = ["mouse", "human"]

rule all:
    input:
        expand("output/exon_incorporation/{species}/statistical_information.txt", species=SPECIES),
        expand("output/exon_gc/{species}_gc_exons.bed", species=SPECIES_TEMP),
        "output/exon_gc/full_gc_analysis.txt"


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


rule gc_exons_analysis:
    input:
        gc=expand("output/exon_gc/{species}_gc_exons.bed", species=SPECIES_TEMP),
        usage=expand("output/exon_incorporation/{species}_exon_usage.pickle", species=SPECIES_TEMP)
    output:
        "output/exon_gc/full_gc_analysis.txt"
    shell:
        "python3 scripts/gc_analysis.py {output} {input.gc} {input.usage}"


rule cds_2_bed:
    input:
        "input/ENSEMBL/{species}.gff3"
    output:
        "output/codon_usage/{species}_cds.bed"
    shell:
        "python3 scripts/parse_cds_2_bed.py {input} {output}"

rule cds_get_fasta:
    input:
        bed="output/codon_usage/{species}_cds.bed",
        fa="input/ENSEMBL/fasta/{species}.fa"
    output:
        "output/codon_usage/{species}_cds_seq.bed"
    shell:
        "bedtools getfasta -s -bedOut -fi {input.fa} -bed {input.bed} -fo {output}"


rule merge_cds_sequences:
    input:
        "output/codon_usage/{species}_cds_seq.bed"
    output:
        "output/codon_usage/{species}_cds.fa"
    rule:
        "python3 scripts/merge_cds_sequences.py {input} {output}"

