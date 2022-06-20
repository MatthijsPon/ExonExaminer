GFF_DIR = "data/ENSEMBL"
SPECIES = ["mouse", "human", "cat", "dog", "bonobo", "chimpanzee", "c.elegans", "zebrafish", "chicken", "arabidopsis_thaliana"]
SPECIES_TEMP = ["mouse", "human"]
OUTDIR = "data/out"

rule all:
    input:
        expand("output/exon_incorporation/{species}/statistical_information.txt", species=SPECIES),
        expand("output/exon_gc/{species}_gc_exons.bed", species=SPECIES_TEMP)


rule exon_incorporation_script:
    input:
        gff3="input/ENSEMBL/{species}.gff3",
    output:
        "output/exon_incorporation/{species}/statistical_information.txt",
    params:
        out_dir="output/exon_incorporation/{species}/",
        sizes="0 100 250 500 1000 2500",
        temp_dir="output/exon_incorporation/{species}/"
    threads: 1
    shell:
        "python3 scripts/exon_incorporation.py {input.gff3} {params.out_dir} --temp_dir {params.temp_dir} {params.sizes}"


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
        "python3 scripts/exons_2_bed.py {input} {output}"


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
