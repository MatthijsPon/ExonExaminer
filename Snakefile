GFF_DIR = "data/ENSEMBL"
GFFS = ["mouse", "human", "cat", "dog", "bonobo", "chimpanzee", "c.elegans", "zebrafish", "chicken", "arabidopsis_thaliana"]
OUTDIR = "data/out"

rule all:
    input:
        expand("{outdir}/{gff}/statistical_information.txt", outdir=OUTDIR, gff=GFFS)


rule exon_incorporation_script:
    input:
        gff3="data/ENSEMBL/{gff}.gff3",
    output:
        "{outdir}/{gff}/statistical_information.txt",
    params:
        out_dir="data/out/{gff}/",
        sizes="0 100 250  500 1000 2500",
        temp_dir="data/temp/{gff}"
    threads: 1
    shell:
        "python3 exon_incorporation.py {input.gff3} {params.out_dir} --temp_dir {params.temp_dir} {params.sizes}"
