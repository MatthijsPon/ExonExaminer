GFF_DIR = "data/ENSEMBL"
GFFS = ["mouse", "human", "cat", "dog", "bonobo", "chimpanzee"]
OUTDIR = "data/out"

rule all:
    input:
        expand("{outdir}/{gff}/statistical_information.txt", outdir=OUTDIR, gff=GFFS)


rule exon_incorporation_script:
    input:
        gff3="data/ENSEMBL/{gff}.gff3",
    output:
        "{outdir}/{gff}/statistical_information.txt"
    shell:
        "echo {input.gff3}"
