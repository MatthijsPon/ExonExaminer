SPECIES = ["mouse", "human", "cat", "dog", "bonobo", "chimpanzee", "c.elegans", "zebrafish", "chicken", "arabidopsis_thaliana"]
SPECIES_TEMP = ["mouse", "human"]
CDS_SIZES = ["(0,50)", "(50,300)", "(300,10000)"]

rule all:
    input:
        expand("output/exon_incorporation/{species}/statistical_information.txt", species=SPECIES),
        expand("output/exon_gc/{species}_gc_exons.bed", species=SPECIES_TEMP),
        "output/exon_gc/full_gc_analysis.txt",
        expand("output/codon_usage/human/hbar_graph_size_{cds}.png", cds=CDS_SIZES)


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
        "python3 scripts/full_gc_analysis.py {output} {input.gc} {input.usage} && touch {output}"


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


rule cds_divide_into_groups:
    input:
        fa="output/codon_usage/{species}_cds_seq_phase_aware.fa",
        pickle="output/parsed_gff/{species}.pickle"
    output:
        ["output/codon_usage/{species}/group_" + cds + ".fa" for cds in CDS_SIZES]
    params:
        cds=CDS_SIZES,
        outdir="output/codon_usage/{species}/"
    shell:
        "python3 scripts/cds_split_fa_2_size_groups.py {input.fa} {input.pickle} --out_dir {params.outdir} {params.cds}"

rule cusp:
    input:
        "output/codon_usage/{species}/group_{cds}.fa"
    output:
        "output/codon_usage/{species}/codon_usage_group_{cds}.fa"
    params:
        cusp_install="/exports/humgen/mnpon/EMBOSS/bin/cusp"
    shell:
        "{params.cusp_install} {input} {output}"

rule analyze_codon_usage:
    input:
        ["output/codon_usage/{species}/codon_usage_group_" + cds + ".fa" for cds in CDS_SIZES]
    output:
        ["output/codon_usage/{species}/hbar_graph_size_" + cds + ".png" for cds in CDS_SIZES]
    params:
        out_dir="output/codon_usage/{species}/"
    shell:
        "python3 scripts/codon_usage.py --out_dir {params.out_dir} {input}"

