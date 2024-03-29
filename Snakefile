SPECIES = ["mouse", "human", "cat", "dog", "bonobo", "chimpanzee", "c.elegans", "zebrafish", "chicken", "arabidopsis_thaliana"]
SPECIES_SMALL = ["mouse", "human", "chicken"]
CSD_COMPARISON = ["49,259_vs_0,48", "49,259_vs_260,1000000", "69,198_vs_0,68", "69,198_vs_199,100000"]

# Rule to gather all output files from rules
rule all:
    input:
        # Gather exon statistics for all species
        expand("output/statistics/{species}/exon_statistics.txt", species=SPECIES),
        # Gather exon incorporation for all species
        expand("output/exon_incorporation/{species}/cumulative_barplot.png", species=SPECIES),
        # Gather codon usage of species small
        expand("output/codon_usage/{species}/aa_diff_hbar_{cds}.png",species=SPECIES_SMALL,cds=CSD_COMPARISON),
        expand("output/codon_usage/{species}/aa_over_size_all_types.png", species=SPECIES_SMALL),
        # Gather GC content of all species
        expand("output/exon_gc/{species}/scatter_roll_avg_20_by_size.png", species=SPECIES_SMALL),
        # Get Expression data for human
        "output/expression/human/scatterplot_low_usage.png",


rule parse_gff3:
    input:
        "input/ENSEMBL/{species}.gff3"
    output:
        "output/parsed_gff/{species}.pickle"
    shell:
        "python3 scripts/parse_gff3.py {input} {output}"

rule gff_features:
    input:
        "output/parsed_gff/{species}.pickle"
    output:
        "output/statistics/{species}/exon_statistics.txt",
        "output/statistics/{species}/5_95_quartiles.txt"
    params:
        out_dir="output/statistics/{species}/"
    shell:
        "python3 scripts/gff_features.py {input} {params.out_dir}"


rule exon_incorporation_pickle:
    input:
        "input/ENSEMBL/{species}.gff3"
    output:
        "output/exon_incorporation/{species}_exon_usage.pickle"
    shell:
        "python3 scripts/exon_incorporation_pickle.py {input} {output}"

rule exon_incorporation_script:
    input:
        pickle="output/exon_incorporation/{species}_exon_usage.pickle",
        quartiles="output/statistics/{species}/5_95_quartiles.txt"
    output:
        "output/exon_incorporation/{species}/cumulative_barplot.png",
        "output/exon_incorporation/{species}/exon_incorporation_positive_borders.txt",
        "output/exon_incorporation/{species}/genes_ratio_smaller_0.csv",
        "output/exon_incorporation/{species}/genes_ratio_bigger_0.csv",
    params:
        out_dir="output/exon_incorporation/{species}/",
    shell:
        "python3 scripts/exon_incorporation_analysis.py {input.pickle} {input.quartiles} {params.out_dir}"


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
        "bedtools nuc -fi {input.fa} -bed {input.bed} > {output}"

rule gc_exons_analysis:
    input:
        pickle="output/exon_incorporation/{species}_exon_usage.pickle",
        bed="output/exon_gc/{species}_gc_exons.bed",
        quartiles="output/statistics/{species}/5_95_quartiles.txt",
        inc_borders="output/exon_incorporation/{species}/exon_incorporation_positive_borders.txt"
    output:
        "output/exon_gc/{species}/scatter_roll_avg_20_by_size.png"
    params:
        outdir="output/exon_gc/{species}/"
    shell:
        "python3 scripts/gc_analysis.py {input.pickle} {input.bed} "
        "{input.quartiles} {input.inc_borders} {params.outdir}"


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

rule codon_usage_over_size:
    input:
        fa="output/codon_usage/{species}_cds_seq_renamed.fa",
        pickle="output/parsed_gff/{species}.pickle"
    output:
        "output/codon_usage/{species}/aa_over_size_all_types.png"
    params:
        outdir="output/codon_usage/{species}/"
    shell:
        "python3 scripts/codon_usage_over_size.py {input.fa} {input.pickle} {params.outdir}"


rule cds_divide_into_groups:
    input:
        fa="output/codon_usage/{species}_cds_seq_renamed.fa",
        pickle="output/parsed_gff/{species}.pickle"
    output:
        "output/codon_usage/{species}/group_{cds}.fa"
    shell:
        "python3 scripts/cds_split_fa_2_size_groups.py {input.fa} {input.pickle} {output}"

rule analyze_codon_usage:
    input:
        cds1="output/codon_usage/{species}/group_{cds1}.fa",
        cds2="output/codon_usage/{species}/group_{cds2}.fa",
    output:
        "output/codon_usage/{species}/aa_diff_hbar_{cds1}_vs_{cds2}.png",
        "output/codon_usage/{species}/codon_usage_hbar_{cds1}_vs_{cds2}.png"

    params:
        out_dir="output/codon_usage/{species}/",
        cds="{cds1}_vs_{cds2}"
    shell:
        "python3 scripts/codon_usage.py {input} {params.out_dir} {params.cds}"


rule transcript_expression:
    input:
        gff="input/ENSEMBL/human.gff3",
        trans_tpm="input/GTEx/transcript_tpm.gct.gz",
        partner="input/GTEx/transcript_annotations.txt",
        smaller_genes="output/exon_incorporation/human/genes_ratio_smaller_0.csv",
        bigger_genes="output/exon_incorporation/human/genes_ratio_bigger_0.csv"
    output:
        "output/expression/human/scatterplot_low_usage.png",
        "output/expression/human/scatterplot_high_usage.png"
    params:
        out_dir="output/expression/human/"
    shell:
        "python3 scripts/transcript_expression.py {input.gff} {input.trans_tpm} {input.partner} "
        "{input.smaller_genes} {input.bigger_genes} -o {params.out_dir}"

