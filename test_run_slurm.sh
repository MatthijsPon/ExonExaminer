#!/bin/bash
#SBATCH --job-name=mnpon_exonexaminerL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=00:05:00
#SBATCH --output=logs/exon_examiner_%j.log

module load tools/miniconda/python3.8/4.9.2
module load genomics/ngs/bedtools2/2.29.1
conda activate exonexaminer

snakemake -np
