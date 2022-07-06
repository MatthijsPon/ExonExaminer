#!/bin/bash
#SBATCH --job-name=exonexaminer
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=m.n.pon@lumc.nl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12gb
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/exon_examiner_%j.log

module load tools/miniconda/python3.8/4.9.2
module load genomics/ngs/bedtools2/2.29.1
conda activate exonexaminer

snakemake -j4

