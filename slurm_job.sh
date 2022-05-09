#!/bin/bash
#SBATCH --job-name=mnpon_exonexaminer
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=m.n.pon@lumc.nl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12gb
#SBATCH --time=01:00:00
#SBATCH --output=exon_examiner_%j.log

module load tools/miniconda/python3.8/4.9.2
conda activate /home/mnpon/conda/exonexaminer

snakemake -n4

