#!/bin/bash
#SBATCH --job-name=Viral_dereplication_processing
#SBATCH --output=Viral_dereplication_processing.out
#SBATCH --mem=40gb
#SBATCH --time=00:19:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

viral_fasta_dir=$1 # directory with all the FASTA files to be used as input for dereplication

# Clean environment, load modules
module purge; ml R; module list

# Concatenate all the viral sequences to be used for dereplication 
cat $viral_fasta_dir/*fna > Dereplication_input_sequences.fa

# Execute the R script
Rscript Viral_dereplication_processing.R Dereplication_input_sequences.fa
