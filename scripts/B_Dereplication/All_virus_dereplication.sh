#!/bin/bash
#SBATCH --job-name=All_virus_dereplication
#SBATCH --output=All_virus_dereplication.out
#SBATCH --mem=20gb
#SBATCH --time=04:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with all viral sequences
blast_db=$2 #path to blast+ database generated from all viral sequences

mkdir -p DEREPLICATION_RESULTS

# Clean environment, load modules
module purge; ml Python/3.10.8-GCCcore-12.2.0 CheckV seqtk; module list

echo -e '\n-------------------- DEREPLICATING VIRAL GENOMES INTO vOTUs --------------------'

#Using the blast+ database as input, calculate pairwise ANI by combining local alignments between sequence pairs:
python anicalc.py \
    -i BLASTN_RESULTS/viruses_blast.tsv  \
    -o DEREPLICATION_RESULTS/viruses_ani.tsv
    
# Perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
python aniclust.py \
    --fna $contig_file \
    --ani DEREPLICATION_RESULTS/viruses_ani.tsv \
    --out DEREPLICATION_RESULTS/viral_clusters.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0

# Generate a TXT and a FASTA file with the representative sequences of each vOTU
cat DEREPLICATION_RESULTS/viral_clusters.tsv | cut -f1 > rep_seqs_vOTUs.txt
seqtk subseq $contig_file  rep_seqs_vOTUs.txt > rep_seqs_vOTUs.fa

# Set permissions 
chmod 440 rep_seqs_vOTUs.txt rep_seqs_vOTUs.fa
