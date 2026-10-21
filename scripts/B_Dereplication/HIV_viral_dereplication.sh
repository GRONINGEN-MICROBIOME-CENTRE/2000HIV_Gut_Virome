#!/bin/bash
#SBATCH --job-name=HIV_viral_dereplication
#SBATCH --output=HIV_viral_dereplication.out
#SBATCH --mem=20gb
#SBATCH --time=09:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #file with all the renamed viral genomes
blast_db=$2 #path to blast+ database generated from all viral sequences

mkdir -p HIV_DEREPLICATION

# Load SeqKit
module purge; ml SeqKit Python/3.10.8-GCCcore-12.2.0 CheckV seqtk; ml list

# Extract only HIV genomes from the FASTA file with the renamed viral sequences
seqkit grep -r -p "HIV_" $contig_file > HIV_viral_genomes.fa

# Extract only comparisons between HIV genomes 
grep -v "UHGV" BLASTN_RESULTS/viruses_blast.tsv | grep -v "RefSeq" > BLASTN_RESULTS/HIV_viruses.tsv

echo -e '\n-------------------- DEDUPLICATING HIV VIRAL GENOMES --------------------'

# Run deduplication for contigs with CheckV scripts

# Using the blast+ database as input, calculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/p304845/HIV2000/VIROME/3_DEREPLICATION/anicalc.py \
    -i BLASTN_RESULTS/HIV_viruses.tsv  \
    -o HIV_DEREPLICATION/HIV_viruses_ani.tsv

# Perform UCLUST-like clustering (99% ANI + 95% AF):
python /scratch/p304845/HIV2000/VIROME/3_DEREPLICATION/aniclust.py \
    --fna HIV_viral_genomes.fa \
    --ani HIV_DEREPLICATION/HIV_viruses_ani.tsv \
    --out HIV_DEREPLICATION/HIV_viral_genome_clusters.tsv \
    --min_ani 99 \
    --min_tcov 95 \
    --min_qcov 0

# Extract the representatives viral genomes (deduplicated genomes)
cat HIV_DEREPLICATION/HIV_viral_genome_clusters.tsv  | cut -f1 > HIV_dedup_viral_genomes.txt
seqtk subseq HIV_viral_genomes.fa HIV_dedup_viral_genomes.txt > HIV_dedup_viral_genomes.fa

echo -e '\n-------------------- DEREPLICATING HIV VIRAL GENOMES INTO vOTUs --------------------'

# Perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
# Use NON-deduplicated genomes as input
python /scratch/p304845/HIV2000/VIROME/3_DEREPLICATION/aniclust.py \
    --fna HIV_viral_genomes.fa \
    --ani HIV_DEREPLICATION/HIV_viruses_ani.tsv \
    --out HIV_DEREPLICATION/HIV_vOTU_clusters.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0

# Generate a TXT and a FASTA file with the representative sequences of each LLNEXT vOTU (only LLNEXT)
cat HIV_DEREPLICATION/HIV_vOTU_clusters.tsv | cut -f1 > rep_seqs_HIV_vOTUs.txt
seqtk subseq HIV_viral_genomes.fa rep_seqs_HIV_vOTUs.txt > rep_seqs_HIV_vOTUs.fa

# Set permissions 
chmod 440 rep_seqs_HIV_vOTUs.txt rep_seqs_HIV_vOTUs.fa
