#!/bin/bash
#SBATCH --job-name=HMM_search_markers
#SBATCH --output=HMM_search_markers.out
#SBATCH --mem=64gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

hmm_dir=$1          # directory with *.h3m HMMs
proteins=$2         # protein fasta (HIV2000 proteins)
output_dir=$3       # output directory
tooldir="/scratch/p304845//HIV2000/2_PHYLOGENY/diversify_pfam"

mkdir -p "$output_dir"

# TerL markers (a lower threshold will be used)
terl_markers="PF03237 PF03354 PF17289 PF20441 PF20454"

# Load HMMER and SeqKit
module purge; ml HMMER SeqKit; ml list

# Loop over each HMM
for hmm in "$hmm_dir"/*.h3m; do

    pfam=$(basename "$hmm" | cut -d'_' -f1)

    echo "Processing HMM: $pfam"

    # Set score threshold
    threshold=25
    for terl in $terl_markers; do
        if [[ "$pfam" == "$terl" ]]; then
            threshold=10
        fi
    done

    # Output files
    tblout="$output_dir/${pfam}.tblout"
    domtblout="$output_dir/${pfam}.domtblout"
    fasta_out="$output_dir/${pfam}.faa"

    # Run hmmsearch
    hmmsearch \
        --cpu $SLURM_CPUS_PER_TASK \
        --tblout "$tblout" \
        --domtblout "$domtblout" \
        "$hmm" \
        "$proteins" > /dev/null

    # Filter: coverage ≥ 60% + score threshold
    awk -v thr="$threshold" '
    BEGIN { OFS="\t" }
    !/^#/ {
        score = $8
        hmm_from = $16
        hmm_to = $17
        hmm_len = $6

        coverage = (hmm_to - hmm_from + 1) / hmm_len

        if (score >= thr && coverage >= 0.6) {
            print $1
        }
    }' "$domtblout" | sort | uniq > "$output_dir/${pfam}_hits.txt"

    # Extract sequences
    seqkit grep -f "$output_dir/${pfam}_hits.txt" "$proteins" > "$fasta_out"

done
