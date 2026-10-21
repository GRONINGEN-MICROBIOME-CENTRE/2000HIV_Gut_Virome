#!/bin/bash
#SBATCH --job-name=Muscle_alignment_markers
#SBATCH --output=Muscle_alignment_markers.out
#SBATCH --mem=64gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

hits=$1          # HMM_best_hits.tsv
proteins=$2      # protein fasta
output_dir=$3    # output directory

mkdir -p "$output_dir/MSAs"
mkdir -p "$output_dir/tmp"

# Load  MUSCLE and Seqkit
module purge; ml MUSCLE Seqkit; ml list

# Get unique PFAMs
cut -f2 "$hits" | tail -n +2 | sort | uniq > "$output_dir/pfam_list.txt"

# Loop over PFAMs
while read -r pfam; do

    echo "Processing PFAM: $pfam"

    tmp_regions="$output_dir/tmp/${pfam}_regions.txt"
    trimmed_fasta="$output_dir/tmp/${pfam}_trimmed.faa"
    msa_raw="$output_dir/tmp/${pfam}_raw.aln"
    msa_out="$output_dir/MSAs/${pfam}.aln"

    # Extract protein + coordinates (ali_from, ali_to)
    awk -v pf="$pfam" 'BEGIN{OFS="\t"} NR>1 && $2==pf {print $1, int($6), int($7)}' "$hits" > "$tmp_regions"

    # Clear previous file
    > "$trimmed_fasta"

    # Extract subsequences
    while read -r id start end; do
        echo "  Extracting $id:$start-$end"

        seqkit grep -p "$id" "$proteins" | \
        seqkit subseq -r ${start}:${end} >> "$trimmed_fasta"

    done < "$tmp_regions"

    # Run MUSCLE
    muscle -super5 "$trimmed_fasta" -output "$msa_raw"

    # Fix headers: keep only genome names (instead of proteins)
    sed -E 's/^>([^ ]+).*/>\1/; s/_[^_]+$//' "$msa_raw" > "$msa_out"

done < "$output_dir/pfam_list.txt"
