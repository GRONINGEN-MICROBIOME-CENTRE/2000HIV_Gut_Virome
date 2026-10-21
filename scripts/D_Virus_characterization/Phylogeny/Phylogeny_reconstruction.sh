#!/bin/bash
#SBATCH --job-name=Final_alignments_tree
#SBATCH --output=Final_alignments_tree.out
#SBATCH --mem=64gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

msa_dir=$1        # Directory with PFAM MSAs (*.aln)
output_dir=$2     # Output directory

mkdir -p "$output_dir/trimmed"
mkdir -p "$output_dir/concat"

# Activate conda environment containing ClipKIT, AMAS, FastTree
# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /home2/p304845/Conda_envs/Fasttree_conda/; conda list

# 1. Trim MSAs with ClipKIT
echo "=== Trimming MSAs ==="

for msa in "$msa_dir"/*.aln; do
    base=$(basename "$msa" .aln)
    echo "Trimming $base"
    clipkit "$msa" -o "$output_dir/trimmed/${base}.aln"
done

# 2. Concatenate MSAs with AMAS
echo "=== Concatenating MSAs with AMAS ==="

python $CONDA_PREFIX/lib/python3.11/site-packages/amas/AMAS.py concat \
    -i "$output_dir"/trimmed/*.aln \
    -f fasta \
    -d aa \
    -u fasta \
    -t "$output_dir/concat/concatenated.aln" \
    -p "$output_dir/concat/partitions.txt"

# 3. Filter sequences (>90% gaps)
echo "=== Filtering sequences with >90% gaps ==="

filtered="$output_dir/concat/concatenated.filtered.aln"

awk '
BEGIN { RS=">"; FS="\n" }
NR>1 {
    header=$1
    seq=""
    for(i=2;i<=NF;i++) seq=seq $i

    tmp=seq
    gsub("-", "", tmp)

    if(length(tmp)/length(seq) >= 0.1) {
        print ">"header
        for(i=2;i<=NF;i++) print $i
    }
}' "$output_dir/concat/concatenated.aln" > "$filtered"

# 4. Build phylogenetic tree with FastTree
echo "=== Building phylogenetic tree ==="

fasttree -wag -gamma "$filtered" > "$output_dir/tree.nwk"

echo "=== DONE ==="
