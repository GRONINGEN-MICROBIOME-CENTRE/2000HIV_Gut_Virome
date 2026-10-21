#!/bin/bash
#SBATCH --job-name=inStrain_compare
#SBATCH --output=inStrain_compare_%A_%a.out
#SBATCH --mem=70gb
#SBATCH --time=3-0
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

IS_profiles=$1      # directory with inStrain profiles
vOTU_list=$2        # file with vOTU IDs
results_dir=$3      # output directory

vOTU_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${vOTU_list})

echo "---- Running inStrain compare for ${vOTU_ID} ----"

# Load environment
ml purge; module load Python/3.10.8-GCCcore-12.2.0
source /scratch/hb-llnext/python_venvs/inStrain/bin/activate

# Temp dir
TMP_VOTU_DIR="${TMPDIR}/${vOTU_ID}"
mkdir -p "${TMP_VOTU_DIR}"

# Get profiles with breadth ≥ 0.74
grep -rH -w ${vOTU_ID} $IS_profiles --include="*inStrain_scaffold_info.tsv" \
  | awk -F'\t' '$4 > 0.74 {print $1}' \
  | cut -d ':' -f1 \
  | sed 's|\(.*_inStrain\).*|\1|' \
  | sort -u > "${TMP_VOTU_DIR}/profiles.txt"

profile_count=$(wc -l < "${TMP_VOTU_DIR}/profiles.txt")

# Require at least 2 samples
if (( profile_count < 2 )); then
  echo "${vOTU_ID}" >> "Skipped_vOTUs_few_profiles.txt"
  echo "Only ${profile_count} profiles found for ${vOTU_ID}, skipping."
  exit 0
fi

# Extract sample names 
sed 's|.*/||' "${TMP_VOTU_DIR}/profiles.txt" > "${TMP_VOTU_DIR}/samples.txt"

# Output folder
mkdir -p "${TMP_VOTU_DIR}/${vOTU_ID}_output"
echo "${vOTU_ID}" > "${TMP_VOTU_DIR}/${vOTU_ID}.txt"

# All-vs-all pairwise comparisons
samples=($(cat "${TMP_VOTU_DIR}/samples.txt"))
total=${#samples[@]}
total_comparisons=$(( total * (total - 1) / 2 ))
counter=0

echo "Total comparisons to run: ${total_comparisons}"

for ((i=0; i<${#samples[@]}; i++)); do
  for ((j=i+1; j<${#samples[@]}; j++)); do

    s1=${samples[$i]}
    s2=${samples[$j]}

    ((counter++))
    echo "[${counter}/${total_comparisons}] Comparing $s1 vs $s2 for $vOTU_ID"

    outdir="${TMP_VOTU_DIR}/${vOTU_ID}_output/${s1}_${s2}"
    mkdir -p "$outdir"

    inStrain compare \
      -i ${IS_profiles}/${s1} ${IS_profiles}/${s2} \
      -o "$outdir" \
      --scaffolds "${TMP_VOTU_DIR}/${vOTU_ID}.txt" \
      --min_cov 5 \
      --min_freq 0.05 \
      --fdr 1e-06 \
      -p ${SLURM_CPUS_PER_TASK}

  done
done

# Merge results
comp_files=$(find "${TMP_VOTU_DIR}/${vOTU_ID}_output/" -type f -name '*comparisonsTable.tsv')

if [[ -n "$comp_files" ]]; then
    mkdir -p "${results_dir}"
    out_file="${results_dir}/${vOTU_ID}_comparisons.tsv"

    first=1
    for f in $comp_files; do
        if [[ $first -eq 1 ]]; then
            cat "$f" > "$out_file"
            first=0
        else
            tail -n +2 "$f" >> "$out_file"
        fi
    done

    echo "Results written to $out_file"
else
    echo "No comparison tables found for ${vOTU_ID}"
fi

# Clean up
rm -rf "${TMP_VOTU_DIR}"

echo "---- DONE for ${vOTU_ID} ----"
