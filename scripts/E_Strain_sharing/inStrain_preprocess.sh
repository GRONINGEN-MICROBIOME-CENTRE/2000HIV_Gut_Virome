#!/bin/bash
#SBATCH --job-name=inStrain_preprocess
#SBATCH --output=inStrain_preprocess.out
#SBATCH --mem=8gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

IS_profiles=$1     # Path to directory with inStrain profile results
vOTU_ID=$2         # vOTU ID (in this case IMG_VR_22092)
Meta_file=$3       # Metadata file with at least: sampleID, sample_type (mother/infant)

echo "-------------------- Preprocessing vOTU ${vOTU_ID} --------------------"

# Activate conda env
ml purge; module load Python/3.10.8-GCCcore-12.2.0; ml list
source /scratch/hb-llnext/python_venvs/inStrain/bin/activate

# Select samples with sufficient breadth for this vOTU
grep -rH -w ${vOTU_ID} $IS_profiles --include="*inStrain_scaffold_info.tsv" | \
awk -F'\t' '$4 > 0.74 {print $1}' | cut -d '/' -f1-8 | sort | uniq > list_profiles_${vOTU_ID}.txt
sed 's|.*/||;s/_inStrain$//' list_profiles_${vOTU_ID}.txt | sort | uniq > list_cleaned_${vOTU_ID}.txt

# Filter metadata to get mother and infant sample IDs present in this vOTU
cut -f1,2 ${Meta_file} | grep -wf list_cleaned_${vOTU_ID}.txt | sort | uniq > filtered_meta_${vOTU_ID}.tsv

# Split into mother and infant sample lists
awk '$2 == "mother" {print $1}' filtered_meta_${vOTU_ID}.tsv > mother_samples_${vOTU_ID}.txt
awk '$2 == "infant" {print $1}' filtered_meta_${vOTU_ID}.tsv > infant_samples_${vOTU_ID}.txt

# Due to the high number of maternal samples, split list of maternal samples into several lists
mkdir -p LISTS
split -d -a 3 -n l/100 "mother_samples_${vOTU_ID}.txt" "LISTS/mother_samples_${vOTU_ID}_part_"
ls LISTS > LISTS/maternal_list.txt
cp infant_samples_${vOTU_ID}.txt LISTS

# Generate a file with the vOTU name to be compared
echo ${vOTU_ID} > ${vOTU_ID}.txt # Generate a file with the vOTU name to be compared

