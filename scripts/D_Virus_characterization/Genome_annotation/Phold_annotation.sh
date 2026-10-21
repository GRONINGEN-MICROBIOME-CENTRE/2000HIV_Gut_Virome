#!/bin/bash
#SBATCH --job-name=Phold_annotation
#SBATCH --output=Phold_annotation_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

vOTU_dir=$1 #path to directory with FASTA files with viral sequences
vOTU_list=$2 #path to file with the list of vOTU files
output_dir=$3 # output directory
vOTU_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" "${vOTU_list}" | sed 's/\.fa$//')

echo '-------------------- RUNNING Phold on vOTU '${vOTU_ID}' --------------------'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
conda activate /scratch/hb-llnext/conda_envs/Phold/; conda list 

mkdir -p $output_dir

# Run Phold
phold run \
	-i ${vOTU_dir}/$vOTU_ID.fa \
	-o ${output_dir}/$vOTU_ID \
	-p $vOTU_ID \
	-t ${SLURM_CPUS_PER_TASK}

# DB default: /scratch/hb-llnext/conda_envs/Phold/lib/python3.13/site-packages/phold/database
