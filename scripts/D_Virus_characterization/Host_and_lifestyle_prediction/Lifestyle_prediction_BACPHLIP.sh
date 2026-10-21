#!/bin/bash
#SBATCH --job-name=BACPHLIP_present_vOTUs
#SBATCH --output=BACPHLIP_present_vOTUs.out
#SBATCH --mem=90gb
#SBATCH --time=1-0
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

vOTU_seqs=$1 #FASTA file of representative sequences of vOTUs present in the samples

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/BACPHLIP_env; conda list

# Run BACPHLIP
bacphlip \
    -i $vOTU_seqs \
    --multi_fasta
conda deactivate
