#!/bin/bash
#SBATCH --job-name=iPHOP_host_assignment
#SBATCH --output=iPHOP_host_assignment_%A_%a.out
#SBATCH --time=02:59:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_dir=$1 #path to directory with the split FASTA files with predicted viral contigs
file_list=${contig_dir}/$2 #file with the list of all split FASTA filenames in the directory
contig_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file_list})
contig_file_name=$(basename "$contig_file" | cut -d . -f1)

echo '-------------------- WORKING WITH '${contig_file}' FILE --------------------'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/iphop133_env; conda list

export PATH=/scratch/hb-llnext/conda_envs/iphop133_env/bin:$PATH
echo "Python: $(which python)"
python -c "import pandas; print('Pandas:', pandas.__file__)"

mkdir -p RESULTS ${TMPDIR}/${contig_file_name}

iphop predict \
	--fa_file ${contig_dir}/${contig_file} \
	--db_dir /scratch/hb-llnext/databases/iphop133_db/Aug_2023_pub_rw \
	--out_dir ${TMPDIR}/${contig_file_name} \
	--num_threads ${SLURM_CPUS_PER_TASK}

mv ${TMPDIR}/${contig_file_name}/Host_prediction_to_genome_m90.csv RESULTS/Host_prediction_to_genome_m90_${contig_file_name}.csv
mv ${TMPDIR}/${contig_file_name}/Host_prediction_to_genus_m90.csv RESULTS/Host_prediction_to_genus_m90_${contig_file_name}.csv
mv ${TMPDIR}/${contig_file_name}/Detailed_output_by_tool.csv RESULTS/Detailed_output_by_tool_${contig_file_name}.csv

# To merge results
# awk 'FNR==1 && NR!=1 {next} {print}' Host_prediction_to_genome_* > Host_prediction_to_genome_all.txt
# awk 'FNR==1 && NR!=1 {next} {print}' Host_prediction_to_genus_* > Host_prediction_to_genus_all.txt
