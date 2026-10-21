#!/bin/bash
#SBATCH --job-name=Summary_viral_detection
#SBATCH --output=Summary_viral_detection.out
#SBATCH --mem=4gb
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --export=NONE
#SBATCH --partition=himem

results_folder=$1 #Folder with Viral detection results per sample

ls $results_folder > list_samples.txt
while read i
do
	sample=$(echo $i)
	echo ${sample} | tr '\n' '\t'

# Check geNomad result
	if [[ -f "${results_folder}/${sample}/virome_discovery/geNomad/${sample}_metaspades_contigs_virus_summary.tsv" ]]; then #check if geNomad output file is present
		echo "yes" | tr '\n' '\t'
		awk 'NR>1' ${results_folder}/${sample}/virome_discovery/geNomad/${sample}_metaspades_contigs_virus_summary.tsv | awk -F '\t' '$2 >= 10000' | wc -l 
    	else
      		echo -e 'no' 
    	fi

done < list_samples.txt > Viral_detection_results_summary.txt

# Remove intermediate files
rm list_samples.txt

# Add header to results file
sed -i $'1 i\\\nSample\tgeNomad_file\tgeNomad_Contigs > 10kb' Viral_detection_results_summary.txt
