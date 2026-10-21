#!/bin/bash
#SBATCH --job-name=HMM_generation_markers_pfam
#SBATCH --output=HMM_generation_markers_pfam.out
#SBATCH --mem=60gb
#SBATCH --time=04:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

markerfile=$1    # TXT file with list of PFAM/NCBIFAM IDs (one per line)
databases=$2     # directory with DBs (e.g., UniRef30, BFD)
output_dir=$3    # output directory
tooldir="/scratch/p304845/HIV2000/2_PHYLOGENY/diversify_pfam"

# Markers to be used:
# TerL ATPase domain:
# PF03237 PF03354
# TerL nuclease domain:
# PF17289 PF20441 PF20454
# Major capsid protein:
# PF05065 PF19821 PF21703
# Portal protein:
# PF04860 PF05133 PF16510

# Crreate a subdirectory to store the HMMs
mkdir -p $output_dir/HMMs

# Loop over each marker ID
while read -r pfam; do

    echo "Processing marker: $pfam"

    # Run diversify-pfam via Pixi
    cd "$tooldir" || exit 1
    
    pixi run python ./diversify-pfam \
        "$pfam" \
        $output_dir \
        "${output_dir}/HMMs/${pfam}_div.h3m" \
        -d "$databases/UniRef30_2023_02" \
        --verbose

done < "$markerfile"
