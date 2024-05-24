#!/usr/bin/bash
#SBATCH --job-name=Solo_PhasingBarcodes
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source activate linkedread_env

# Process the phaseblock files to generate variant-specific phase information
for bc in *.bc.txt ; do 
    ~/miniconda3/envs/linkedread_env/bin/R --slave --args $bc < /home/ahgillmo/master_scripts_slurm/nuclearBarcoding_scripts/phaseblock_scripts/phasingMethod_perblock.r ; 
done
