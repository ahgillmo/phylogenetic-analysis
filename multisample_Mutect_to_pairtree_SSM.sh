#!/usr/bin/bash

# SBATCH directives for job scheduling and resource allocation
#SBATCH --partition=single,lattice,parallel
#SBATCH --job-name=pairtree_calc
#SBATCH --mem=9G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Pre Steps
# Step 1: Processing ControlFreec into a format that includes diploid regions 
# (also filters based on p-value and genotype certainty)
# Step 2: Split Multisample Mutect SNV/INDEL VCF by Sample

# Purpose: Process single samples CNV and SNV/INDEL VCF information from Tumor, Xeno, or CellLine.

# SNV/INDEL files are created using the following method:

# Step 1: Split the initial MuTect file into N samples
# For splitting the XXXX.Merged.Filtered.leftnormalized.Mutect2.vcf.gz file
# This loop processes each VCF file and splits it into sample
 for file in *.vcf*; do 
     for sample in `bcftools query -l $file`; do 
         bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
     done 
 done

# Echoing the input files
echo -e "The Mutect2 SNV/INDEL variant file is: $1"
echo -e "The ControlFreec p-value segment file is: $2"

# Process the ControlFreec segments to produce the most certain segments and add in unaffected (diploid) regions
bash process_controlfreec_segments_for_pairtree.sh $2

# Create a variant list with labels "Variant" "ALT-count" "Depth" "CN_A" "CN_B"
bedtools intersect -wa -wb -a <(zless $1 | egrep -v '#' | awk '$7=="PASS"' | sed 's/ /\t/g' | awk '{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$10}') -b <(zless $2".MostCertain_and_Significant_Segments_with_NormalRegionAdded.seg.tsv" | egrep -v Start | sed 's/ /\t/g') | \
awk '{print $1"_"$2"_"$4"_"$5"\t"$10"\t"$11"\t"$6}' | sed 's/:/\t/g' | \
awk '{print $1"\t"$5"\t"$2"\t"$3}' | sed 's/,/\t/g' | \
awk '{print $1"\t"$3"\t"($3+$2)"\t"$4"\t"$5}' | sort -u > $1".temp.PairTree.txt"

# Calculate the variant read probability required for PairTree
bash Calc_read_probability.sh $1".temp.PairTree.txt"

# Clean up the temporary files
rm $2".MostCertain_and_Significant_Segments_with_NormalRegionAdded.seg.tsv"
