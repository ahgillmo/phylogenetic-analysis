#!/usr/bin/bash
#SBATCH --job-name=S1_PhasingBarcodes
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo -e "The 10x-phased (somatic) vcf is:"'\t'$1
echo -e "The Somatic variants vcf is:"'\t'$2
echo -e "The patientID/SampleID is:"'\t'$3
echo -e "The working directory is:"'\t'$PWD

# Generate the location of somatic variants
bioawk -tc vcf '$filter == "PASS"' $2 | awk '{print $1"\t"$2}' > $3"_somaticvariants.list"

# Reduce the phase variants based on somatic variants
bcftools view -f "PASS","10X_ALLELE_FRACTION_FILTER","10X_PHASING_INCONSISTENT" -Ov -R $3"_somaticvariants.list" -o $1"_10Xphased_somatic.vcf" $1 &&

# Break the somatic phased vcf into chromosome fragments
less -S $1"_10Xphased_somatic.vcf" | egrep -v "#" | awk -F\\\t -v SampleID=$3 '{print>SampleID"_"$1".txt"}'

# Process each chromosome file to extract important phasing information
for file in $3"_chr"*".txt"; do bash /home/ahgillmo/master_scripts_slurm/nuclearBarcoding_scripts/phaseblock_scripts/chromosome_processing_step1.sh $file ; done

# Recombine chromosome phased variants and gather only phase blocks with more than 20 variants
cat *_processed.txt | awk 'NF==7' | sed 's/ /\t/g' | awk -F "\t" '{count[$6]++}END{for (i in count) if (count[i] > 20) {print i}}' | egrep -wv 1 > $3"_phaseblockFile.txt"

# Join processed variants by phase blocks with 20 or more variants, producing phase block files
join -i -1 6 -2 1 <(cat *_processed.txt | awk 'NF==7' | sed 's/ /\t/g' | sort -fk 6) <( less $3"_phaseblockFile.txt" | egrep -v "\-1" | sort) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$7}' | awk '{print >> $6".bc.txt"; close ($6)}'

source activate linkedread_env

# Process phase block files to generate variant-specific phase information
for bc in *.bc.txt ; do ~/miniconda3/envs/linkedread_env/bin/R --slave --args $bc < /home/ahgillmo/master_scripts_slurm/nuclearBarcoding_scripts/phaseblock_scripts/phasingMethod_perblock.r ; done

# Clean up temporary files
# rm *.bc.txt
# rm *_chr*.txt
