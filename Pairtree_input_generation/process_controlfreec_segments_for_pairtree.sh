#!/usr/bin/bash
#SBATCH --job-name=process_controlfreec
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

####
# Input: The input is the CNV p-value txt from ControlFreeC
####

# Echo the input ControlFreeC segment file
# echo -e "The input ControlFreeC segment file is: $1"

# Merge BED based on the gain
mergeBed -c 4,5,6,7,8,9 -o collapse,distinct,collapse,collapse,collapse,collapse -i <(less -S $1 | egrep -v "uncertainty" | awk '{print "chr"$1"_"$5 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' | sortBed) | \
sed 's/_/\t/g' | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' | sort -V -k1,1 -k2,2 > $1".modifySegment.txt"

# Process each line in the modified segment file
while read -r line; do
    # Extract uncertainty column and sort values
    v1=$(echo $line | sed 's/ /\t/g' | cut -f 7 | sed -e $'s/,/\\\n/g')
    v2=($v1)
    IFS=$'\n' sorted=($(sort <<<"${v2[*]}")); unset IFS
    Best_Certainty=${sorted[0]}
    
    # Find the index of the best certainty value
    Best_Certainty_index=0
    for i in "${!v2[@]}"; do
       if [[ "${v2[$i]}" = "${Best_Certainty}" ]]; then
           Best_Certainty_index=$i
       fi
    done

    # Extract the genotype column and find the most certain genotype
    GT_col=$(echo $line | sed 's/ /\t/g' | cut -f 6 | sed -e $'s/,/\\\n/g')
    GT=($GT_col)
    MostCertain_Genotype=${GT[$Best_Certainty_index]}

    # Extract the Wilcox statistics and find the most certain value
    Wilcox_stats=$(echo $line | sed 's/ /\t/g' | cut -f 8 | sed -e $'s/,/\\\n/g')
    Wilcox=($Wilcox_stats)
    MostCertain_Wilcox=${Wilcox[$Best_Certainty_index]}

    # Extract the Kolmogrov statistics and find the most certain value
    Kolmogrov_stats=$(echo $line | sed 's/ /\t/g' | cut -f 9 | sed -e $'s/,/\\\n/g')
    Kolmogrov=($Kolmogrov_stats)
    MostCertain_Kolmogrov=${Kolmogrov[$Best_Certainty_index]}

    # Extract the raw copy number and find the most certain value
    Raw_copynumber=$(echo $line | sed 's/ /\t/g' | cut -f 4 | sed -e $'s/,/\\\n/g')
    copynumber=($Raw_copynumber)
    MostCertain_copynumber=${copynumber[$Best_Certainty_index]}

    # Append the most certain segments to the output file
    echo -e $line '\t'$Best_Certainty'\t'$MostCertain_Genotype'\t'$MostCertain_Wilcox'\t'$MostCertain_Kolmogrov | \
    sed 's/ /\t/g' | sort -V -k1,1 -k2,2 | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $10 "\t" $12 "\t" $13 "\t" $11}' | \
    awk '$5<=1.0' | awk '$7<=0.05' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8}' | \
    awk -F 'B' '{print $1 "\t" NF-1}' | awk '{print $1 "\t" $2 "\t" $3 "\t" length($6) "\t" $7}' >> $1".MostCertain_Genotype.txt"

done < $1".modifySegment.txt"

# Alternative to bedops to remove duplicated regions
bedtools subtract -a /home/ahgillmo/references/wgs_callingRegions/10X_hg38_primary_sex_chromosomes_centromereRegions_removed.bed -b $1".MostCertain_Genotype.txt" | \
awk '{print $1 "\t" $2 "\t" $3 "\t" 1 "\t" 1}' | grep -v chrY | grep -v chrX | sed 's/ /\t/g' > $1".normal.diploidRegions.txt"

# Output the header for the significant segments file
echo -e Chromosome'\t'Start_position'\t'End_Position'\t'A1_CN'\t'A2_CN > $1".MostCertain_and_Significant_Segments_with_NormalRegionAdded.seg.tsv"

# Combine the normal segments with the most certain CNV modified segments
cat $1".MostCertain_Genotype.txt" $1".normal.diploidRegions.txt" | sort -V -k1,1 -k2,2 | \
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> $1".MostCertain_and_Significant_Segments_with_NormalRegionAdded.seg.tsv"

# Clean up temporary files
rm $1".modifySegment.txt"
rm $1".MostCertain_Genotype.txt"
rm $1".normal.diploidRegions.txt"
