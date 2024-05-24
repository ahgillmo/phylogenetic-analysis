#!/usr/bin/bash
#SBATCH --job-name=PhasingBam_and_Barcodes
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo "Initial input is $1"

while read -r line ; 
do 
    # Extract basic information
    basics=$(echo $line | sed 's/ /\t/g' | cut -f 1,2,4,5) 
    # Extract AO (Alternate Observations), RO (Reference Observations), and DP (Depth)
    ao=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f 2 | sed 's/,/\t/g' | cut -f 2)
    ro=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f 2 | sed 's/,/\t/g' | cut -f 1)
    dp=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f 3)
    # Extract block data
    blockdat=$(echo $line | sed 's/ /\t/g' | cut -f 9 | sed 's/:/\t/g')
    # Calculate PS (Phase Set) and BX (Barcode Index)
    PS=$(echo ${blockdat[@]/PS//} | cut -d/ -f1 | wc -w | tr -d ' ' | awk '{print$1+1}')
    BX=$(echo ${blockdat[@]/BX//} | cut -d/ -f1 | wc -w | tr -d ' ' | awk '{print$1+1'})
    # Extract barcode information
    BARCODE_All=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f $BX)
    BARCODE_Ref=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f $BX | sed 's/,/\t/g' | cut -f 1)
    BARCODE_Alt=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f $BX | sed 's/,/\t/g' | cut -f 2)
    # Extract phase block and genotype
    PHASEBLOCK=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f $PS)
    GENTY=$(echo $line | sed 's/ /\t/g' | cut -f 10 | cut -d ":" -f 1)
    
    # Print processed information
    echo $basics $ro $ao $dp $GENTY $PHASEBLOCK $BARCODE_All | awk '{print $1"_"$2"_"$3"_"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ; 
done < $1 > $1"_processed.txt"
