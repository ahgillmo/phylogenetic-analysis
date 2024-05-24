#!/usr/bin/bash
#SBATCH --job-name=calc_varReadProb
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Read each line from the input file
while read -r line; do 

    # Extract the variant ID, alt reads, total reads, major CN, and minor CN
    VarID=$(echo $line | sed 's/ /\t/g' | cut -f 1)
    AltReads=$(echo $line | sed 's/ /\t/g' | cut -f 2)
    TotalReads=$(echo $line | sed 's/ /\t/g' | cut -f 3)
    MajorCN=$(echo $line | sed 's/ /\t/g' | cut -f 4)
    MinorCN=$(echo $line | sed 's/ /\t/g' | cut -f 5)

    # Initialize VarReadProb to 0
    VarReadProb=0

    # Calculate IfGain and NotGain probabilities
    IfGain=$(echo $line | sed 's/ /\t/g' | cut -f 4,5 | awk '{print $1/($1+$2)}')
    NotGain=$(echo $line | sed 's/ /\t/g' | cut -f 4,5 | awk '{print $2/($1+$2)}')

    # Normalize read count by total CN
    NormalizedReadCount_By_TotalCN=$(echo $line | sed 's/ /\t/g' | cut -f 3,4,5 | awk '{print $1/($2+$3)}')

    # Determine the variant read probability based on the CN values and read counts
    if [ $MajorCN == "1" ] && [ $MinorCN == "1" ]; then
        VarReadProb=0.5
    elif (( $(echo "$AltReads < $NormalizedReadCount_By_TotalCN" | bc -l) )); then
        VarReadProb=$NotGain
    elif (( $(echo "$AltReads >= $NormalizedReadCount_By_TotalCN" | bc -l) )); then
        VarReadProb=$IfGain
    else
        echo "Unexpected case encountered."
    fi

    # Output the result
    echo -e $VarID'\t'$AltReads'\t'$TotalReads'\t'$VarReadProb

done < $1 | sort > $1".calcvarreadprob.txt"
