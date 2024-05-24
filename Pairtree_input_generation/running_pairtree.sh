#!/usr/bin/bash
#SBATCH --partition=cpu2019,cpu2021,sherlock
##SBATCH --partition=cpu2019-bf05,cpu2017-bf05,cpu2022-bf24,cpu2021-bf24
#SBATCH --mem=100G
#SBATCH --job-name=pairtree_RunTime
##SBATCH --time=5:00:00
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=10

#echo -e "The multisample VCF is:"$1
#echo -e "The PatientID is:"$2


#Need to ensure the starting VCF is index Only required for splitting by -R (region file)
#tabix $1

#<optional for subsetting massive files 200000 + variants> : 
#for sample in `bcftools query -l $1` ; do  bcftools view -R /home/ahgillmo/references/exome_capture_kits/xgen-exome-research-panel-v2-targets-hg38.bed.gz -c1 -Oz -s $sample -o $sample".by.singlesample.vcf.gz" $1 ; done


#Split the multisample VCF into single samples based on sample Name
for sample in `bcftools query -l $1` ; do  bcftools view -c1 -Oz -s $sample -o $sample".by.singlesample.vcf.gz" $1 ; done

#Split the multisample VCF into single samples based on sample Name and remove variants from Xenograft regions that might be due to poor alignment (Mouse homology is high)
#for sample in `bcftools query -l $1` ; do  bcftools view -c1 -Oz -s $sample -T "^"/home/ahgillmo/master_scripts_slurm/pairtree_scripts/Xenograft_commonlyReocurring_variants.bed -o $sample".by.singlesample.vcf.gz" $1 ; done


#Reduce variants by filtering of variants that don't PASS and don't have a minimum depth of 10.
for ssVCF in *.by.singlesample.vcf.gz ; do nn=$(echo $ssVCF | sed 's/.by.singlesample.vcf.gz/.variantfiltered.singlesample.vcf.gz/g') && bcftools filter  -i 'FILTER == "PASS" && FORMAT/DP >= 10' -Oz -o $nn $ssVCF && mv $nn $ssVCF ; done


#Extract the sample names from the multi-sample VCF and form for json parameters file --> Stuck here for some reason
Samples=$(zless $1 | head -n 1000 | grep '#' | egrep "#CHROM" | head -n 1 | cut -f 10- | sed 's/ /\t/g' | sed 's/\t/\n/g' | grep -v BLOOD | grep -v Germline | grep -v BL | grep -v NORMAL | grep -v SM3925_updated_name_LongerRunTime | grep -v CD45 | grep -v SM2819 | grep -v SM2907 | grep -v SM3375 | grep -v SM3787 | grep -v SM3623 | sed 's/_LongerRunTime//g' | grep -v "_N_") 

PerfectSamples=$(echo $Samples | sed 's/ /","/g' | awk '{print "\""$0"\""}')

#Take the single sample mutect2 vcfs and process
for x in $Samples ; do 

echo $x ;
SNV=$(ls $x*.by.singlesample.vcf.gz)
CNV=$(ls $x*.p.value.txt)
#SNV=$(ls $x"_LongerRunTime"*".by.singlesample.vcf.gz") #For SM4218 only
#CNV=$(ls $x"_phased_possorted_bam.bam_CNVs.p.value.txt") #For SM4218 only

bash /home/ahgillmo/master_scripts_slurm/pairtree_scripts/multisample_Mutect_to_pairtree_SSM.sh $SNV $CNV ; 

done

#Combine the data using R or using Join on the command line 
cat ~/master_scripts_slurm/pairtree_scripts/processing_variantReadProbability.r | R --slave --args $PWD $2 && 

#Reconfigure the ssm to match pairtree specifcation (column style)
less $2".ssm" | grep -v NA | sort -Vk 1,1 | cat -n | awk '{print "s"($1-1)"\t"$2"\t"$3"\t"$4"\t"$5}' | sed '1i id name var_reads total_reads var_read_prob' | sed 's/ /\t/g' > $2".complete.ssm"

#create the JSON file using the R script as the order that samples were processed in
json_samples=$(less $2".temporary.json" | sed 's/, /", "/g') &&
echo '{"samples": ['$json_samples'], "clusters": [], "garbage": []}' > $2".complete.json" 

#run pairtree solo
sbatch /home/ahgillmo/master_scripts_slurm/pairtree_scripts/Solo_pairtree_only.sh $2 $3

#Activate the proper script
##source activate pairtree

#JSON format && cluster SNV using pairtree native clustering && run pairtree && plottrees

#Run pairtree cluster-variants
#python3 ~/pairtree/bin/clustervars --parallel 8 $2".complete.ssm" $2".complete.json" $2".filenameout.json"

#Run PairTree
#python3 ~/pairtree/bin/pairtree --params $2".filenameout.json" $2".complete.ssm" $2".results.npz"

#Run plotTree
#python3 ~/pairtree/bin/plottree --runid $2 $2".complete.ssm" $2".filenameout.json" $2".results.npz" $2".results.html"

#Plot unique tree's
#python3 ~/pairtree/bin/summposterior --runid $2 $2".complete.ssm" $2".filenameout.json" $2".results.npz" $2".results.summposterior.html"

#clean up section
#rm $2"_temporary.json"
#rm $2".ssm"

