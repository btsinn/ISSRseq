#!/bin/bash

echo "        

ISSRseq -- AnalyzeBAMs
                       
development version 0.2
use help for usage 
    
"

if [[ $1 = help ]]
  then
  
echo "This script takes the output of CreateBAMs and calls and analyzes SNP variants.

Usage is as follows:

REQUIRED:

-O <AssembleReference output prefix -- DO NOT end with /> 
-T <number of parallel threads>
-S <min samples for SNP inclusion in fasta,nexus,phylip matrices>


Dependencies: GATK, vcf2phylip

"

exit 1
fi

while getopts "O:T:M:K:D:S:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        M) MAX_MISSING_PERC=$OPTARG ;;
        K) FAST_K=$OPTARG ;;
        S) MIN_SAMPLES=$OPTARG ;;
       esac
done

#########################################################################

mkdir $PREFIX/gvcfs
mkdir $PREFIX/matrices

REF_DB=$PREFIX/reference/reference_assembly.fa
SAMPLE_LIST=$PREFIX/samples.txt
FAST_STRUC_DIR=/usr/local/src/fastStructure

#########################################################################

#create interval file

grep "^>" $REF_DB | sed 's/^>//' > $PREFIX/reference/list.intervals

#run GATK first step -- HaplotypeCaller

while read -r sample
do

    gatk --java-options "-Xmx100g" HaplotypeCaller --min-base-quality-score 20 -R $REF_DB -L $PREFIX/reference/list.intervals -I $PREFIX/bams/${sample}_sorted_RG.bam -O $PREFIX/gvcfs/${sample}.g.vcf -ERC GVCF >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

done < $SAMPLE_LIST

#run GATK second step -- CombineGVCFs

options=()
for file in $PREFIX/gvcfs/*.g.vcf; do    
    options+=(-V "${file}")    # If you want the full path, leave off ##*/
done

gatk --java-options "-Xmx100g" CombineGVCFs -L $PREFIX/reference/list.intervals -R $REF_DB "${options[@]}" -O $PREFIX/matrices/combined_gvcfs.g.vcf >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

#run GATK third step -- GenotypeGVCFs

gatk --java-options "-Xmx100g" GenotypeGVCFs -L $PREFIX/reference/list.intervals -R $REF_DB -V $PREFIX/matrices/combined_gvcfs.g.vcf -ploidy 2 -O $PREFIX/matrices/raw_SNPs.vcf >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

#select variants from VCF that are characterized by certain attributes

gatk --java-options "-Xmx100g" SelectVariants -R $REF_DB -V $PREFIX/matrices/raw_SNPs.vcf -O $PREFIX/matrices/filtered_SNPs.vcf --select-type-to-exclude INDEL -select-type-to-include SNP --restrict-alleles-to BIALLELIC --selectExpressions "AF > 0.01 && AF < 0.99 && QD > 2.0 && MQ > 40.0 && FS < 60.0 && SOR < 3.0 && MQRankSum > -5.0 && ReadPosRankSum > -4.0" >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

#VCF to phylip, nexus, and fasta formats using python program available from: https://github.com/edgardomortiz/vcf2phylip

vcf2phylip -i $PREFIX/matrices/filtered_SNPs.vcf -m $MIN_SAMPLES -f -n -b >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1
