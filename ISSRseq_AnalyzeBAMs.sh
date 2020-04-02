#!/bin/bash

echo "        

ISSRseq -- AnalyzeBAMs
                       
development version 0.4
use help for usage 
    
"

if [[ $1 = help ]]
  then
  
echo "This script takes the output of CreateBAMs and calls and analyzes SNP variants.

Usage is as follows:

REQUIRED:

-O <AssembleReference output prefix -- DO NOT end with /> 
-T <number of parallel threads>

Dependencies: GATK

0.3 -- 	removed fastStructure
	added minimum base quality for HaplotypeCaller variant scoring
	removed from SelectVariants --restrict-alleles-to BIALLELIC
	removed minor allele filter AF < 0.01 && AF > 0.99 from SelectVariants due to inclusion of multiallelic variants

"

exit 1
fi

while getopts "O:T:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
       esac
done

#########################################################################

mkdir $PREFIX/gvcfs
mkdir $PREFIX/variants

REF_DB=$PREFIX/reference/final_reference_assembly.fa
SAMPLE_LIST=$PREFIX/samples.txt

#########################################################################

#create interval file

grep "^>" $REF_DB | sed 's/^>//' > $PREFIX/reference/list.intervals

#run GATK first step -- HaplotypeCaller

while read -r sample
do

    gatk --java-options "-Xmx100g" HaplotypeCaller --linked-de-bruijn-graph --native-pair-hmm-threads $THREADS --native-pair-hmm-use-double-precision true -R $REF_DB -L $PREFIX/reference/list.intervals -I $PREFIX/bams/${sample}_sorted_RG.bam -O $PREFIX/gvcfs/${sample}.g.vcf -ERC GVCF >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

done < $SAMPLE_LIST

#run GATK second step -- CombineGVCFs

options=()
for file in $PREFIX/gvcfs/*.g.vcf; do    
    options+=(-V "${file}")    # If you want the full path, leave off ##*/
done

gatk --java-options "-Xmx100g" CombineGVCFs -L $PREFIX/reference/list.intervals -R $REF_DB "${options[@]}" -O $PREFIX/variants/combined_gvcfs.g.vcf >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

#run GATK third step -- GenotypeGVCFs

gatk --java-options "-Xmx100g" GenotypeGVCFs -L $PREFIX/reference/list.intervals -R $REF_DB -V $PREFIX/variants/combined_gvcfs.g.vcf -O $PREFIX/variants/raw_SNPs.vcf >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

#select variants from VCF that are characterized by GATK-recommended hard filters

gatk --java-options "-Xmx100g" SelectVariants -R $REF_DB -V $PREFIX/variants/raw_SNPs.vcf -O $PREFIX/variants/filtered_SNPs.vcf --select-type-to-include SNP --selectExpressions "QD > 2.0 && MQ > 40.0 && FS < 60.0 && SOR < 3.0 && MQRankSum > -5.0 && ReadPosRankSum > -4.0" >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1
