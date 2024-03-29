#!/bin/bash

echo "        

ISSRseq -- AnalyzeBAMs
                       
version 1.0
use help for usage 
    
"

set -o errexit
set -o pipefail
set -o nounset

if [[ $1 = help ]]
  then
  
echo "This script takes the output of CreateBAMs and calls and analyzes variants.

Usage is as follows:

REQUIRED:

-O <AssembleReference output prefix -- DO NOT end with /> 
-T <number of parallel threads>
-P <expected ploidy>

Dependencies: GATK
	
	added user input for ploidy expected by HaplotypeCaller


"

exit 1
fi

while getopts "O:T:P:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        P) PLOIDY=$OPTARG ;;
       esac
done

#########################################################################

mkdir $PREFIX/gvcfs
mkdir $PREFIX/variants
mkdir $PREFIX/variants/haplotypecallerBAMs

REF_DB=$PREFIX/reference/final_reference_assembly.fa
SAMPLE_LIST=$PREFIX/samples.txt

#########################################################################

#create interval file

grep "^>" $REF_DB | sed 's/^>//' > $PREFIX/reference/list.intervals

#run GATK first step -- HaplotypeCaller

while read -r sample
do

    gatk --java-options "-Xmx100g" HaplotypeCaller --linked-de-bruijn-graph --native-pair-hmm-threads $THREADS --sample-ploidy $PLOIDY --native-pair-hmm-use-double-precision true -R $REF_DB -L $PREFIX/reference/list.intervals -I $PREFIX/bams/${sample}_sorted_RG_DM.bam -O $PREFIX/gvcfs/${sample}.g.vcf -ERC GVCF -bamout $PREFIX/variants/haplotypecallerBAMs/${sample}_HaplotypeCaller_out.bam >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

done < $SAMPLE_LIST

#run GATK second step -- CombineGVCFs

options=()
for file in $PREFIX/gvcfs/*.g.vcf; do    
    options+=(-V "${file}")    
done

gatk --java-options "-Xmx100g" CombineGVCFs -L $PREFIX/reference/list.intervals -R $REF_DB "${options[@]}" -O $PREFIX/variants/combined_gvcfs.g.vcf >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

#run GATK third step -- GenotypeGVCFs

gatk --java-options "-Xmx100g" GenotypeGVCFs -L $PREFIX/reference/list.intervals -R $REF_DB -V $PREFIX/variants/combined_gvcfs.g.vcf -O $PREFIX/variants/raw_variants.vcf >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

#select variants from VCF that are characterized by GATK-recommended hard filters
#users can modify the number of "--select-type-to-include" flags included below 
#for example, to keep just SNPs, remove "--select-type-to-include INDEL"

gatk --java-options "-Xmx100g" SelectVariants -R $REF_DB -V $PREFIX/variants/raw_variants.vcf -O $PREFIX/variants/filtered_variants.vcf --restrict-alleles-to BIALLELIC --select-type-to-include SNP --select-type-to-include INDEL --selectExpressions "AF > 0.01 && AF < 0.99 && QD > 2.0 && MQ > 40.0 && FS < 60.0 && SOR < 3.0 && ReadPosRankSum > -8.0 && MQRankSum > -12.5 && QUAL > 30.0" >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1
