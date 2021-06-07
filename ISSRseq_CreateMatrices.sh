#!/bin/bash

echo "        

ISSRseq -- CreateMatrices
                       
development version 0.3
use help for usage 
    
"

set -o errexit
set -o pipefail
set -o nounset

if [[ $1 = help ]]
  then
  
echo "This script takes the output of AnalyzeBAMs as input, filters scored variants, and outputs matrices for downstream analyses.

Usage is as follows:

REQUIRED:

-O <AssembleReference output prefix -- DO NOT end with /> 
-T <number of parallel threads [integer]>
-S <min samples [integer] for SNP inclusion in fasta,nexus,phylip matrices>
-D <minimum distance [integer] allowed between SNPs -- set this value to higher than the longest de novo contig to obtain thinned matrices with one SNP per locus>


Dependencies: vcf2phylip

"

exit 1
fi

while getopts "O:T:S:D:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        S) MIN_SAMPLES=$OPTARG ;;
		    D) THINNING_DISTANCE=$OPTARG ;;
       esac
done

if [ ! -d "$PREFIX/matrices" ]; then
    mkdir $PREFIX/matrices
fi

#########################################################################

#create a VCF file thinned to user-specified minimum distance, set to longer than maximum contig length for de novo assembled reference
vcftools --vcf $PREFIX/variants/filtered_SNPs.vcf --thin $THINNING_DISTANCE --recode --recode-INFO-all --temp $PREFIX/matrices --out $PREFIX-thinned_filtered_SNPs
mv $PREFIX-thinned_filtered_SNPs.recode.vcf $PREFIX/variants/thinned_filtered_SNPs.vcf

#VCF to phylip, nexus, and fasta formats using python program available from: https://github.com/edgardomortiz/vcf2phylip
vcf2phylip -i $PREFIX/variants/filtered_SNPs.vcf -m $MIN_SAMPLES -f -n -b >>$PREFIX/ISSRseq_CreateMatrices.log 2>&1
mv $PREFIX/variants/filtered_SNPs.min* $PREFIX/matrices

#thinned VCF to phylip, nexus, and fasta formats using python program available from: https://github.com/edgardomortiz/vcf2phylip
vcf2phylip -i $PREFIX/variants/thinned_filtered_SNPs.vcf -m $MIN_SAMPLES -f -n -b >>$PREFIX/ISSRseq_CreateMatrices.log 2>&1
mv $PREFIX/variants/thinned_filtered_SNPs.min* $PREFIX/matrices
rm $PREFIX-thinned_filtered_SNPs.log

