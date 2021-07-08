#!/bin/bash

echo "        

ISSRseq -- CreateMatrices
                       
development version 0.4
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
-S <min samples [integer] for variants inclusion in fasta,nexus,phylip matrices>
-D <minimum physical distance [integer] allowed between variants -- set this value to higher than the longest de novo contig to obtain thinned matrices with one variant per locus>
-M <maximum percent missing data [floating point value between 0 & 1] allowed per variant, variants are removed for which missing data is larger than this value>

Dependencies: vcf2phylip

"

exit 1
fi

while getopts "O:T:S:D:M:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        S) MIN_SAMPLES=$OPTARG ;;
        D) THINNING_DISTANCE=$OPTARG ;;
        M) MAX_MISSING_DATA=$OPTARG ;;
       esac
done

if [ ! -d "$PREFIX/matrices" ]; then
    mkdir $PREFIX/matrices
fi

#########################################################################

#create a VCF file filtered to user-specified maximum percentage missing data
vcftools --vcf $PREFIX/variants/filtered_variants.vcf --max-missing $MAX_MISSING_DATA --recode --recode-INFO-all --temp $PREFIX/matrices --out $PREFIX-missing_filtered_variants >>$PREFIX/ISSRseq_CreateMatrices.log 2>&1
mv $PREFIX-missing_filtered_variants.recode.vcf $PREFIX/variants/missing_filtered_variants.vcf

#create a VCF file thinned to user-specified minimum distance, set to longer than maximum contig length for de novo assembled reference
vcftools --vcf $PREFIX/variants/missing_filtered_variants.vcf --thin $THINNING_DISTANCE --recode --recode-INFO-all --temp $PREFIX/matrices --out $PREFIX-thinned_filtered_variants >>$PREFIX/ISSRseq_CreateMatrices.log 2>&1
mv $PREFIX-thinned_filtered_variants.recode.vcf $PREFIX/variants/missing_thinned_filtered_variants.vcf

#VCF to phylip, nexus, and fasta formats using python program available from: https://github.com/edgardomortiz/vcf2phylip
vcf2phylip -i $PREFIX/variants/missing_filtered_variants.vcf -m $MIN_SAMPLES -f -n -b >>$PREFIX/ISSRseq_CreateMatrices.log 2>&1
mv $PREFIX/variants/missing_filtered_variants.min* $PREFIX/matrices

#thinned VCF to phylip, nexus, and fasta formats using python program available from: https://github.com/edgardomortiz/vcf2phylip
vcf2phylip -i $PREFIX/variants/missing_thinned_filtered_variants.vcf -m $MIN_SAMPLES -f -n -b >>$PREFIX/ISSRseq_CreateMatrices.log 2>&1
mv $PREFIX/variants/missing_thinned_filtered_variants.min* $PREFIX/matrices
