#!/bin/bash

echo "        

ISSRseq -- CreateMatrices
                       
development version 0.2
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
-T <number of parallel threads>
-S <min samples for SNP inclusion in fasta,nexus,phylip matrices>


Dependencies: vcf2phylip

"

exit 1
fi

while getopts "O:T:S:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        S) MIN_SAMPLES=$OPTARG ;;
       esac
done

if [ ! -d "$PREFIX/matrices" ]; then
    mkdir $PREFIX/matrices
fi

#########################################################################

#VCF to phylip, nexus, and fasta formats using python program available from: https://github.com/edgardomortiz/vcf2phylip

vcf2phylip -i $PREFIX/variants/filtered_SNPs.vcf -m $MIN_SAMPLES -f -n -b >>$PREFIX/ISSRseq_CreateMatrices.log 2>&1
mv $PREFIX/variants/filtered_SNPs.min* $PREFIX/matrices
