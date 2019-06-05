#!/bin/bash

echo "        

ISSRseq -- AnalyzeBAMs
                       
development version 0.1
use help for usage 
    
"

if [[ $1 = help ]]
  then
  
echo "This script takes the output of CreateBAMs and calls and analyzes SNP variants.

Usage is as follows:

REQUIRED:

-O <AssembleReference output prefix -- DO NOT end with /> 
-T <number of parallel threads>
-M <max percent missing for fastStructure>
-K <max K for fastStructure>


Dependencies: vcftools ,fastStructure, plink

"

exit 1
fi

while getopts "O:T:M:K:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        M) MAX_MISSING_PERC=$OPTARG ;;
        K) FAST_K=$OPTARG ;;
       esac
done

#########################################################################

REF_DB=$PREFIX/reference/reference_assembly.fa
SAMPLE_LIST=$PREFIX/samples.txt
FAST_STRUC_DIR=/usr/local/src/fastStructure

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I/usr/local/include"
export LDFLAGS="-L/usr/local/lib"

#########################################################################


#conduct fastStructure analysis
vcftools --vcf $PREFIX/matrices/filtered_SNPs.recode.vcf --out $PREFIX/faststructure/faststructure_input --plink >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

plink --file $PREFIX/faststructure/faststructure_input --geno $MAX_MISSING_PERC --make-bed --out $PREFIX/faststructure/faststructure_input_plink --noweb >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

for (( i = 1; i <= $FAST_K; i++ )) 
  do

python $FAST_STRUC_DIR/structure.py -K $i --input=$PREFIX/faststructure/faststructure_input_plink --output=$PREFIX/faststructure/faststructure_output >>$PREFIX/ISSRseq_AnalyzeBAMs.log 2>&1

done

python $FAST_STRUC_DIR/chooseK.py --input=$PREFIX/fastStructure/faststructure_output > $PREFIX/faststructure/bestK_result.txt

#execute R script

