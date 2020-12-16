#!/bin/bash

echo "        

ISSRseq -- ReferenceBased
                       
development version 0.8
use help for usage 
    
"

set -o errexit
set -o pipefail
set -o nounset

if [[ $1 = help ]]
  then
  
echo "This script uses processes input reads and prepares the pipeline for use with a pre-existing reference.

Usage is as follows:

REQUIRED:

-O <output_prefix> 
-I <read directory>
-S <sample list> 
-R <path to reference assembly fasta file>
-T <number of parallel threads>
-M <minimum post-trim read length>
-P <fasta of ISSR primers used>
-H <N bases to hard trim at each end of reads>
-X <bbduk trimming kmer, equal to or longer than shortest primer used>


Dependencies: bbmap, bbduk

"

exit 1
fi

while getopts "O:I:S:T:M:H:P:R:X:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        I) READ_DIR=$OPTARG ;;
        S) SAMPLE_LIST=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        M) MIN_LENGTH=$OPTARG ;;
        H) HARD_TRIM=$OPTARG ;;
        P) ISSR_MOTIF=$OPTARG ;;
        R) REF_GENOME=$OPTARG;;
        X) TRIM_K=$OPTARG ;;
       esac
done

#setup necessary reference directory structure
TIMESTAMP=$(date '+%Y_%m_%d_%H_%M')
OUTPUT_DIR="$PREFIX"_"$TIMESTAMP"
SRC=$(pwd)

mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/trimmed_reads
mkdir $OUTPUT_DIR/reference
cp $REF_GENOME $OUTPUT_DIR/reference/final_reference_assembly.fa
cp $SAMPLE_LIST $OUTPUT_DIR/samples.txt

#########################################################################


#process reads
while read -r sample
do

    bbduk in=$READ_DIR/${sample}_R1.fastq in2=$READ_DIR/${sample}_R2.fastq mingc=0.1 maxgc=0.9 forcetrimleft=$HARD_TRIM forcetrimright2=$HARD_TRIM qtrim=rl trimq=10 k=$TRIM_K tbo=t tpe=t ktrim=r ktrim=l mink=8 ref=$ISSR_MOTIF minlength=$MIN_LENGTH threads=$THREADS out=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R1.fastq out2=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_read_trimming.log 2>&1
   
echo ""${sample}" reads processed"

done < $SAMPLE_LIST