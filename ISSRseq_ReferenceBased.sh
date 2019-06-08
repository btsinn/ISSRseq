#!/bin/bash

echo "        

ISSRseq -- ReferenceBased
                       
development version 0.2
use help for usage 
    
"

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
-H <N bases to hard trim at each end of reads>


Dependencies: ABYSS-PE, bbduk

"

exit 1
fi

while getopts "O:I:S:R:T:M:H:P:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        I) READ_DIR=$OPTARG ;;
        S) SAMPLE_LIST=$OPTARG ;;
        R) REF_GENOME=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        M) MIN_LENGTH=$OPTARG ;;
        H) HARD_TRIM=$OPTARG ;;
       esac
done      

#create output directory naming variables
TIMESTAMP=$(date '+%Y_%m_%d_%H_%M')
OUTPUT_DIR="$PREFIX"_"$TIMESTAMP"
SRC=$(pwd)

mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/trimmed_reads
mkdir $OUTPUT_DIR/reference

cp $SAMPLE_LIST $OUTPUT_DIR/samples.txt

#########################################################################


#process reads
while read -r sample
do

    bbduk in=$READ_DIR/${sample}_R1.fastq in2=$READ_DIR/${sample}_R2.fastq mingc=0.1 maxgc=0.9 forcetrimleft=$HARD_TRIM forcetrimright2=$HARD_TRIM qtrim=t trimq=20 k=18 ktrim=l mink=10 ref=$ISSR_MOTIF minlength=$MIN_LENGTH threads=$THREADS out=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R1.fastq out2=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_read_trimming.log 2>&1
    
    bbduk in=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R1.fastq in2=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R2.fastq minlength=$MIN_LENGTH qtrim=t trimq=20 k=18 ktrim=r mink=10 ref=$ISSR_MOTIF threads=$THREADS mingc=0.1 maxgc=0.9 out=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R1.fastq out2=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1
    
    rm $OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R1.fastq $OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R2.fastq
    
echo ""${sample}" reads processed"

done < $SAMPLE_LIST

#setup necessary reference directory structure

cp $REF_GENOME $OUTPUT_DIR/reference/reference_assembly.fa
