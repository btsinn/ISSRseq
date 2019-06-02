#!/bin/bash

echo "        

ISSRseq -- AssembleReference
                       
development version 0.1
use help for usage 
    
"

if [[ $1 = help ]]
  then
  
echo "This script uses processes input reads and uses ABYSS to generate a reference de novo from a designated ISSRseq sample.

Usage is as follows:

REQUIRED:

-O <output_prefix> 
-I <read directory>
-S <sample list> 
-R <read prefix for sample to be used for reference assembly>
-T <number of parallel threads>
-M <minimum post-trim read length>
-H <N bases to hard trim at each end of reads>
-P <fasta file of ISSR motifs used>
-K <length of shortest primer used>


Dependencies: ABYSS-PE, bbduk

"

exit 1
fi

while getopts "O:I:S:R:T:M:H:P:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        I) READ_DIR=$OPTARG ;;
        S) SAMPLE_LIST=$OPTARG ;;
        R) REF_SAMPLE=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        M) MIN_LENGTH=$OPTARG ;;
        H) HARD_TRIM=$OPTARG ;;
        P) ISSR_MOTIF=$OPTARG ;;
       esac
done      

#create output directory naming variables
TIMESTAMP=$(date '+%Y_%m_%d_%H_%M')
OUTPUT_DIR="$PREFIX"_AssembleReference_"$TIMESTAMP"
SRC=$(pwd)

mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/trimmed_reads
mkdir $OUTPUT_DIR/reference

cp $SAMPLE_LIST $OUTPUT_DIR

#########################################################################


#process reads
while read -r sample
do

    bbduk in=$READ_DIR/${sample}_R1.fastq in2=$READ_DIR/${sample}_R2.fastq forcetrimleft=$HARD_TRIM forcetrimright2=$HARD_TRIM qtrim=t trimq=20 minlength=$MIN_LENGTH threads=$THREADS out=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R1.fastq out2=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_read_trimming.log 2>&1
    
echo ""${sample}" reads processed"

done < $SAMPLE_LIST

#kmer trim reference sample reads and assemble reference amplicons

echo "

starting de novo reference assembly using $REF_SAMPLE
"

wait

bbduk in=$OUTPUT_DIR/trimmed_reads/''$REF_SAMPLE''_trimmed_R1.fastq in2=$OUTPUT_DIR/trimmed_reads/''$REF_SAMPLE''_trimmed_R2.fastq minlength=$MIN_LENGTH restrictleft=30 k=18 ktrim=l mink=10 ref=$ISSR_MOTIF threads=$THREADS mingc=0.1 maxgc=0.9 out=$OUTPUT_DIR/reference/Ltrimmed_R1.fastq out2=$OUTPUT_DIR/reference/Ltrimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1

bbduk in=$OUTPUT_DIR/reference/Ltrimmed_R1.fastq in2=$OUTPUT_DIR/reference/Ltrimmed_R2.fastq minlength=$MIN_LENGTH restrictright=30 k=18 ktrim=r mink=10 ref=$ISSR_MOTIF threads=$THREADS mingc=0.1 maxgc=0.9 out=$OUTPUT_DIR/reference/Ktrimmed_R1.fastq out2=$OUTPUT_DIR/reference/Ktrimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1

wait 

abyss-pe -C $OUTPUT_DIR/reference name=reference k=41 np=$THREADS lib='paired' paired='Ktrimmed_R1.fastq Ktrimmed_R2.fastq' >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1

wait

bbduk in=$OUTPUT_DIR/reference/reference-contigs.fa k=15 minlength=250 ktrim=r mink=10 ref=$ISSR_MOTIF threads=$THREADS mingc=0.1 maxgc=0.9 out=$OUTPUT_DIR/reference/reference-contigs_R250bpMIN.fa >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1

bbduk in=$OUTPUT_DIR/reference/reference-contigs_R250bpMIN.fa k=15 minlength=250 ktrim=l mink=10 ref=$ISSR_MOTIF threads=$THREADS mingc=0.1 maxgc=0.9 trd=t out=$OUTPUT_DIR/reference/reference_assembly.fa >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1
