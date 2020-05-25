#!/bin/bash

echo "        

ISSRseq -- SPadesReference
                       
development version 0.1
use help for usage 
    
"

if [[ $1 = help ]]
  then
  
echo "This script processes input reads and uses SPades to generate a reference de novo from a designated ISSRseq sample, then remove loci assembled from common contaminants.

Usage is as follows:

REQUIRED:

-O <output prefix> 
-I <read directory>
-S <sample list> 
-R <read prefix for sample to be used for reference assembly>
-T <number of parallel threads>
-M <minimum post-trim read length>
-H <N bases to hard trim at each end of reads>
-P <fasta file of ISSR motifs used>
-K <single kmer, or range seperated by commas, for SPades assembly ie. 33,55,77,99,127>
-L <minimum contig length>
-N <negative reference in fasta format to filter reads against, ex. sequenced plastome or specific contaminant>


Dependencies: SPades, bbduk, blastn

added common plant contaminant locus filter
contaminant loci are saved to a seperate fasta (reference/contaminant_loci.fa)
introduced entropy filtering -- loci must have entropy of 0.85 within a sliding window of 25 bp, kmer of 5
introduced GC filtering -- loci must contain between 35% and 65% GC content
enabled tbo and tpe trimming during BBDUK read trimming
removed secondary trimming of reference sample reads prior to assembly 
reduced BBDUK mink setting to 8 and trimming kmer to 18
enabled a GC content filtering of reads below 0.1 and above 0.9 
decreased the contaminant filter evalue to 0.00001
added the N flag to negative filter each read pool against a reference sequence

"

exit 1
fi

while getopts "O:I:S:R:T:M:H:P:L:K:N:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        I) READ_DIR=$OPTARG ;;
        S) SAMPLE_LIST=$OPTARG ;;
        R) REF_SAMPLE=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        M) MIN_LENGTH=$OPTARG ;;
        H) HARD_TRIM=$OPTARG ;;
        P) ISSR_MOTIF=$OPTARG ;;
        K) SPADES_K=$OPTARG ;;
        L) MIN_CONTIG=$OPTARG ;;
        N) NEG_REF=$OPTARG ;;
       esac
done      

#create output directory naming variables
TIMESTAMP=$(date '+%Y_%m_%d_%H_%M')
OUTPUT_DIR="$PREFIX"_"$TIMESTAMP"
SRC=$(pwd)
BLASTDB=/usr/local/src/ISSRseq/contam_filters/default_contam_filter_BLASTDB/default_contam_filter.fasta

mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/trimmed_reads
mkdir $OUTPUT_DIR/reference

cp $SAMPLE_LIST $OUTPUT_DIR/samples.txt

#########################################################################


#process reads
while read -r sample
do

    bbduk in=$READ_DIR/${sample}_R1.fastq in2=$READ_DIR/${sample}_R2.fastq mingc=0.1 maxgc=0.9 forcetrimleft=$HARD_TRIM forcetrimright2=$HARD_TRIM qtrim=t trimq=10 k=18 tbo=t tpe=t ktrim=l mink=8 mingc=0.1 maxgc=0.9 ref=$ISSR_MOTIF minlength=$MIN_LENGTH threads=$THREADS out=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R1.fastq out2=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_read_trimming.log 2>&1
    
    bbduk in=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R1.fastq in2=$OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R2.fastq minlength=$MIN_LENGTH qtrim=t trimq=10 k=18 tbo=t tpe=t  ktrim=r mink=8 mingc=0.1 maxgc=0.9 ref=$ISSR_MOTIF threads=$THREADS mingc=0.1 maxgc=0.9 out=$OUTPUT_DIR/trimmed_reads/${sample}_LRtrimmed_R1.fastq out2=$OUTPUT_DIR/trimmed_reads/${sample}_LRtrimmed_R2.fastq >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1
    
    bbmap threads=$THREADS nodisk=f killbadpairs=t slow=t ref=$NEG_REF in=$OUTPUT_DIR/trimmed_reads/${sample}_LRtrimmed_R1.fastq in2=$OUTPUT_DIR/trimmed_reads/${sample}_LRtrimmed_R2.fastq outu=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R1.fastq outu2=$OUTPUT_DIR/trimmed_reads/${sample}_trimmed_R2.fastq statsfile=stderr >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1

    rm $OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R1.fastq $OUTPUT_DIR/trimmed_reads/${sample}_Ltrimmed_R2.fastq $OUTPUT_DIR/trimmed_reads/${sample}_LRtrimmed_R1.fastq $OUTPUT_DIR/trimmed_reads/${sample}_LRtrimmed_R2.fastq
    
echo ""${sample}" reads processed"

done < $SAMPLE_LIST

#kmer trim reference sample reads and assemble reference amplicons

echo "

starting de novo reference assembly using $REF_SAMPLE
"

cp $OUTPUT_DIR/trimmed_reads/""$REF_SAMPLE""_trimmed_R1.fastq $OUTPUT_DIR/reference/trimmed_R1.fastq

cp $OUTPUT_DIR/trimmed_reads/""$REF_SAMPLE""_trimmed_R2.fastq $OUTPUT_DIR/reference/trimmed_R2.fastq

spades.py --isolate -1 $OUTPUT_DIR/reference/trimmed_R1.fastq -2 $OUTPUT_DIR/reference/trimmed_R2.fastq -t $THREADS -k $SPADES_K -o $OUTPUT_DIR/reference >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1

bbduk in=$OUTPUT_DIR/reference/contigs.fasta k=18 minlength=$MIN_CONTIG entropy=0.85 entropywindow=25 entropyk=5 mingc=0.35 maxgc=0.65 ktrim=r mink=8 ref=$ISSR_MOTIF threads=$THREADS out=$OUTPUT_DIR/reference/reference-contigs_R""$MIN_CONTIG""bpMIN.fa >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1
 
bbduk in=$OUTPUT_DIR/reference/reference-contigs_R""$MIN_CONTIG""bpMIN.fa k=18 minlength=$MIN_CONTIG entropy=0.85 entropywindow=25 entropyk=5 mingc=0.35 maxgc=0.65 ktrim=l mink=8 ref=$ISSR_MOTIF threads=$THREADS trd=t out=$OUTPUT_DIR/reference/reference_assembly.fa >>$OUTPUT_DIR/ISSRseq_reference_assembly.log 2>&1

grep "^>" $OUTPUT_DIR/reference/reference_assembly.fa | sed 's/>//' > $OUTPUT_DIR/reference/reference_loci_list.txt

makeblastdb -in $OUTPUT_DIR/reference/reference_assembly.fa -input_type fasta -dbtype nucl -parse_seqids -out $OUTPUT_DIR/reference/reference_assemblyDB

#filter contaminants from assembled reference

blastn -db $BLASTDB -query $OUTPUT_DIR/reference/reference_assembly.fa -num_threads $THREADS -word_size 9 -out $OUTPUT_DIR/reference/assembly_contam_blasthits.out -evalue 0.00001 -outfmt "6 qseqid" -max_target_seqs 1

#need to be able to continue if no contaminants are identified ...
sort -u -n -o $OUTPUT_DIR/reference/contam_loci_list.txt $OUTPUT_DIR/reference/assembly_contam_blasthits.out

grep -v -f $OUTPUT_DIR/reference/contam_loci_list.txt $OUTPUT_DIR/reference/reference_loci_list.txt > $OUTPUT_DIR/reference/filtered_loci_list.txt

blastdbcmd -db $OUTPUT_DIR/reference/reference_assemblyDB -dbtype nucl -entry_batch $OUTPUT_DIR/reference/contam_loci_list.txt -outfmt %f -out $OUTPUT_DIR/reference/contaminant_loci.fa

blastdbcmd -db $OUTPUT_DIR/reference/reference_assemblyDB -dbtype nucl -entry_batch $OUTPUT_DIR/reference/filtered_loci_list.txt -outfmt %f -out $OUTPUT_DIR/reference/final_reference_assembly.fa
