#!/bin/bash

echo "        

ISSRseq -- CreateBAMs
                       
development version 0.3
use help for usage 
    
"

if [[ $1 = help ]]
  then
  
echo "This script takes the output of AssembleReference and creates BAM files for AnalyzeBAMs.

Usage is as follows:

REQUIRED:

-O <AssembleReference output prefix> 
-T <number of parallel threads>
-M <minimum percent identity for read mapping[0-1]>
-I <maximum indel length to allow during mapping>


Dependencies: samtools, picard, bbmap

"

exit 1
fi

while getopts "O:T:M:I:" opt; do

      case $opt in 
        O) PREFIX=$OPTARG ;;
        T) THREADS=$OPTARG ;;
        M) MIN_MAPPING_ID=$OPTARG ;;
        I) MAX_INDEL=$OPTARG ;;
       esac
done

#########################################################################

REF_DB=$PREFIX/reference/final_reference_assembly.fa
SAMPLE_LIST=$PREFIX/samples.txt

mkdir $PREFIX/bams


#create reference dictionary and index

samtools faidx $REF_DB >>$PREFIX/ISSRseq_CreateBAMs.log 2>&1

java -jar /usr/local/src/picard/build/libs/picard.jar CreateSequenceDictionary R=$REF_DB O=$REF_DB.dict >>$PREFIX/ISSRseq_CreateBAMs.log 2>&1

rename 's/.fa.dict/.dict/' $REF_DB.dict

#create use bbmap to create BAM files, name sample groups with picard, and index the modified bams with samtools 

while read -r sample
do
    
    bbmap threads=$THREADS minid=$MIN_MAPPING_ID nodisk=t maxindel=$MAX_INDEL strictmaxindel=t killbadpairs=t ref=$REF_DB in=$PREFIX/trimmed_reads/${sample}_trimmed_R1.fastq in2=$PREFIX/trimmed_reads/${sample}_trimmed_R2.fastq out=$PREFIX/bams/${sample}.bam bamscript=sort.sh >>$PREFIX/ISSRseq_CreateBAMs.log 2>&1
       
     sh sort.sh >>$PREFIX/ISSRseq_CreateBAMs.log 2>&1
        
    java -XX:ParallelGCThreads=$THREADS -jar /usr/local/src/picard/build/libs/picard.jar AddOrReplaceReadGroups I=$PREFIX/bams/${sample}_sorted.bam O=$PREFIX/bams/${sample}_sorted_RG.bam RGID=1 RGLB=${sample} RGPL=illumina RGPU=combined_runs RGSM=${sample} >>$PREFIX/ISSRseq_CreateBAMs.log 2>&1
    
    samtools index $PREFIX/bams/${sample}_sorted_RG.bam >>$PREFIX/ISSRseq_CreateBAMs.log 2>&1

    rm sort.sh $PREFIX/bams/${sample}.bam $PREFIX/bams/${sample}_sorted.bam $PREFIX/bams/${sample}_sorted.bam.bai

echo ""${sample}" complete"
    
done < $SAMPLE_LIST

