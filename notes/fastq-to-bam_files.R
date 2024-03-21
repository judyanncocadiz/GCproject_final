## Do everything in bash

##Check the quality of fastq files##
## pwd
## cd data
## cd GC110post
## for file in * 
## >do
## >echo $file
## >fastqc $file -o ~/Documents/GC_DELFIproject/output/fastqc/GC121post
## >done


## pwd
## cd ..
## cd ..
## cd output
## cd fastqc
## pwd
## ls
## multiqc . -o ~/Documents/GC_DELFIproject/output/multiqc/GC121post

##Trim adapters using TRIMMOMATIC ##
## cd ..
## cd ..
## cd data
## cd GC110post

## for infile in ./*_R1_001.fastq.gz
## >do
## >base=$(basename ${infile} _R1_001.fastq.gz)
## >java -jar /Volumes/userdata/student_users/judyanncocadiz/Documents/R/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${infile} ${base}_R2_001.fastq.gz \
## >${base}_R1.trim.fastq.gz ${base}_R1un.trim.fastq.gz \
## >${base}_R2.trim.fastq.gz ${base}_R2un.trim.fastq.gz \
## >ILLUMINACLIP:TruSeq3-PE.fa:1:10:5:9:true MINLEN:20
## >done

## mkdir trimmed_fastq
## mv *.trim* trimmed_fastq
## cd trimmed_fastq
## pwd
## ls

## Run fastqc again 
## for file in * 
## >do
## >echo $file
## >fastqc $file -o ~/Documents/GC_DELFIproject/output/trimmed-fastqc/GC121post
## >done


## Run multiqc
## pwd
## cd ..
## cd ..
## cd output
## cd trimmed-fastqc
## pwd
## ls
## multiqc . -o ~/Documents/GC_DELFIproject/output/trimmed-multiqc/GC121post

## Align with bwa
## mkdir bwa_aligned
## cd trimmed_fastq
## for infile in *1.trim.*
## > do
## > base=$(basename ${infile} 1.trim.fastq.gz)
## > bwa mem -t 5 /Volumes/userdata/student_users/judyanncocadiz/Documents/R/align/bwa-0.7.17/hg19index \
## > ${infile} ${base}2.trim.fastq.gz >/Volumes/userdata/student_users/judyanncocadiz/Documents/GC_DELFIproject/data/GC121post/bwa_aligned/${base}.aligned.sam
## > done


## Convert sam to bam ##
## cd ..
## cd bwa_aligned

##OPTION 1
## for infile in *.sam
## > do
## > base=$(basename ${infile} .aligned.sam)
## > samtools view -bo ${base}.aligned.bam ${infile}
## > done

## or 

##OPTION 2 (from the bash couse notes; both work fine)
## for filename in *.sam
## > do
## > base=$(basename ${filename} .sam)
## > samtools view -S -b ${filename} -o ${base}.bam
## > done


## Sort bam files ##
## mkdir aligned_bam
## mv *.bam aligned_bam
## cd aligned_bam
## for infile in *.aligned.bam
## > do
## > base=$(basename ${infile} .aligned.bam)
## > samtools sort -o ${base}.sorted.bam ${infile}
## > done

## cd ..
## mkdir sorted_bam
## cd sorted_bam


## merge files from the same library ##
## samtools merge -O BAM ./GC121post.merged.bam *


## Flagstat
## samtools flagstat GC121post.merged.bam



