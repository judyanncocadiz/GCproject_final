#!/bin/bash

## NOTE: This script has several requirements to be able to run...
# Need to be in a directory containing four folders
  # a directory containing the bam files you want to split.
    # bam file names must match the pattern [A-Z]-[A-Z][...].merged.bam
    # eg. A-M011-H-T2.merged.bam
  # an empty directory named sorted
  # an empty directory named short
  # an empty directory named long

# This results in files being split into two bams
# containing short (100-150) and long (151-220) fragments



# sort all bam files
for infile in ./[A-Z]-[A-Z]*bam
  do
  base=$(basename ${infile} .merged.bam)
  echo "sorting "$base.merged.bam
  samtools sort $base.merged.bam > $base.sorted.bam
  echo "saved as "$base.sorted.bam
  done

mv *sorted.bam ../sorted_bams
cd ../sorted_bams

## split to short


for file in *
  do
   base=`basename $file .sorted.bam`
   echo "getting short reads for "$base
   samtools view -h $file |\
   awk 'substr($0,1,1)=="@" || ($9>= 100 && $9<=149) || ($9<=-100 && $9>=-149)' | \
   samtools view -b > ./${base}.short.bam
  done
mv *short.bam ../short

## split to long

for file in *
  do
   base=`basename $file .sorted.bam`
   echo "getting long reads for "$base
   samtools view -h $file |\
   awk 'substr($0,1,1)=="@" || ($9>= 150 && $9<=220) || ($9<=-150 && $9>=-220)' | \
   samtools view -b > ./${base}.long.bam
  done
mv *long.bam ../long


## create index for each file

cd ../short

for file in *
do
echo "indexing "$file
samtools index ${file} ${file}.bai
done

cd ../long

for file in *
do
echo "indexing "$file
samtools index ${file} ${file}.bai
done


