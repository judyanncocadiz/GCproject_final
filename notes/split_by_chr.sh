#Splitting bam files by chromosome
bamtools split -in GC115pre.merged.bam -reference


## create index for each file
for file in *
do
echo "indexing "$file
samtools index ${file} ${file}.bai
done