## subsampled bams (75%)
for infile in *bam
do
base=$(basename ${infile} .bam)
echo "subsampling "$base.bam
samtools view -bs 0.75 $base.bam > subsampled_75p_$base.bam
echo "saved as "subsampled_75p_$base.bam
done

mv subsampled_*.bam ~/Documents/GCproject_final/data/bams/subsampled_bams/

## subsampled bams (50%)
for infile in *bam
do
base=$(basename ${infile} .bam)
echo "subsampling "$base.bam
samtools view -bs 0.5 $base.bam > subsampled_50p_$base.bam
echo "saved as "subsampled_50p_$base.bam
done

mv subsampled_*.bam ~/Documents/GCproject_final/data/bams/subsampled_bams/
  
cd ~/Documents/GCproject_final/data/bams/subsampled_bams/
  
## create index for each file
for file in *.bam
do
echo "indexing "$file
samtools index ${file} ${file}.bai
done

# sort all bam files
# cd to bam directory
for infile in *bam
do
base=$(basename ${infile} .bam)
echo "sorting "$base.bam
samtools sort $base.bam > $base.sorted.bam
echo "saved as "$base.sorted.bam
done