#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time 24:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH -D ANOTHERPATH
#SBATCH -J trim
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=FAIL
#SBATCH -p batch
#SBATCH --account=USERNAME
#SBATCH -o PATH/THEFILE_trim_%A.log
#SBATCH -e PATH/THEFILE_trim_%A.err

module load anaconda/3_5.0.1_20180125
module load trim_galore/0.6.10
module load cutadapt/3.4

trim_galore --output_dir ../trim/ --paired --fastqc RAWPATHTHE2FILES1 RAWPATHTHE2FILES2

# remove human sequences 
module load bowtie2/2.5.2
module load samtools/1.18
module load bedtools2/2.31.0

bowtie2 -x ../genome_human/GRCh38_noalt_as/GRCh38_noalt_as -1 ../trim/TRIM1 -2 ../trim/TRIM2 -p 2 -S ../trim/THEFILE.sam
rm ../trim/TRIM1
rm ../trim/TRIM2
samtools view -bS ../trim/THEFILE.sam > ../trim/THEFILE.bam -@ 2
rm ../trim/THEFILE.sam
samtools view -b -f 12 -F 256 ../trim/THEFILE.bam > ../trim/THEFILE.filtered.bam -@ 2
rm ../trim/THEFILE.bam
samtools sort -n ../trim/THEFILE.filtered.bam -o ../trim/THEFILE.filtered.sorted.bam -@ 2 -m 20G
rm ../trim/THEFILE.filtered.bam
bedtools bamtofastq -i ../trim/THEFILE.filtered.sorted.bam -fq ../trim/TRIM1_clean.fastq -fq2 ../trim/TRIM2_clean.fastq
rm ../trim/THEFILE.filtered.sorted.bam
cat ../trim/TRIM1_clean.fastq ../trim/TRIM2_clean.fastq > ../trim/TRIMBASE_clean.fastq
gzip -f -1 ../trim/TRIMBASE_clean.fastq
# 
rm ../trim/TRIM1_clean.fastq 
rm ../trim/TRIM2_clean.fastq
