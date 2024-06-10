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

trim_galore --output_dir ../trim/ --fastqc RAWPATHTHEFILE

# remove human sequences 
module load bowtie2/2.5.2
module load samtools/1.18
module load bedtools2/2.31.0

bowtie2 -x ../genome_human/GRCh38_noalt_as/GRCh38_noalt_as -U ../trim/TRIM -p 2 -S ../trim/THEFILE.sam
rm ../trim/TRIM
samtools view -bS ../trim/THEFILE.sam > ../trim/THEFILE.bam -@ 2
rm ../trim/THEFILE.sam
samtools view -b -f 4 -F 256 ../trim/THEFILE.bam > ../trim/THEFILE.filtered.bam -@ 2
rm ../trim/THEFILE.bam
samtools sort -n ../trim/THEFILE.filtered.bam -o ../trim/THEFILE.filtered.sorted.bam -@ 2 -m 20G
rm ../trim/THEFILE.filtered.bam
bedtools bamtofastq -i ../trim/THEFILE.filtered.sorted.bam -fq ../trim/TRIM_clean.fastq 
gzip -f -1 ../trim/TRIM_clean.fastq
# 
rm ../trim/THEFILE.filtered.sorted.bam
