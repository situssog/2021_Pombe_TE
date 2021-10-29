#!/bin/bash
#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 3:00:00
#SBATCH -J BWA_PB

module load bioinfo-tools bwa/0.7.17
module load picard/2.20.4


sample=$3
ReadGroup="@RG\tID:$sample\tSM:Sample\tPL:PacBio\tLB:lib1\tPU:unit1" 

echo $sample

bwa mem -x pacbio -t 8 -R $ReadGroup $1 $2 | java -jar /sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar SortSam \
  INPUT=/dev/stdin \
  OUTPUT="$sample.R2.bwa.picardSort.bam" \
  SORT_ORDER=coordinate



# module load samtools
# samtools view -b $3"_BWA.sam" > $3"_BWA.bam"


