#!/bin/bash
#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 5:00:00
#SBATCH -J pbmm2_PB

sample=$3
echo $sample


module load conda 
source conda_init.sh
conda activate /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/all21_samples/DeepVariant/test_conda/t2/pbmm2

# pbmm2 index Pomberef_ed.fasta Pomberef_ed.mmi
pbmm2 align $1 $2 "job_"$3"_"conda_pbmm2_sorted.bam --sort -j 8 -J 8 --preset CCS --rg '@RG\tID:$sample\tSM:mysample_$sample'

pbmm2 align $1 $2 "job_"_$3"_"conda_pbmm2_sorted_sample.bam --sort -j 8 -J 8 --preset CCS --sample $sample

