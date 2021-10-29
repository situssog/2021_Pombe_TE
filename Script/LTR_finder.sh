#!/bin/bash

#SBATCH -A snic2018-3-658
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 6:00:00
#SBATCH -J LTR_FINDER_parallel



sample=$1
module load bioinfo-tools LTR_Finder/1.0.7

# $LTR_FINDER_GTRNADB/GtRNAdb-all-tRNAs.fa

perl LTR_FINDER_parallel/LTR_FINDER_parallel -seq $sample.fasta -threads 3 -s GtRNAdb-all-tRNAs.fa  > $sample.out



