#!/bin/bash

#SBATCH -A snic2018-8-12 
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 48:00:00
#SBATCH -J PCA_byWindow_popStat

##snic2018-8-12

module load R/3.4.3 MariaDB/10.2.11
R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4

windowSize=$1
overlap=$2
chromosome=$3
vcf_file=$4

mkdir PCA_byWindow_$1"_"$2

Rscript --vanilla /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/PCA_plots_byWindows_withPlots.R $chromosome $windowSize $overlap $vcf_file

