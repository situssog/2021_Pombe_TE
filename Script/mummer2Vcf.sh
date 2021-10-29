#!/bin/bash
#SBATCH -A b2014286
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -J snps_fromMUMmer

# this script takes the reference seq and aligns it to a pacbio assembly. 

module load bioinfo-tools BioPerl/1.7.2_Perl5.26.2

perl /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/mummer2Vcf.pl -f $1 $2 > $3 



