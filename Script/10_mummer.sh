#!/bin/bash
#SBATCH -A b2014286
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -J mummer

# this script takes the reference seq and aligns it to a pacbio assembly. 

module load bioinfo-tools MUMmer/3.23
nucmer --maxmatch --nosimplify --prefix=$3 $1 $2 



