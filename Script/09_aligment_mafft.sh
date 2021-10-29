#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:00:00
#SBATCH -J MAFFT

module load bioinfo-tools MAFFT/7.407

sample=$1
#name=$2
dir_run=$PWD

#mkdir $name
#cd $name

# mafft --thread 10 --threadtb  --threadit 0 --reorder --adjustdirection --auto input > output

mafft  --adjustdirection --auto --reorder --thread 8 $dir_run"/"$sample > alig_$sample

