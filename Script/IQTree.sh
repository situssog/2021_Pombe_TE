#!/bin/bash
#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 144:00:00
#SBATCH -J IQTree


module load bioinfo-tools iqtree/1.6.10-omp-mpi
module load gcc/6.3.0 openmpi/2.1.0

input=$1
output=$2

echo "input file=" $1
echo "output file=" $1
echo "directory=" $PWD

iqtree-omp -s $input -pre $output -b 1000 -nt 16

echo "Done!"

date
echo "--------------------------"
