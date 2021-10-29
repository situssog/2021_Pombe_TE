#!/bin/bash
#SBATCH -A snic2018-3-655
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 140:00:00
#SBATCH -J RaxMl_NG

# ================== Molecular systematics =================
module load bioinfo-tools
module load RAxML-NG

echo "--------------------------"
date

alignment=$1
model_used=$2
CORES=10  # No. of cores used. On Milou there are 16 cores per node. An example, specifying 2 nodes, and thus 32 (2 * 16) cpus, would be -n 32 (above)
#RAXML=raxmlHPC-PTHREADS-AVX   # raxmlHPC
#searches=30
#boots=200

#raxml-ng --parse --msa $alignment --model $model_used --threads $CORES --force 

# if there is a model:
# raxml-ng --all --msa $alignment --model $model_used --bs-trees autoMRE --threads $CORES --force

raxml-ng --all --msa $alignment --bs-trees autoMRE --threads $CORES --force

# TIM3+F+R10
# TN+F+R3

echo "Done!"

date
echo "--------------------------"
