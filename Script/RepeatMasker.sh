#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 4:00:00
#SBATCH -J RepeatMasker



module load bioinfo-tools RepeatMasker/4.0.8
module load blast/2.9.0+

# to create library:
#makeblastdb -in pombe_TF_refrepeats.fa -dbtype nucl

RepeatMasker -pa 4 -a -nolow -norna -no_is  -dir ./ -lib $1 $2



