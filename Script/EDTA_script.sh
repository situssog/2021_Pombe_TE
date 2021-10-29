#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=4:00:00
#SBATCH -J EDTA

# module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux* ; module load user_spack


sample=$1
perl /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/EDTA/EDTA.pl --genome $sample.fasta --cds cds.fa --curatedlib ./pombe_repeats_withTF1_2_only.fa --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10 --force 1



