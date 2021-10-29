#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=12:00:00
#SBATCH -J RepeatMarsk_LTR

# conda activate env_others


sample=$1
echo $sample
mkdir repeatMask_$sample"_"FullLtr
cd repeatMask_$sample"_"FullLtr
# identify full-length LTRs
RepeatMasker -pa 1 -a -nolow -norna -dir ./ -lib ../tf_simpleID.fasta ../$sample"_allSeq_EDTA_masked_conSeq_minLen"0.fasta
mv $sample"_"allSeq_EDTA_masked_conSeq_minLen0.fasta.out ../$sample"_"allSeq_EDTA_masked_conSeq_minLen0.fullltr





