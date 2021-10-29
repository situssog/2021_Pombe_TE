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
reference_seq=tf_simpleID.fasta


echo $sample
mkdir repeatMask_$sample"_"FullLtr
cd repeatMask_$sample"_"FullLtr
# identify full-length LTRs and solo LTRs per sequence:
RepeatMasker -pa 1 -a -nolow -norna -dir ./ -lib ../$reference_seq ../$sample"_"tf_masked_conSeq_minLen100.fasta
mv $sample"_"tf_masked_conSeq_minLen100.fasta.out ../$sample"_"tf_masked_conSeq_minLen100.fullltr



