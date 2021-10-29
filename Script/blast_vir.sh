#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH -J blast_vir


sample=$1
cd dir_$sample/results_classify

module load bioinfo-tools blast/2.7.1+
tblastx -db /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/sequence_Fungi_plusRetroViruses.fasta -query notKnown.fa -max_hsps 1 -seg no -evalue 0.00001 -num_threads 8 -max_target_seqs 1 -word_size 2 -outfmt 6 -out notKnown.fa.ervwb.ncbi
awk '{print $1"\t""blast""\t""hit""\t"$7"\t"$8"\t"$11"\t"".""\t"".""\t""Target sp|"$2" "$9" "$10}' notKnown.fa.ervwb.ncbi > tmp3
awk '{if($4>$5) print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12; else print $0}' tmp3 > notKnown.fa.ervwb.gff



