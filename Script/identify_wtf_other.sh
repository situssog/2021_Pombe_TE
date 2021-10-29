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



mkdir repeatMask_$sample"_"wtf
cd repeatMask_$sample"_"wtf
blastn -query ../$sample"_allSeq_EDTA_masked_conSeq_minLen"0.fasta -db ../all_wtf_pom_kam.fasta -evalue 0.00000001 -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send slen qlen evalue bitscore" > blast_wtfInTFseq.$sample.txt
cat blast_wtfInTFseq.$sample.txt | awk '{if($4>700) print $1}' | sort | uniq > wtfSeq_ID.$sample.txt
# makeblastdb -in SPNCRNA.fasta -dbtype nucl
blastn -query ../$sample"_allSeq_EDTA_masked_conSeq_minLen"0.fasta -db ../SPNCRNA.fasta -evalue 0.00000001 -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send slen qlen evalue bitscore" > blast_SPNCRNAInTFseq.$sample.txt
cat blast_SPNCRNAInTFseq.$sample.txt | awk '{if($4>700) print $1}' | sort | uniq > SPNCRNASeq_ID.$sample.txt
cat wtfSeq_ID.$sample.txt SPNCRNASeq_ID.$sample.txt | sort | uniq > ../wtfSeq_SPNCRNASeq_ID.$sample.txt









