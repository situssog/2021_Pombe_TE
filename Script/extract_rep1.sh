#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 6:00:00
#SBATCH -J extract_repeats

module load bioinfo-tools
module load blast/2.2.26
#module load blast/2.9.0+
module load go/1.11.5
module load git/2.16.1
export PATH=$PATH:$(go env GOPATH)/bin
module load  BioPerl/1.7.2_Perl5.24.1
#export PATH=$PATH:/sw/apps/bioinfo/blast/2.9.0+/rackham/bin
export PATH=$PATH:/sw/apps/bioinfo/blast/2.2.26/rackham/bin
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/censor-4.2.29/bin
module load muscle/3.8.31


sample=$1
mkdir dir_$sample
cd dir_$sample
cp ../../$sample"_ed.fa" ./

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $sample"_ed.fa"

matrix -threads=16 -krishnaflags="-tmp=./ -threads=15 -log -filtid=0.80 -filtlen=100" myseq*.fa

find ./ -maxdepth 1 -name 'myseq[!.]*.gff' -print0 | xargs -r0 cat > hg_krishna.gff

igor -in hg_krishna.gff -out hg94_krishna.json

gffer < hg94_krishna.json > hg94_krishna.igor.gff

cat myseq*.fa > hg19v37.mfa

seqer -aligner=muscle -dir=consensus -fasta=true -maxFam=100 -subsample=true -minLen=0.95 -threads=15 -ref=hg19v37.mfa hg94_krishna.igor.gff

module load bioinfo-tools blast/2.2.26
#blast/2.9.0+
module load java

find ./consensus -maxdepth 1 -name '[!.]*.fq' -print0 | xargs -r0 cat > ConsensusSequences.fa

cp /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/fungi.fa /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/pombe_repeats_withTF1.fa ./

censor.ncbi -debug -lib ./fungi.fa -lib ./pombe_repeats_withTF1.fa ConsensusSequences.fa

mkdir results_classify

cp /proj/uppstore2017159/b2014286/private/tools/CARP/carp_for_raylab/code/ClassifyConsensusSequences.java ./

sed -i 's/annotationfiles\/Vertebrate_use.fa/fungi.fa/g' ClassifyConsensusSequences.java

sed -i 's/annotationfiles\/our_known_reps_20130520.fasta/pombe_repeats_withTF1.fa/g' ClassifyConsensusSequences.java

javac ClassifyConsensusSequences.java

java ClassifyConsensusSequences

