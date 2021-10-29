#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 6:00:00
#SBATCH -J final_annotation

module load java

sample=$1
cd dir_$sample/results_classify

# get protain annotation:
echo "get protains"
cp /proj/uppstore2017159/b2014286/private/tools/CARP/carp_for_raylab/code/GetProteins.java ./
sed -i 's/ProteinReport/./g' GetProteins.java
javac GetProteins.java
java GetProteins

# produce compiled annotation:
echo "Identify annotated regions in the genome"

mkdir final
cd final
mkdir results_classify
mkdir finallibrary
mkdir annotationfiles
cp ../../ConsensusSequences.fa ../../ConsensusSequences.fa.map ../notKnown.fa.ervwb.gff ../protein.txt ../known.txt ../notKnown.fa.tewb.gff ./results_classify/ 
cp ../../ConsensusSequences.fa ../../ConsensusSequences.fa.map ./
cp /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/GB_TE.fa ./annotationfiles/GB_TE.fa
cp /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/sequence_Fungi_plusRetroViruses.fa ./annotationfiles/all_retrovirus.fasta
cp /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/pombe_repeats_format.fa ./annotationfiles/pombe_repeats.fa
touch ./results_classify/SSR.txt


cp /proj/uppstore2017159/b2014286/private/tools/CARP/carp_for_raylab/code/GenerateAnnotatedLibrary.java ./
sed -i 's/\/home\/a1635743\/RepBase20.04.fasta/\.\//g' ./GenerateAnnotatedLibrary.java

javac GenerateAnnotatedLibrary.java
java GenerateAnnotatedLibrary

# Identify annotated regions in the genome:
echo "Identify annotated regions in the genome"

module load bioinfo-tools RepeatMasker/4.0.8

cat /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/pombe_repeats.fa ./finallibrary/Denovo_TE_Library.fasta /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/wtf_ref.fasta /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/all_wtf_CBS5557.fa > combined_library_Withwtf.fa
makeblastdb -in combined_library_Withwtf.fa -dbtype nucl
cp ../../hg19v37.mfa $sample"_"annotation.fa
RepeatMasker -pa 4 -a -nolow -norna -dir ./ -lib combined_library_Withwtf.fa $sample"_"annotation.fa


cp ../../hg19v37.mfa $sample"_"general_annotation.mfa
RepeatMasker -s -species ascomycetes -pa 4 $sample"_"general_annotation.mfa

