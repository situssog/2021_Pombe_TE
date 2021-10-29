#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 6:00:00
#SBATCH -J final_annotation

module load java


sample=$1
ref_seq=$2

cd $ref_seq
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
cp /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/all21_samples/CARP_TF/$ref_seq.fa ./annotationfiles/pombe_repeats.fa
touch ./results_classify/SSR.txt

cp /proj/uppstore2017159/b2014286/private/tools/CARP/carp_for_raylab/code/GenerateAnnotatedLibrary.java ./
sed -i 's/\/home\/a1635743\/RepBase20.04.fasta/\.\//g' ./GenerateAnnotatedLibrary.java

javac GenerateAnnotatedLibrary.java
java GenerateAnnotatedLibrary

# Identify annotated regions in the genome:
echo "Identify annotated regions in the genome"

module load bioinfo-tools RepeatMasker/4.0.8

cat /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/all21_samples/CARP_TF/$ref_seq.fa ./finallibrary/Denovo_TE_Library.fasta > combined_library_Withwtf.fa

makeblastdb -in combined_library_Withwtf.fa -dbtype nucl
cp ../../hg19v37.mfa $sample"_"annotation.fa

sed -i 's/_unitig_/_un_/g' $sample"_"annotation.fa
sed -i 's/|quiver_/_qu_/g' $sample"_"annotation.fa
sed -i 's/_tig000000/_t/g' $sample"_"annotation.fa
sed -i 's/_Segkk/_S/g' $sample"_"annotation.fa
sed -i 's/|arrow_/_a_/g' $sample"_"annotation.fa

RepeatMasker -pa 4 -a -nolow -norna -dir ./ -lib combined_library_Withwtf.fa $sample"_"annotation.fa

sed -i 's/_un_/_unitig_/g' $sample"_"annotation.fa.out
sed -i 's/_qu_/|quiver_/g' $sample"_"annotation.fa.out
sed -i 's/_t/_tig000000/g' $sample"_"annotation.fa.out
sed -i 's/_S/_Segkk/g' $sample"_"annotation.fa.out
sed -i 's/_a_/|arrow_/g' $sample"_"annotation.fa.out



cp ../../hg19v37.mfa $sample"_"general_annotation.fa

sed -i 's/_unitig_/_un_/g' $sample"_"general_annotation.fa
sed -i 's/|quiver_/_qu_/g' $sample"_"general_annotation.fa
sed -i 's/_tig000000/_t/g' $sample"_"general_annotation.fa
sed -i 's/_Segkk/_S/g' $sample"_"general_annotation.fa
sed -i 's/|arrow_/_a_/g' $sample"_"general_annotation.fa

RepeatMasker -s -species ascomycetes -pa 4 $sample"_"general_annotation.fa

sed -i 's/_un_/_unitig_/g' $sample"_"general_annotation.fa.out
sed -i 's/_qu_/|quiver_/g' $sample"_"general_annotation.fa.out
sed -i 's/_t/_tig000000/g' $sample"_"general_annotation.fa.out
sed -i 's/_S/_Segkk/g' $sample"_"general_annotation.fa.out
sed -i 's/_a_/|arrow_/g' $sample"_"general_annotation.fa.out


