#!/bin/bash
#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 24:00:00
#SBATCH -J canu_1_8

# this script will run canu, using PacBio reads and assambly a de-novo genome
# it first create a forder to put the output
# example: /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/01_assembly_canu_general_2.sh ebc7_13_8 13.8m /proj/b2014286/private/raw-data/pb_237/filtered_subreads/pb_237_filtered_subreads.fastq.gz 0.013

sample=$1
output_directory=$sample"_dir"
genome_size=$2
reads=$3
#errorRate=$4

module unload java
module load bioinfo-tools canu/1.8

mkdir $output_directory
cd $output_directory 
#/proj/b2014286/private/tools/canu/Linux-amd64/bin/canu -p $sample -d $output_directory genomeSize=$genome_size correctedErrorRate=$errorRate useGrid=false -pacbio-raw $reads

canu -p $sample -d $output_directory genomeSize=$genome_size useGrid=false -pacbio-raw $reads

