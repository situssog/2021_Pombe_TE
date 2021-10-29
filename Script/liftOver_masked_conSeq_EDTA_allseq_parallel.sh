#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=12:00:00
#SBATCH -J liftOver

sample=$1
annotation_file=$2

#mkdir dir_ext_masked_conSeq_$sample"_"liftover
sourcefa=$PWD"/"$sample.fasta
targetfa=$PWD"/"ref_masked_conSeq.fa
cd dir_ext_masked_conSeq_$sample"_"liftover

echo "translate coordenates per scafold"

grep -v "^#" $annotation_file | awk '{print $1}' | sort | uniq > $sample'_'scafolds.txt
grep -v "^#" $annotation_file > $sample'_'repeats_annotation_file.gff3

cat $sample'_'scafolds.txt | parallel -j 10 -I{} grep -P "'\t'"{}"__" chain_table_$sample"_"refCoor.bed ">" {}"_test__"chain_table_$sample"_"refCoor.bed
cat $sample'_'scafolds.txt | parallel -j 10 -I{} Rscript --vanilla /dss/dsshome1/lxc03/di36guz2/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/annotate_coordenates_parallel.R $sample {}"__"chain_table_$sample"_"refCoor.bed $sample'_'repeats_annotation_file.gff3

echo "Done!"
echo "Merging"
cat parallel_final_* > all_parallel_final_$sample"_"refCoor.txt
echo "Done!"


