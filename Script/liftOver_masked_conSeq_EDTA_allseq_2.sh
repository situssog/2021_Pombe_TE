#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=23:00:00
#SBATCH -J liftOver

sample=$1
annotation_file=$2

#mkdir dir_ext_masked_conSeq_$sample"_"liftover
sourcefa=$PWD"/"$sample.fasta
targetfa=$PWD"/"ref_masked_conSeq.fa
cd dir_ext_masked_conSeq_$sample"_"liftover

echo "translate coordenates per scafold"

# #ltr_finder_table=/dss/dsslegfs01/pr53da/pr53da-dss-0022/nobackup/private/pac_bio/02_TE/all21_samples/LTR_finder/$sample.fasta.finder.combine.scn 
# #repeatmasker_table=/dss/dsslegfs01/pr53da/pr53da-dss-0022/nobackup/private/pac_bio/02_TE/all21_samples/LTR_RepeatMasker/$sample.fasta.ori.out 
# grep "^\[" $ltr_finder_table | awk -v OFS='\t' '{print $2, $3}' | sed 's/-/\t/g' > ltr_finder_table_$sample.bed
# awk -v OFS='\t' '{print $5, $6, $7}' $repeatmasker_table > tem_ltr_repeatMasker_table_$sample.bed
# module load BEDTools/2.27.1
# bedtools merge -i tem_ltr_repeatMasker_table_$sample.bed > ltr_repeatMasker_table_$sample.bed
# rm tem_ltr_repeatMasker_table_$sample.bed
# bedtools intersect -v -a ltr_repeatMasker_table_$sample.bed -b ltr_finder_table_$sample.bed > $sample"_"annotation1.bed
# # bedtools intersect -wa -a ltr_finder_table_$sample.bed -b ltr_repeatMasker_table_$sample.bed > $sample"_"annotation2.bed
# bedtools intersect -wo -a ltr_finder_table_$sample.bed -b ltr_repeatMasker_table_$sample.bed | awk -v OFS='\t' '{if($7>500) print $1, $2, $3}' | sort | uniq > $sample"_"annotation2.bed
# cat $sample"_"annotation2.bed $sample"_"annotation1.bed | awk -v OFS='\t' '{print $0, NR"_""'$sample'""_"$1"_"$2"_"$3}' > $sample"_"annotation_ed.bed 
# # cat ltr_finder_table_$sample.bed $sample"_"annotation.bed | awk -v OFS='\t' '{print $0, NR"_""'$sample'""_"$1"_"$2"_"$3}' > $sample"_"annotation_ed.bed 

grep -v "^#" $annotation_file | awk '{print $1}' | sort | uniq > $sample'_'scafolds.txt
grep -v "^#" $annotation_file > $sample'_'repeats_annotation_file.gff3

#awk '{print $1}' $sample"_"annotation_ed.bed | sort | uniq > $sample'_'scafolds.txt

#module load R/3.4.3 MariaDB/10.2.11
#R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4

for scafold_ID in $( cat $sample'_'scafolds.txt)
do
# grep $scafold_ID"__" chain_table_$sample"_"refCoor.bed > $scafold_ID"__"chain_table_$sample"_"refCoor.bed
grep -P "\t"$scafold_ID"__" chain_table_$sample"_"refCoor.bed > $scafold_ID"__"chain_table_$sample"_"refCoor.bed
done

for scafold_ID in $( cat $sample'_'scafolds.txt)
do
Rscript --vanilla /dss/dsshome1/lxc03/di36guz2/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/annotate_coordenates.R $sample $scafold_ID"__"chain_table_$sample"_"refCoor.bed $sample'_'repeats_annotation_file.gff3 
done
echo "Done!"


