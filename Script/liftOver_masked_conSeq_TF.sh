#!/bin/bash
#SBATCH -A snic2018-3-658
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J liftOver

module load ruby/2.6.2
module load bioinfo-tools blat/36

export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/bin
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/gems
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/gems/bio-1.5.2
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/gems/bio-1.5.2/bin
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems

sample=$1
ltr_finder_table=$2
repeatmasker_table=$3

mkdir dir_ext_masked_conSeq_$sample"_"liftover
sourcefa=$PWD"/"$sample.fasta
targetfa=$PWD"/"ref_masked_conSeq.fa
cd dir_ext_masked_conSeq_$sample"_"liftover
echo "Produce input files for liftOver"
cp /proj/uppstore2017159/b2014286/private/tools/flo/flo/opts_example_SITG.yaml ./flo_opts.yaml
sed -i 's/:source_fa: .unset./\:source_fa\: '"'"''$(echo $sourcefa | sed 's/\//\\\//g')''"'"'/' flo_opts.yaml
sed -i 's/:target_fa: .unset./\:target_fa\: '"'"''$(echo $targetfa | sed 's/\//\\\//g')''"'"'/' flo_opts.yaml
sed -i 's/:processes: .unset./\:processes\: '"'"'1'"'"'/' flo_opts.yaml
rake -f /proj/uppstore2017159/b2014286/private/tools/flo/flo/Rakefile
echo "Done!"

echo "Produce table with all changes in coordenates"
module load bioinfo-tools python/2.7.15 biopython/1.73
python /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/list_chr_len.py $sourcefa $sample
echo "Done!"

echo "Run LiftOver"
module load bioinfo-tools liftOver/2017-03-14
liftOver -minMatch=0.7 $sample"_"chr_pos_table.bed run/liftover.chn chain_table_$sample"_"refCoor.bed chain_table_$sample"_"unMapped.txt
echo "Done!"

echo "translate coordenates per scafold"

#ltr_finder_table=/proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/all21_samples/LTR_finder/$sample.fasta.finder.combine.scn 
#repeatmasker_table=/proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/all21_samples/LTR_RepeatMasker/$sample.fasta.ori.out 

grep "^\[" $ltr_finder_table | awk -v OFS='\t' '{print $2, $3}' | sed 's/-/\t/g' > ltr_finder_table_$sample.bed
awk -v OFS='\t' '{print $5, $6, $7}' $repeatmasker_table > ltr_repeatMasker_table_$sample.bed

module load BEDTools/2.27.1
bedtools intersect -v -a ltr_repeatMasker_table_$sample.bed -b ltr_finder_table_$sample.bed > $sample"_"annotation1.bed
# bedtools intersect -wa -a ltr_finder_table_$sample.bed -b ltr_repeatMasker_table_$sample.bed > $sample"_"annotation2.bed
bedtools intersect -wo -a ltr_finder_table_$sample.bed -b ltr_repeatMasker_table_$sample.bed | awk -v OFS='\t' '{if($7>500) print $1, $2, $3}' | sort | uniq > $sample"_"annotation2.bed

cat $sample"_"annotation2.bed $sample"_"annotation1.bed | awk -v OFS='\t' '{print $0, NR"_""'$sample'""_"$1"_"$2"_"$3}' > $sample"_"annotation_ed.bed 

# cat ltr_finder_table_$sample.bed $sample"_"annotation.bed | awk -v OFS='\t' '{print $0, NR"_""'$sample'""_"$1"_"$2"_"$3}' > $sample"_"annotation_ed.bed 

awk '{print $1}' $sample"_"annotation_ed.bed | sort | uniq > $sample'_'scafolds.txt

module load R/3.4.3 MariaDB/10.2.11
R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4

for scafold_ID in $( cat $sample'_'scafolds.txt)
do
# grep $scafold_ID"__" chain_table_$sample"_"refCoor.bed > $scafold_ID"__"chain_table_$sample"_"refCoor.bed
grep -P "\t"$scafold_ID"__" chain_table_$sample"_"refCoor.bed > $scafold_ID"__"chain_table_$sample"_"refCoor.bed
Rscript --vanilla /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/annotate_coordenates.R $sample $scafold_ID"__"chain_table_$sample"_"refCoor.bed $sample"_"annotation_ed.bed 
done
echo "Done!"

