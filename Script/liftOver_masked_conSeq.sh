#!/bin/bash
#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J liftOver

module load ruby/2.6.2

export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/bin
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/gems
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/gems/bio-1.5.2
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems/gems/bio-1.5.2/bin
export PATH=$PATH:/proj/uppstore2017159/b2014286/private/tools/gems

sample=$1

mkdir dir_ext_masked_conSeq_$sample"_"liftover
sourcefa=$PWD"/"$sample"_"ed.fa
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
grep -vw "family.*_consensus " /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/02_TE/all21_samples/CARP/dir_$sample/results_classify/final/$sample""_annotation.fa.out | grep -v "^$" | grep -v "class/family" | grep -v "repeat" | awk '{print $5, "\t",$6-2, "\t", $7+2, "\t", $10"_""'$sample'""_"NR}' | sed 's/family......_consensus\://g' > $sample"_"annotation.bed
awk '{print $1}' $sample"_"annotation.bed | sort | uniq > $sample'_'scafolds.txt
module load R/3.4.3 MariaDB/10.2.11
R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4
for scafold_ID in $( cat $sample'_'scafolds.txt)
do
grep $scafold_ID"__" chain_table_$sample"_"refCoor.bed > $scafold_ID"__"chain_table_$sample"_"refCoor.bed
Rscript --vanilla /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/annotate_coordenates.R $sample $scafold_ID"__"chain_table_$sample"_"refCoor.bed
done
echo "Done!"

