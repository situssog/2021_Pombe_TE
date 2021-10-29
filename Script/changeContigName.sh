#This script change take a draft de-novo assembly, identify the mitochondrial contig, trims it leaving only one copy, and change the names of all contigs. The new name has percentage con coverage of reference genome, identifying reference chromosome. 

# The script uses the reference genome in /proj/b2014286/private/nobackup/pac_bio/03_denovo_assembly/natural_strains_group1/
# /proj/b2014286/private/nobackup/pac_bio/03_denovo_assembly/natural_strains_group1/MTRef.fa
# /proj/b2014286/private/nobackup/pac_bio/03_denovo_assembly/natural_strains_group1/PombeRef_withAB325691.fa

# use: example: /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/11_trim_MT_changeContigName.sh EBC069

# Blast of contigs with reference to identify MT contig
sample=$1
# module load bioinfo-tools blast/2.2.29+ biopython 
module load bioinfo-tools blast/2.9.0+ biopython
blastn -query $sample".fasta" -db /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/03_denovo_assembly/natural_strains_group1/MTRef.fa -outfmt 6 | awk '{ if ($9 <= 10 || $10 <= 10) print $0 }'  > $sample"_MTblast_scf".txt

# extract MT contig and nuclear contigs in different files:
#python /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/08_extract_MT_break_singlerepeat.py $sample"_consensus".fa $sample"_MTblast_scf".txt MT > $sample"_MT".fa 
#sed -i "s/>/>$sample"_"/g" $sample"_MT".fa
#python /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/08_extract_MT_break_singlerepeat.py $sample"_consensus".fa $sample"_MTblast_scf".txt Nucl > $sample"_Nucl".fa 

# re-name nuclear contigs:
module load bioinfo-tools MUMmer/3.23
#nucmer --mum  -c 100 --prefix=ref_align_$sample"_forName" /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/03_denovo_assembly/natural_strains_group1/PombeRef_withAB325691.fa $sample"_Nucl".fa
nucmer --mum  -c 100 --prefix=ref_align_$sample"_forName" /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/03_denovo_assembly/natural_strains_group1/PombeRef_withAB325691.fa $sample".fasta"
show-coords -c ref_align_$sample"_forName".delta > ref_align_$sample"_forName".coords
# I filtered the output leaving only alignment fragments largers than 10000 bp
awk '{ if ($7 >= 10000) print $13, $15, $16 }' ref_align_$sample"_forName".coords | grep -e "^[0-9]" - > $sample"_contigs_assigment".txt
# the next script changes the name of the contigs according to the coverage on the reference genome. 
# python 07_assigment_contings.py "assembly" "EBC069_contigs_assigment.txt" "sample_name"
python /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/09_TE/scripts_all21Samples/07_assigment_contings.py $sample".fasta" $sample"_contigs_assigment".txt $sample | sed '/^$/d' > $sample"_ed".fa


# merge MT and nuclear contigs in a single file:
#cat $sample"_Nucl_ed".fa $sample"_MT".fa > $sample"_consensus_ed".fa


