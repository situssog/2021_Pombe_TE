# I run it under conda environment
# conda activate

# list.txt is formated like this
# JB22	JB22_finisherSC.fasta
# JB864	JB864_finisherSC.fasta

# Actually, it is the same as 2.longreads_polishing1_with_gcpp.sh

while read sample assembly
do

    echo "-----------------$sample---------------"
    echo "-----------------$assembly"

    ### input
	contig=$assembly  # assembled contig fasta file
    subreads=$sample.bam     # raw long reads file
	
	### output
    mapbam=$sample.pbmm2.bam    # pbmm2 mapping file
    outname=$sample             # output base name

    echo $contig
    echo $subreads
    echo $mapbam
    echo $outname

    samtools faidx $contig
    pbmm2 align $contig $subreads $mapbam --sort -j 8 -J 8 --preset SUBREAD
    gcpp -j 8 -r $contig -o $sample.consensus.fasta,$sample.consensus.vcf,$sample.consensus.gff $mapbam

    sed 's/|//' $sample.consensus.fasta > $outname.gcpp.fasta

done < list.txt

echo "gcpp is done"


