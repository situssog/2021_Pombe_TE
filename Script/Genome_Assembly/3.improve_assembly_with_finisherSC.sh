# This script is used to run finisherSC

samplelist=(
DY34373
DY39827
JB1180
JB1206
JB22
JB4
JB758
JB760
JB837
JB840
JB854
JB864
JB869
JB873
JB874
JB879
JB929
JB938
JB943
JB953
)

########### for each of strain, do this
for i in ${samplelist[@]}
do

mkdir $i
cd $i

# following the finisherSC github description 
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' /data/LongReads/$i.fasta.1000.fasta > raw_reads.fasta
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' /data/Assembly/Gcpp/$i.gcpp.fasta > contigs.fasta
cd ..

# run finisherSC, the output file is "improved3.fasta".
python /home/suofang/Software/finishingTool-2.1/finisherSC.py $i /home/suofang/Software/MUMmer3.23

done

