# This script is used to run Canu assembler.

canuroot=/home/suofang/software/canu-1.8/Linux-amd64/bin

samplelist=(
JB1180.longreads.lenght1000.fasta
JB1206.longreads.lenght1000.fasta
JB22.longreads.lenght1000.fasta
JB4.longreads.lenght1000.fasta
JB758.longreads.lenght1000.fasta
JB760.longreads.lenght1000.fasta
JB837.longreads.lenght1000.fasta
JB840.longreads.lenght1000.fasta
JB854.longreads.lenght1000.fasta
JB864.longreads.lenght1000.fasta
JB869.longreads.lenght1000.fasta
JB873.longreads.lenght1000.fasta
JB874.longreads.lenght1000.fasta
JB879.longreads.lenght1000.fasta
JB929.longreads.lenght1000.fasta
JB938.longreads.lenght1000.fasta
JB943.longreads.lenght1000.fasta
JB953.longreads.lenght1000.fasta
DY34373.longreads.lenght1000.fasta
DY39827.longreads.lenght1000.fasta
)

for i in ${samplelist[@]}
do
        echo "-------------------$i-----------------"
        sample=$i
        p=$sample"_canu"
        d=$sample"_canu"
        $canuroot/canu -p $p -d $d genomeSize=12.5m useGrid=false -pacbio-raw $i

done

echo "canu is done"
