# This script is used to run wtdbg assembler.

wtdbgroot=/home/suofang/Software/wtdbg-2.4_x64_linux


###### for strains sequenced by PacBio Sequel, wtdbg options "-x sq"
samplelist=(
JB1206.longreads.lenght1000.fasta
JB22.longreads.lenght1000.fasta
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

       $wtdbgroot/wtdbg2 -x sq -L 3000 -g 12.5m -t 8 -i $i -fo wtdbg_3000.$i
       $wtdbgroot/wtpoa-cns -t 8 -i wtdbg_3000.$i.ctg.lay.gz -fo wtdbg_3000.$i.ctg.fa
 
done

###### for strains sequenced by PacBio RS II, wtdbg options "-x rs"
samplelist=(
JB1180.longreads.len1000.fasta
JB4.longreads.len1000.fasta
)

for i in ${samplelist[@]}
do
       echo "-------------------$i-----------------"
       sample=$i

       $wtdbgroot/wtdbg2 -x rs -L 3000 -g 12.5m -t 8 -i $i -fo wtdbg_3000.$i
       $wtdbgroot/wtpoa-cns -t 8 -i wtdbg_3000.$i.ctg.lay.gz -fo wtdbg_3000.$i.ctg.fa
 
done

echo "wtdbg is done"




