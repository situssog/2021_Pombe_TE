​### This script is used to run quast for assembly evaluation.

# quast list is formated like this
# JB22.canu.gcpp.finisherSC.gcpp.fasta	pombe_reference.fasta	JB22.canu.gcpp.finisherSC.gcpp
# JB864.canu.gcpp.finisherSC.gcpp.fasta	pombe_reference.fasta	JB864.canu.gcpp.finisherSC.gcpp


while read ass ref out
do 
    echo "===$ass"
    echo "===$ref"
    echo "===$out"


        /home/suofang/Software/quast-5.0.2/quast.py $ass -r $ref -o $out.quast


done < quast.list


echo "finished"
