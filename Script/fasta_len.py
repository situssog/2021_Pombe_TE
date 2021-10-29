#!/usr/bin/python

import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'r')

for rec in SeqIO.parse(FastaFile, 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    print(name+"\t"+str(seqLen))

#    print name, seqLen
#FastaFile.close()




