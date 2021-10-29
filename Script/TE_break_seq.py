#!/usr/bin/python

import os
import sys
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from Bio.Blast.Applications import NcbiblastxCommandline



fasta_contigs=open(sys.argv[1],"r")
annotation_table=open(sys.argv[2],"r")
extra_seq=open(sys.argv[3],"r").readlines()
minlen=int(sys.argv[5])



sequences_dic={}
for seq_record in SeqIO.parse(open(sys.argv[1],"r"), "fasta"):
	sequences_dic[seq_record.id]=seq_record.seq

final_seqs=[]
for line in annotation_table:
	sedID = line.split("\t")[0]
	new_sedID = line.split("\t")[1]
	new_start = int(line.split("\t")[2])-1
	new_end = int(line.split("\t")[3])-1
	new_len = line.split("\t")[4].split("\n")[0]
	if (int(new_len)>minlen):
		record = SeqRecord(Seq(str(sequences_dic[sedID][new_start:new_end])), id=str(new_sedID)+"_"+new_len, description="")
		#print(str(new_sedID)+"_"+new_len)
		final_seqs.append(record)

#print("done1")

for sedID in extra_seq:
	sedID_ed = sedID.split("\n")[0]
	if (len(str(sequences_dic[sedID_ed]))>minlen):
		record = SeqRecord(Seq(str(sequences_dic[sedID_ed])), id=str(sedID_ed)+"_1_"+str(len(sequences_dic[sedID_ed])), description="")
		final_seqs.append(record)
	#print(str(sedID_ed)+"_1_"+str(len(sequences_dic[sedID_ed])))

output_handle = open(sys.argv[4], "w")
SeqIO.write(final_seqs, sys.argv[4], "fasta")
output_handle.close()

