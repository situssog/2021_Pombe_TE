#!/usr/bin/python

# import argparse # For the fancy options
from Bio import SeqIO # to manipulate fasta sequences/files
from Bio.SeqRecord import SeqRecord # to manipulate fasta sequences/files
import sys # To exit the script 
import os # For the input name
import numpy 
import matplotlib.pyplot as plt


# example:
# python extract_seq_withCoor.py file.fasta chr start end out 
# ------------------------------------------------------
fasta_contigs=open(sys.argv[1],"r")

chr_name=str(sys.argv[2])

start_pos=int(sys.argv[3])

end_pos=int(sys.argv[4])

out=str(sys.argv[5])


sequences_dic={}
final_seqs = []
for seq_record in SeqIO.parse(fasta_contigs, "fasta"):
	if (seq_record.id==chr_name):
		final_seq = seq_record.seq[start_pos-1:end_pos-1]
		seq_name=seq_record.id+'_start_'+str(start_pos)+'_end_'+str(end_pos)+'_len_'+str(end_pos-start_pos)
		record = SeqRecord(final_seq, seq_name, '', '')
		final_seqs.append(record)

SeqIO.write(final_seqs, out+".fasta", "fasta")



