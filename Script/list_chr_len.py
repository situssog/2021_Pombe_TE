#!/usr/bin/python

# import argparse # For the fancy options
from Bio import SeqIO # to manipulate fasta sequences/files
from Bio.SeqRecord import SeqRecord # to manipulate fasta sequences/files
import sys # To exit the script
import os # For the input name
import numpy
import matplotlib.pyplot as plt


# example:
# python extract_seq.py assembly_sample.fasta annotation.txt wtf_found.txt all_annotation_refCoor.bed EBCXXX wtf 500
# ------------------------------------------------------
fasta_contigs=open(sys.argv[1],"r")
sample=sys.argv[2]

table_summary = open(sample+'_chr_pos_table.bed', "w")

sequences_dic={}
for seq_record in SeqIO.parse(fasta_contigs, "fasta"):
	print(seq_record.id)
	for i in range(len(seq_record.seq)-1):
		table_summary.write(seq_record.id+'\t'+str(i+1)+'\t'+str(i+2)+'\t'+seq_record.id+'__'+str(i+1)+'\n')

table_summary.close()
