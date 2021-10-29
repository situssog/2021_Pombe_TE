#!/usr/bin/python

# import argparse # For the fancy options
from Bio import SeqIO # to manipulate fasta sequences/files
from Bio.SeqRecord import SeqRecord # to manipulate fasta sequences/files
import sys # To exit the script 
import os # For the input name
import numpy 
import matplotlib.pyplot as plt


# example:
# python extract_seq.py assembly_sample.fasta annotation.txt EBCXXX allSeq 0 
# ------------------------------------------------------
fasta_contigs=open(sys.argv[1],"r")

table_annotation=open(sys.argv[2],"r")
# table_annotation_line=table_annotation.readlines()

sample=sys.argv[3]
name=sys.argv[4]
min_len=sys.argv[5]


sequences_dic={}
for seq_record in SeqIO.parse(fasta_contigs, "fasta"):
	sequences_dic[seq_record.id] = seq_record.seq
	#print(seq_record.id)


table_summary = open(name+'_minLen'+min_len+ '_' + sys.argv[2], 'w')

final_seqs = []
id_series = 1
for record in table_annotation:
	chromosome = record.split("\t")[0]
	start = int(record.split("\t")[3])-1
	end = int(record.split("\t")[4])-1
	annotation_family=record.split("\t")[2]
	direction=record.split("\t")[6]
	if direction=="+":
		ID_seq=str(sample + "_" + str(id_series)) + '_plus_len'+str(end-start)
	else:
		ID_seq=str(sample + "_" + str(id_series)) + '_comp_len'+str(end-start)
	id_series += 1
	table_summary.write(record.split("\n")[0]+'\t'+ ID_seq+ '\n')
	if (end-start>int(min_len)):
		if direction=="+":
			final_seq = sequences_dic[chromosome][start:end]
			seq_name=ID_seq
			record = SeqRecord(final_seq, seq_name, '', '')
			final_seqs.append(record)
		else:
			final_seq = sequences_dic[chromosome][start:end]
			final_seq_rc = final_seq.reverse_complement()
			seq_name=ID_seq
			record = SeqRecord(final_seq_rc, seq_name, '', '')
			final_seqs.append(record)


SeqIO.write(final_seqs, sample+'_'+name+'_minLen'+min_len+'.fasta', "fasta")

table_summary.close()

