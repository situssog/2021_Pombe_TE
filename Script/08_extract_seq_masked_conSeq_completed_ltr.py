#!/usr/bin/python

# import argparse # For the fancy options
from Bio import SeqIO # to manipulate fasta sequences/files
from Bio.SeqRecord import SeqRecord # to manipulate fasta sequences/files
import sys # To exit the script 
import os # For the input name
import numpy 
import matplotlib.pyplot as plt

# ------------------------------------------------------
fasta_contigs=open(sys.argv[1],"r")

table_annotation=open(sys.argv[2],"r")
table_annotation_line=table_annotation.readlines()

locus=sys.argv[3]
name=sys.argv[4]
min_len=sys.argv[5]


table_summary = open(locus+'_'+name+'_minLen'+min_len+'_table.txt', 'wb')


ID_loc_dic={}
ID_loc=[]

for record in table_annotation_line:
	#print(record)
	sample = record.split("\t")[0]
	start = int(record.split("\t")[1])-1
	end = int(record.split("\t")[2])-1
	length = end-start
	#best_match = record.split("\t")[3]
	if length > int(min_len):
		if sample in ID_loc:
			seqID += 1
			#print(seqID)
			ID_loc_dic[sample+'_'+str(seqID)] = str(start)+'\t'+str(end)+'\t'+str(length)
		else:
			seqID = 1
			ID_loc_dic[sample+'_1'] = str(start)+'\t'+str(end)+'\t'+str(length)
			ID_loc.append(sample)


for sample in ID_loc_dic:
	table_summary.write(sample+'\t'+str(int(ID_loc_dic[sample].split("\t")[0])+1)+'\t'+str(int(ID_loc_dic[sample].split("\t")[1])+1)+'\t'+str(ID_loc_dic[sample].split("\t")[2]))


sequences_dic= []
for seq_record in SeqIO.parse(fasta_contigs, "fasta"):
	seqID = 1
	while seq_record.id+'_'+str(seqID) in ID_loc_dic:
		#print(seq_record.id+'_'+str(seqID))
		start = int(ID_loc_dic[seq_record.id+'_'+str(seqID)].split("\t")[0])
		end = int(ID_loc_dic[seq_record.id+'_'+str(seqID)].split("\t")[1])
		record = SeqRecord(seq_record.seq[start:end], seq_record.id+'_'+str(seqID)+'_'+name, '', '')
		sequences_dic.append(record)
		seqID += 1
		#print(seqID)


#
SeqIO.write(sequences_dic, locus+'_'+name+'_minLen'+min_len+'.fasta', "fasta")
#
table_summary.close()
#

