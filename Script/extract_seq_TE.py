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

table_annotation=open(sys.argv[2],"r")
table_annotation_line=table_annotation.readlines()

## table_annotation_included=open(sys.argv[3],"r")
## table_annotation_included_list=table_annotation_included.readlines()

table_know_coor=open(sys.argv[3],"r")
table_know_coor_list=table_know_coor.readlines()

sample=sys.argv[4]
name=sys.argv[5]
min_len=sys.argv[6]

## table_summary = open(sample+'_'+name+'_minLen'+min_len+'_table.txt', 'wb')


## annotation_included = []
## for record in table_annotation_included_list:
	## annotation_included.append(record.split("\n")[0])
## #print(annotation_included)

sequences_dic={}
for seq_record in SeqIO.parse(fasta_contigs, "fasta"):
	sequences_dic[seq_record.id] = seq_record.seq
	#print(seq_record.id)


ID_loc_dic={}
for record in table_know_coor_list:
	line = record.split("\n")[0]
	ID=line.split("\t")[0]
	location=line.split("\t")[1]
	ID_loc_dic[ID] = location


final_seqs = []
for record in table_annotation_line:
	chromosome = record.split("\t")[13]
	start = int(record.split("\t")[14])-1
	end = int(record.split("\t")[15].split("\n")[0])-1
	seq_chr = record.split("\t")[0]
	seq_start = int(record.split("\t")[1])
	seq_end = int(record.split("\t")[2])
	## annotation_family=record.split("\t")[3]
	## direction=record.split("\t")[11]
	ID_seq=record.split("\t")[3]
	seq_num=ID_seq.split("_")[0]
	if (seq_end-seq_start>int(min_len)):
		## if annotation_family in annotation_included:
		## if direction=="+":
		final_seq = sequences_dic[seq_chr][seq_start-1:seq_end-1]
		if (ID_seq in ID_loc_dic):
			seq_name=seq_num+"_"+sample+'_'+str(ID_loc_dic[ID_seq])
			## table_summary.write(seq_num+"_"+sample+'\t'+str(ID_loc_dic[ID_seq]).split("_")[0]+'\t'+str(ID_loc_dic[ID_seq]).split("_")[1]+'\t'+'\tplus\t'+str(end-start)+'\n')
		else:
			print("No in sim annotation: "+ID_seq)
		record = SeqRecord(final_seq, seq_name, '', '')
		final_seqs.append(record)


SeqIO.write(final_seqs, sample+'_'+name+'_minLen'+min_len+'.fasta', "fasta")

## table_summary.close()

