#!/usr/bin/python
import sys
fasta_contigs=open(sys.argv[1],"r")
as_list=open(sys.argv[2],"r")
assigment_list = as_list.readlines()
sample=str(sys.argv[3])

# this script compares contigs with reference genomes a change the name of contigs according aligment coverage

from Bio import SeqIO

contigs_dictionary = {}
list_names = []

for record in assigment_list:
	contig = record.split(" ")[2].split("\n")[0]
	chromosome = record.split(" ")[1]
	coverage = record.split(" ")[0]
	if contig in contigs_dictionary:
		if chromosome in contigs_dictionary[contig]:
			contigs_dictionary[contig][chromosome] += float(coverage)
		else:
			contigs_dictionary[contig][chromosome] = float(coverage)
	else:
		contigs_dictionary[contig] = {chromosome: float(coverage)}


for seq_record in SeqIO.parse(fasta_contigs, "fasta"):
	if seq_record.id in contigs_dictionary:
		new_name = ""
		for chromosome in contigs_dictionary[seq_record.id]:
			if int(contigs_dictionary[seq_record.id][chromosome]) > 2:
				new_name = new_name + '_' + chromosome + '_' + str(int(contigs_dictionary[seq_record.id][chromosome]))
		if new_name == "":
			print ('>' + sample + '_' + seq_record.id + '_' + str(len(seq_record.seq)) + "\n" + str(seq_record.seq) + "\n")
		else:
			print ('>' + sample + '_' + seq_record.id + new_name + '_' + str(len(seq_record.seq)) + "\n" + str(seq_record.seq) + "\n")
	else:
		print ('>' + sample + '_' + seq_record.id + '_' + str(len(seq_record.seq)) + "\n" + str(seq_record.seq) + "\n")


