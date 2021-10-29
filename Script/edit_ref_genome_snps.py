#!/usr/bin/python
import sys
fasta_contigs=open(sys.argv[1],"r")
gaps_info=open(sys.argv[2],"r")
gaps_info_line = gaps_info.readlines()
snps_info=open(sys.argv[3],"r")
snps_info_line = snps_info.readlines()

from Bio import SeqIO
import numpy

edited_seq_dir = {}

for seq_record in SeqIO.parse(fasta_contigs, "fasta"):
	#print (seq_record.id, str(seq_record.seq[0:10]))
	working_sequence = list(str(seq_record.seq))
	for gap in gaps_info_line:
		chromosome = gap.split("	")[0]
		break_min = int(min(gap.split("	")[2:4]))-1
		break_max = int(max(gap.split("	")[2:4]))-1
		#print (chromosome, break_min, break_max)
		#print ("".join(working_sequence)[break_min-5:break_max+5])
		if chromosome == seq_record.id:
			working_sequence[break_min:break_max] = list('N'*(break_max-break_min))
			#print ("".join(working_sequence)[break_min-5:break_max+5])
	for snp in snps_info_line:
		chromosome_snp = snp.split("	")[10]
		pos = int(snp.split("	")[0])-1
		ref = snp.split("	")[1]
		alt = snp.split("	")[2]
		if chromosome_snp == seq_record.id:
			if working_sequence[pos] != 'N':
				if alt == ".":
					#print (working_sequence[pos], ref, alt)
					working_sequence[pos] = 'N'
					#print (working_sequence[pos], ref, alt)
				else:
					working_sequence[pos] = alt
	edited_seq_dir[seq_record.id] = working_sequence
	print ('>' + seq_record.id + "\n" + str("".join(working_sequence)))





