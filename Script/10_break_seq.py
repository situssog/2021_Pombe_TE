#!/usr/bin/python

import os
import sys
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from Bio.Blast.Applications import NcbiblastxCommandline



# example:
# python extract_seq.py assembly_sample.fasta annotation.txt wtf_found.txt all_annotation_refCoor.bed EBCXXX wtf 500 
# ------------------------------------------------------
#fasta_contigs=open(sys.argv[1],"r")
fasta_contigs=sys.argv[1]
consensus_seq=sys.argv[2]

line_reference_db = 'makeblastdb -in ' + consensus_seq + ' -dbtype nucl'
os.system(line_reference_db)
blastn_cline = NcbiblastxCommandline(cmd="blastn", \
            query=fasta_contigs, db=consensus_seq, \
            evalue=0.00000001, outfmt='"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send slen evalue bitscore"')
blast_results = str(blastn_cline()[0]).split('\n')

sequences_dic={}
seqIDs=[]
for seq_record in SeqIO.parse(open(sys.argv[1],"r"), "fasta"):
	sequences_dic[seq_record.id]=seq_record.seq
	seqIDs.append(seq_record.id)


final_seqs=[]
for seq in seqIDs:
	#print(seq)
	#print(seqIDs)
	n=int(1)
	for line in blast_results:
		if line == '':
			continue
		else:
			temporal_line = line.split('\t')
			seqID = temporal_line[0]
			start = int(min([int(temporal_line[6]),int(temporal_line[7])]))
			end = int(max([int(temporal_line[6]),int(temporal_line[7])]))
			lengthBlast = int(temporal_line[3])
			if lengthBlast > 600 and seqID==seq:
				print(seqID)
				#print(seq)
				#print(str(start-1)+'_'+str(end-1))
				record = SeqRecord(Seq(str(sequences_dic[seqID][start-1:end-1])), id=str(seqID)+'_'+'len'+str(end-start)+'_'+str(n), description="")
				#print(sequences_dic[seqID][start-1:end-1])
				print(str(end-start))
				print(str(n))
				final_seqs.append(record)
				n += 1

output_handle = open(sys.argv[3], "w")
SeqIO.write(final_seqs, sys.argv[3], "fasta")
output_handle.close()

