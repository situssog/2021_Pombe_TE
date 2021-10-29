#!/usr/bin/python

import os
import sys

annotation=open(sys.argv[1],"r")
core_repeatMask=open(sys.argv[2],"r")
Fltr_repeatMask=open(sys.argv[3],"r")

seq_Core=[]
for record in core_repeatMask:
	ID = record.split("\t")[4]
	if not ID in seq_Core:
		seq_Core.append(ID)

def overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return (
        start1 <= start2 <= end1 or
        start1 <= end2 <= end1 or
        start2 <= start1 <= end2 or
        start2 <= end1 <= end2
    )

seq_LTR={}
for record in Fltr_repeatMask:
	ID = record.split("\t")[4]
	start_blast = int(record.split("\t")[5])
	end_blast = int(record.split("\t")[6])
	left_seq = int(record.split("\t")[7].replace('(', '').replace(')', ''))
	if ID in seq_LTR.keys():
		LTR_start1=seq_LTR[ID][0]
		LTR_start2=seq_LTR[ID][1]
		LTR_end1=seq_LTR[ID][2]
		LTR_end2=seq_LTR[ID][3]
		if not overlap(LTR_start1, LTR_start2, start_blast, end_blast):
			if start_blast < 100:
				if start_blast < LTR_start1:
					LTR_start1 = start_blast
					LTR_start2 = end_blast
			elif left_seq < 100:
				if end_blast > LTR_end2:
					LTR_end1=start_blast
					LTR_end2=end_blast
		seq_LTR[ID]=[LTR_start1,LTR_start2,LTR_end1, LTR_end2]
	else:
		if start_blast < 100:
			seq_LTR[ID]=[start_blast, end_blast, 0, 0]
		elif left_seq < 100:
			seq_LTR[ID]=[0,0,start_blast, end_blast]


for record in annotation:
	line = record.split("\n")[0]
	ID = line.split("\t")[-1]
	if ID in seq_Core:
		core_found = 1
	else:
		core_found = 0
	if ID in seq_LTR.keys():
		if seq_LTR[ID][0] == 0:
			ltr_start_found = 0
		else:
			ltr_start_found = 1
		if seq_LTR[ID][2] == 0:
			ltr_end_found = 0
		else:
			ltr_end_found = 1
	else:
		ltr_start_found = 0
		ltr_end_found = 0
	print(line + "\t" + str(ltr_start_found) + "\t" + str(core_found) + "\t" + str(ltr_end_found))





