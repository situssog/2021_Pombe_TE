#!/usr/bin/python

# ================== Subset fastas using IDs within sequences names =================
# Script to produce a subset of an input fasta file using one or several
# strings that identify specific sequences. The script also reports the number
# of sequences in the input file. This version of the script uses a class for
# load the sequences into a dictionary.

# Because the script is based on dictionaries, it can be used also to get rid
# of repeated sequences. Just call it:
# python subsetfastaID.py file.fas ''
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014/05/21
# +++++++++++++++++++++++++++++++++++++++++++++++++
# python subsetfastaID.py file.fasta cosito 954

# ------------------------------------------------------
import sys
import os
# ------------------------------------------------------

class FastaDic:
	"""
	A group of functions to pass a fasta file into a dictionary, print it back into a file and subset it using a string.
	To use it, first create an object with an empty dictionary, then use fast2dic to load a fasta file into it. For example:

	trial = FastaDic({})
	print trial
	trial.fast2dic(fastafile)
	print trial
	trial.printfastadic('test1')
	trial2 = trial.subsetfasta('96a')
	print trial2

	"""
	def __init__(self, dictionary):
		self.fasta = dictionary
		self.numberseq = len(self.fasta.keys())
	# ---------------------------------
	def __str__(self):
		return "Number of sequences: %d" % self.numberseq
	# ---------------------------------
	def fast2dic(self, fastafile):
		fasta_multi = [ x.replace("\n","") for x in open(fastafile,"r").readlines() ] # make a list of the lines in the file (without the \n at the end of each)

		fasta = dict()
		for i in range(0,len(fasta_multi)):
			if len(fasta_multi[i]) > 0:   							# Avoid empty lines
				if fasta_multi[i][0]==">": 							# Put sequence name in the dictionary as a key
					fasta_multi[i]=fasta_multi[i].replace(">","")	# Get rid of the > sign first
					fastakey = fasta_multi[i] 						# Keep it as the working key
					fasta[fastakey] = "" 							# Notice that all the names of the sequences have to be different or they get overwritten
				else:
					fasta[fastakey] = fasta[fastakey] + fasta_multi[i] # Concatenate the sequence for each key
		self.fasta = fasta
		self.numberseq = len(self.fasta.keys())
	# ---------------------------------
	def printfastadic(self, outputname):
		ofile = open(outputname+'.fas', 'w')
		for seq in self.fasta.keys():
			ofile.write(('>'+seq+'\n'))
			ofile.write((self.fasta[seq]+'\n'))
		ofile.close()
	# ---------------------------------
	def subsetfasta(self, string_sub):
		subset = {}
		for name in self.fasta.keys():
			# if string_sub in name:		 # This would find "file1" in the keys file1, file10, file100
			if string_sub == name: 			 # The key has to be exact, so it will only find "file1"
				subset[name] = self.fasta[name]
		return FastaDic(subset)
	# ---------------------------------

# Input from console
try:
	fastafile = sys.argv[1]
	print "Input file:", fastafile
except:
	print "Usage: python subsetfastaID.py file.fasta string [ string2 ] ..."
	sys.exit(1)


originalfasta = FastaDic({})
originalfasta.fast2dic(fastafile)
print originalfasta

# Taking out the prefix of the file
prefixlen = len(os.path.basename(fastafile).split(".")[-1])
namelen = len (os.path.basename(fastafile))
input_name = os.path.basename(fastafile)[:namelen - prefixlen - 1]

for string in sys.argv[2:]:
	print "String to subset fasta file:", string
	sub = originalfasta.subsetfasta(string)
	print "\t", sub
	if sub.numberseq != 0:
		sub.printfastadic((input_name + '_' + string))
	else:
		print "\tNo fasta file printed"