#!/bin/bash

# ================== Get fasta sequences out of parsed output =================
# Script to substract a number of sequences out of a fasta file using a file with a list of identifiers inside (one per line)

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014/07/09
# +++++++++++++++++++++++++++++++++++++++++++++++++

listfile=$1
fastafile=$2

if [[ -z "$listfile" ]] ; then 								# -z string, True if the length of string is zero.
	echo "Usage: $0 listfile.txt file.fasta"
	exit 1
fi

if [[ ! -e "$listfile" || ! -s "$listfile" ]] ; then 			# -e file exists, -s file is not zero size
	echo "Could not find list file $listfile or file is size zero"
	exit 1
fi

if [[ ! -e "$fastafile" || ! -s "$fastafile" ]] ; then 			# -e file exists, -s file is not zero size
	echo "Could not find fasta file $fastafile or file is size zero"
	exit 1
fi


for seq in $(cat $listfile)
	do
		python subsetfastaID.py $fastafile $seq
	done
inputname=$(ls $fastafile | awk -F"/" '{print $NF}' - | cut -d'.' -f1)

cat $inputname'_'* > scf'_'mating'_'$fastafile && rm $inputname'_'*