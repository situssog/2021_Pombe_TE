#!/bin/bash
#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 200:00:00
#SBATCH -J RaxMl_pombe

# ================== Molecular systematics =================
# Doing a full ML analysis in Uppmax.s

# This is a simplified version of RAxML_FullAnalyses_partUppmax.sh to do a ML
# analysis for an alignment without partitions and with relatively little
# computation effort. This is meant to be run within a loop when having MANY
# genes that are NOT concatenated in the script RAxML_GeneByGene.sh.

# The script takes an alignment in phylip (or fasta)
# format as input, but it needs the entire absolute path if its in another
# folder than the current one. It makes new folders in the working folder. In
# summary, it searches the best ML tree out of given number of random starting trees and
# randomized Maximum Parsimony starting trees; it produces bootstraps
# and then it maps the bootstraps to the best ML tree found.

# RAxML program should be in the path as raxmlHPC. Remember RAxML doesn't take
# relative paths, only absolute.

# To run RAxML just to get the alignment: the key is -f c
# module load bioinfo-tools
# module load raxml/8.0.20-mpi 
# $ raxmlHPC-PTHREADS-AVX -f c -m GTRGAMMAI -s alignment.fas -N 1000 -b 11111 -n reduced -p 11111 -T 2 # -T is obligatory if using the threaded version and it has to be at least 2

# ==================================================
# Sandra Lorena Ament-Velasquez
# Uppsala, Sweden, 2016/02/13
# +++++++++++++++++++++++++++++++++++++++++++++++++


# set -e # Exit/fail if any of the commands that are run within the script fail (ie. they return a non-zero exit code)

# The modules are loaded in RAxML_GeneByGene.sh
module load bioinfo-tools
module load raxml/8.2.10-gcc-mpi
#module load raxml/8.0.20-mpi 

echo "--------------------------"
date


alignment=$1
CORES=16  # No. of cores used. On Milou there are 16 cores per node. An example, specifying 2 nodes, and thus 32 (2 * 16) cpus, would be -n 32 (above)
RAXML=raxmlHPC-PTHREADS-AVX   # raxmlHPC

if [[ -z "$alignment" ]] ; then 								# -z string, True if the length of string is zero.
	echo "Usage: `basename $0` path/to/alignment.phy"
	echo "--------------------------"
	exit 1
fi

if [[ ! -e "$alignment" || ! -s "$alignment" ]] ; then 			# -e file exists, -s file is not zero size
	echo "Could not find phylip/fasta file $alignment or file is size zero"
	echo "--------------------------"
	exit 1
fi

# Set how much effort do you want to put
searches=30
boots=1000


echo "The program was called like this: $ bash $0 $alignment"
echo
name=$(basename $alignment | cut -d'.' -f1)
workpath=$(pwd)

# ---
# * 5.2.3 Finding the Best-Known Likelihood tree (BKL)
# ---
# In addition to bootstrapping, you should do multiple inferences on the
# original alignment.  So, to execute a multiple inference on the original
# alignment on a single processor just specify:
#mkdir RandMP_inferences
#$RAXML -f d -m GTRGAMMAI -s $alignment -N $searches -n MultipleOriginal -w $workpath/RandMP_inferences -p 1111 -T $CORES

# This computes 100 inferences with the rapid default climbing algorithm. This
# should let me see if there are clusters of LogL (islands)

# Use random starting trees
#mkdir Rand_inferences
#$RAXML -d -m GTRGAMMAI -s $alignment -N $searches -n MultipleOriginal -w $workpath/Rand_inferences -p 1111 -T $CORES

# Produce a file for R to check the distribution
#echo "=== Extracting optimized LogL values from the resulting RAxML_info files ==="
#mkdir LogL_islands
#cat RandMP_inferences/RAxML_info.MultipleOriginal | grep "GAMMA-based Likelihood" | cut -d' ' -f 5 > LogL_islands/LogL_values_RMP_$name.txt
#cat Rand_inferences/RAxML_info.MultipleOriginal | grep "GAMMA-based Likelihood" | cut -d' ' -f 5 > LogL_islands/LogL_values_R_$name.txt

# ---
# * 5.2.4 Bootstrapping with RAxML
# ---
# echo "=== Starting bootstraping ==="
# This is for generate 1000 bootstrap trees that can later be mapped to the
# Best Tree. I decided to use GTRGAMMAI because in my experience using Blackbox
# millions of times, the BT values are very similar, and the search is faster.
mkdir Bootstraps
cd Bootstraps
$RAXML -f d -m GTRGAMMAI -s $alignment -N $boots -b 11111 -n $name'_MultipleBootstrap.RAxML' -p 11111 -T $CORES  
cd ..

# ---
# *5.2.5 Obtaining Confidence Values
# ---
mkdir ML_Best_tree
# Save LogLikelihood values of MP and random starting trees searches.
randML=$(cat Rand_inferences/RAxML_info.MultipleOriginal | grep "Final GAMMA-based Score of best tree" | cut -d' ' -f 7)
randMPML=$(cat RandMP_inferences/RAxML_info.MultipleOriginal | grep "Final GAMMA-based Score of best tree" | cut -d' ' -f 7)

echo "randML logL = "$randML
echo "randMPML logL = "$randMPML

# Check which one has the best LogLikelihood value
if [ $(bc <<< "$randML <= $randMPML") -eq 1 ] # (<<<) notation, as a nice alternative to echo "$result <= $key1" | bc.
	then
		echo "The best tree comes from the MP random trees"
		cp $workpath/RandMP_inferences/RAxML_bestTree.MultipleOriginal ML_Best_tree/RAxML_bestTree_$name.MultipleOriginal_rMP
		cd ML_Best_tree
		$RAXML -f b -m GTRGAMMAI -s $alignment -z $workpath/Bootstraps/RAxML_bootstrap.$name'_MultipleBootstrap.RAxML' -t RAxML_bestTree_$name.MultipleOriginal_rMP -n $name'_MLfinal_Tree.tre' -T $CORES
		cd ..
	else
		echo "The best tree comes from the completely random trees"
		cp $workpath/Rand_inferences/RAxML_bestTree.MultipleOriginal ML_Best_tree/RAxML_bestTree_$name.MultipleOriginal_rand
		cd ML_Best_tree
		$RAXML -f b -m GTRGAMMAI -s $alignment -z $workpath/Bootstraps/RAxML_bootstrap.$name'_MultipleBootstrap.RAxML' -t RAxML_bestTree_$name.MultipleOriginal_rand -n $name'_MLfinal_Tree.tre' -T $CORES 
		cd ..	
fi

echo
echo "Done!"

date
echo "--------------------------"
