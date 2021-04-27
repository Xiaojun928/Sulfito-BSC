#!/bin/bash
#SBATCH -n 16

#	usage
mkdir 00_seq_folder
echo "sh this.sh"
echo "require protein seqs namely .pep in 00_seq_folder\n\n"

#	initialization
threadsT=16 #diamond threads
threadsA=16

#	operation
time /your/path/of/OrthoFinder/OrthoFinder-2.2.1/orthofinder.py -og -1 -t $threadsT -a $threadsA -S diamond -M msa -f 00_seq_folder -n 01_blast_folder
