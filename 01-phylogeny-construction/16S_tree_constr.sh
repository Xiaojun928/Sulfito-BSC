#!/bin/bash

INPUT="16S_rRNA_gene_30.fasta"  #The input can be changed as '16S_rRNA_all.fasta'
einsi  $INPUT > 16S_rRNA.mafft
sed -i "s/|/ /g" 16S_rRNA.mafft
trimal -in 16S_rRNA.mafft -out Sulfitobacter_16S.trimal -automated1 #-resoverlap 0.55 -seqoverlap 60
iqtree -s Sulfitobacter_16S.trimal -nt 1 -m MFP -mrate E,I,G,I+G -wbtl -mfreq F,F -bb 1000 -alrt 1000 -pre Sulfito_16S -redo
