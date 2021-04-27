#!/bin/bash

#	usage
echo "sh this.sh"
echo "protein seqs namely .faa are reuired in 00_genome_info\n\n"
echo "protein seqs of RMs should be downloaded from REBASE http://rebase.neb.com/rebase/rebase.ftp.html"

/your/path/of/diamond/diamond makedb --in /your/path/of/REBASE/protein_seqs.fasta -d RM_db


for i in `cat gnm_order.txt`
do
j=$(echo $i | cut -d "_" -f 3)
/your/path/of/diamond/diamond blastp -q ../00_genome_info/${i}.faa --db /your/path/of/REBASE/RM_db.dmnd --outfmt 6 --max-target-seqs 0 --evalue 1e-5 --id 30 --query-cover 50 -o ${j}_diamond.txt
done
