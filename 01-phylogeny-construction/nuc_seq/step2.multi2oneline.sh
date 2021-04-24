#!/bin/bash

for i in `ls *mafft`
do
j=$(basename $i .mafft)
seqtk seq -l0 $i > $j.msa
mv $j.msa $i
done
