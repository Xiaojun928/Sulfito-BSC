#!/bin/bash

for i in `ls *treefile`
do
j=$(echo $i | cut -d "." -f 1)

out1=$(grep -o 'Sulfitobacter_sp_D7_B5M07_[0-9]*_AYE[0-9]*.1' $i)
out2=$(grep -o 'Sulfitobacter_sp_HI0082_A3753_[0-9]*_KZZ[0-9]*.1' $i)

nw_reroot -l -s $i $out1 $out2 > $j.rooted.nwk
nw_labels -I $i > $j.label
done


