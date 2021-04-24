#!/bin/bash
#SBATCH -n 1


for i in `ls *.mafft`
do
		j=$(basename $i .mafft)
		#time trimal -in $i -out ${j}.trimai -automated1 -resoverlap 0.55 -seqoverlap 60
		trimal -in $i -out ${j}.trimai -automated1 -resoverlap 0.55 -seqoverlap 60
		if [ -f ${j}.trimai ]; then
			seqtk seq -l0 ${j}.trimai > ${j}.trimal
			rm ${j}.trimai
		fi
done

