#!/bin/bash
#SBATCH -n 1

#	usage
echo "sh this.sh"
echo "require nucalign namely .nucalign.msa in ."

#	initialization

#	operation
for i in `ls *.nucalign`
do
		j=$(basename $i .nucalign)
		#time trimal -in $i -out ${j}.trimai -automated1 -resoverlap 0.55 -seqoverlap 60
		trimal -in $i -out ${j}.trimai -automated1 -resoverlap 0.55 -seqoverlap 60
		if [ -f ${j}.trimai ]; then
			seqtk seq -l0 ${j}.trimai > ${j}.trimal
			rm ${j}.trimai
		fi
done

