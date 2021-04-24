#!/bin/bash
#SBATCH -n 16

#	usage
echo "sh this.sh"

#	initialization


time iqtree -nt 16 -m MFP -bb 1000 -alrt 1000 -s concat.phy -redo -mrate E,I,G,I+G -mfreq FU -wbtl -pre iqtree -spp partition_file.txt

#time iqtree -nt 14  -bb 1000 -s concat.phy -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre iqtree -m MFP+MERGE -rcluster 10 -bnni -spp example.partitions.txt
#example.partitions.txt:
#DNA, part1 = 1-100
#DNA, part2 = 101-384

#time iqtree -nt 16 -m MFP -bb 1000 -s concat.phy -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre iqtree -o GNM969
#for 3 genomes
#iqtree -m MFP -s concat.phy -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre iqtree

