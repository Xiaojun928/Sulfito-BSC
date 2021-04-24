#!/usr/bin/perl
system "ls *.trimal  > list.mafft";
open FILE, "<", "list.mafft";

while (<FILE>)
{
	chomp;
	if (/^(\S+)\.trimal$/)
	{
			open SCRIPT, ">", "$1.iqtree.sh";
			print SCRIPT "#!/bin/bash\n\n";
			printf SCRIPT "iqtree -s $1.trimal -nt 1 -m MFP -mrate E,I,G,I+G -wbtl -mfreq FU,F -alrt 1000 -bb 1000 -pre $1\n";
			#printf SCRIPT "/home-user/software/iqtree/iqtree -s $1.trimal -nt 1 -m MFP -mrate E,I,G,I+G -nstop 500 -pers 0.5 -b 100 -pre $1\n";
			close SCRIPT;
			system "chmod 755 $1.iqtree.sh";
			system "sbatch $1.iqtree.sh";
	}
}

close FILE;
system "rm list.mafft";
