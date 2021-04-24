#!/usr/bin/perl
system "ls *.trimal  > list.trimal";
open FILE, "<", "list.trimal";

while (<FILE>)
{
	chomp;
	if (/^(\S+)\.trimal$/)
	{
			open SCRIPT, ">", "$1.iqtree.sh";
			print SCRIPT "#!/bin/bash\n\n";
			printf SCRIPT "/your/path/of/iqtree/iqtree -s $1.trimal -nt 1 -m MFP -mrate E,I,G,I+G -nstop 500 -pers 0.5 -b 100 -pre $1\n";
			close SCRIPT;
			system "chmod 755 $1.iqtree.sh";
			system "sbatch $1.iqtree.sh";
	}
}

close FILE;
system "rm list.trimal";
