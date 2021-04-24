#!/usr/bin/perl
system "ls *.faa > list.mafft";
open FILE, "<", "list.mafft";

while (<FILE>)
{
	chomp;
	if (/^(\S+)\.faa$/)
	{
			open SCRIPT, ">", "$1.mafft.in";
			print SCRIPT "#!/bin/bash\n\n";
			printf SCRIPT "einsi $1.faa > $1.mafft";
			close SCRIPT;
			system "chmod 755 $1.mafft.in";
			system "sbatch $1.mafft.in";
	}
}

close FILE;
system "rm list.mafft";
