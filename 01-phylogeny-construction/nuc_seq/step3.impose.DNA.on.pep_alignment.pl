#!/usr/bin/perl

system "ls *.dna > list.txt";
open FILE, "<", "list.txt";
while(<FILE>)
{
	chomp;
	if(/^(\S+)\.dna$/)
	{
		if(-e "$1.mafft")
		{
			impose_DNA($1);
		}
		else
		{
			print "$1.pepalign does not exist!\n";
		}
	}
}
close FILE;
system "rm list.txt";


sub impose_DNA{
        my($gene)=@_;
        open SUBREF, "<", "$gene.dna";
        open SUBCOMP, "<", "../aa_seq/$gene.mafft";
        open SUBOUT, ">", "$gene.nucalign";
        my %hash=();
        while(<SUBREF>)
        {
                chomp;
                if(/^>(\S+)$/)
                {
			my $id=$1;
                        my $seq=<SUBREF>;
                        chomp $seq;
                        $hash{$id}=$seq;
                }
        }
        while(<SUBCOMP>)
        {
                chomp;
                /^\s*$/ and next;
                if(/^>(\S+)$/)
                {
			my $id=$1;
                        print SUBOUT ">$1\n";
                        my $pep=<SUBCOMP>;
                        chomp $pep;
                        my $len=length $pep;
                        my $nuc=$hash{$id};
                        my $c=0;
                        for(my $x=0; $x<$len; $x++)
                        {
                                my $aa=substr $pep, $x, 1;
                                if($aa eq '-')
                                {
                                        print SUBOUT "---";
                                }
                                else
                                {
                                        my $codon=substr $nuc, $c, 3;
                                        if($codon eq '')
                                        {
                                                print "In family $gene, the gene $id has different DNA and pep lengths\n";
                                        }
                                        print SUBOUT "$codon";
                                        $c+=3;
                                }
                        }
                        print SUBOUT "\n";
                }
        }
        close SUBCOMP;
        close SUBREF;
        close SUBOUT;
}

