#!/usr/bin/perl -w
#==========================================
#get the best hit for query genes 
#==========================================

my %RM; #indexed with subunit name of R-Ms, pointed to the details of this subunit
open IN,"/your/path/of/REBASE/enzymes.txt";  # Downloaded form REBASE http://rebase.neb.com/rebase/rebase.ftp.html
my $head=<IN>;
while(<IN>)
{
	chomp;
	my @b = split(/\t/);
	$RM{$b[0]}=$_;
}
close IN;

my @list = `ls *_diamond.txt`;
foreach my $file (@list)
{
	my %hash;
	open IN,"$file";
	while(<IN>)
	{
		chomp;
		my @a=split(/\t/);
		my ($locus) = $a[0] =~ /\S+\|(\S+)\|\S+/;
		if(!exists($hash{$locus}))
		{
			$hash{$locus} = $a[1];
		}
		else
		{next;}
	}
	close IN;
	
	my ($strain) = $file =~ /(\S+)_diamond.txt/;
	open OUT,">$strain\_RM.tsv";
	print OUT "locus\t".$head;
	foreach my $gene (keys %hash)
	{
		print OUT $gene."\t".$hash{$gene}."\t".$RM{$hash{$gene}}."\n";
	}
	close OUT;
}
