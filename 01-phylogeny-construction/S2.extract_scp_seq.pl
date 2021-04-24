#!/usr/bin/perl -w

`mkdir single_copy_gene`;
##1. get the map of gene to famid, and map of gene to seq
my @list=`ls ../00_genome_info/*faa`;
my @list=`ls ../00_genome_info/*gene`;
my %hash; #indexed with gene name, pointed to famid;
my %seq;
foreach my $file (@list)
{
	open IN,"$file" || die "can't open faa $!";
	my $gene;
	my ($famid) = $file =~ /_info\/(\S+)\.faa/;
	#my ($famid) = $file =~ /_info\/(\S+)\.gene/;
	while(<IN>)
	{
		chomp;
		if(s/>//)
		{
			$gene = (/^(\S+)\s/) ? $1 : $_;
			$hash{$gene} = $famid;
		}
		else
		{
			$seq{$gene} .= $_;
		}
	}
	close IN;
}

##2. 
my %SCO; 
open IN,"SingleCopyOrthogroups.txt" || die "can't open SCO $!";
while(<IN>)
{
	chomp;
	$SCO{$_}=1;
}
close IN;

`mkdir aa_seq nuc_seq`;
open IN, "Orthogroups.txt" || die "can't open Ortho $!";
while(<IN>)
{
	chomp;
	my ($famid) = $_ =~ /(^\S+):/;
	if( ! $SCO{$famid})
	{next;}
	
	open OUT, ">aa_seq/$famid.faa";
	open OUT, ">nuc_seq/$famid.dna";
	my @genes=split(" ");
	shift @genes;
	my %occur;
	foreach my $gene (@genes)
	{
		 print OUT ">".$gene."\n".$seq{$gene}."\n";
	}
	undef %occur;
	close OUT;
}
close IN;
