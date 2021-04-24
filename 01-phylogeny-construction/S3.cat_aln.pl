#!/usr/bin/perl

#Before running this step, trimmed amino acid alignments should be produced 
#using the scripts (step1&2) in aa_seq folder
#===============================================
#  CONCATENAME MSA FROM DIFFERENT PROT FAMILIES
#  Usage: ./S3.cat_aln.pl PATH No.of.strains
#  e.g.  ./S3.cat_aln.pl ./single_copy_genes 12
# ===============================================

use strict;
use warnings;

if(@ARGV != 2)
{
	print "Pls specify a directory of gene alignments and number of strains\n";
}
my $dir = $ARGV[0]; #the directory including the trimmed amino acid alignments
system "ls $dir/*trimal > list.txt";
my $num = $ARGV[1]; #the number of total genomes

my %seq = ();
my $sid = "";


open GENPOS, ">partition_file.txt" or die $!;

open FILE, "<", "list.txt" or die $!;
while(my $input = <FILE>)
{
	chomp($input);
	my $gnm = `grep ">" $input | wc -l`;
	if($gnm == $num) ## the number of strains in your population
	{
		my $flg = 0;
		my ($fam) = $input =~ /(OG[0-9]+)\./;
		print GENPOS "AA, $fam = ";

		open IN, "$input" or die $!;
		while(my $line = <IN>)
		{
			chomp($line);
			if($line =~ /^>(\S+)\|\S+\|\S+/)
			{
				$sid = $1;
				#$sid =~ s/\|.*//g;
			}else
			{
				my $start = (exists $seq{$sid})? length($seq{$sid})+1:1;
				$seq{$sid} .= $line;
				my $end = length($seq{$sid});
				if($flg == 0)
				{
					print GENPOS "$start-$end\n";
					$flg = 1;
				}
			}
		}
		close IN;
	}

}
close FILE;

# Continue: partition_finder.cfg
#print GENPOS "\n## SCHEMES, search: all | user | greedy ##\n";
#print GENPOS "[schemes]\n";
#print GENPOS "search = rcluster;\n\n";
#print GENPOS "#user schemes go here if search=user. See manual for how to define.#\n";
close GENPOS;

# ==============================================================
# print out the concatinated alignment in fasta format,
# as well as phylip format
# =============================================================

my $lid = 0; # longest sid
my $snm = 0;
my $len = 0;

open CONCAT, ">concat.fas" or die $!;
foreach my $sid (sort{lc($a) cmp lc($b)} keys %seq)
{
	print CONCAT ">$sid\n",$seq{$sid},"\n";
	$lid = length($sid) if length($sid) > $lid;
	$snm++;
	$len = length($seq{$sid});
}
close CONCAT;

open PHYLIP, ">concat.phy" or die $!;
print PHYLIP "    $snm    $len\n";
foreach my $sid (sort{lc($a) cmp lc($b)} keys %seq)
{
	my $sps = " " x ($lid+5-length($sid));
	print PHYLIP $sid,$sps,$seq{$sid},"\n";
}
close PHYLIP;

system "rm list.txt";
exit 0;



