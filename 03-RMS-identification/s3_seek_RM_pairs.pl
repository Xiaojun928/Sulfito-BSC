#!/usr//bin/perl -w
#==========================================
#get the R-M systems based on the best hit 
#==========================================

my %RM;  #indexed with subunit name of R-Ms, pointed to the type
open IN,"/your/path/of/REBASE/enzymes.txt";
while(<IN>)
{
	chomp;
	my @b = split(/\t/);
	$RM{$b[0]}=$b[2];
}
close IN;

my @list = `ls *_diamond.txt`;
foreach my $file (@list)
{
	my %l2e;#indexed with locus, pointed to array of enzymes
	my %e2l;
	my %hash; #indexed with locus and hits, pointed to bitscore
	open IN,"$file";
	while(<IN>)
	{
		chomp;
		my @a=split(/\t/);
		my ($locus) = $a[0] =~ /\S+\|(\S+)\|\S+/;
		$hash{$locus}{$a[1]} = $a[11];
		push @{$l2e{$locus}},$a[1];
		push @{$e2l{$a[1]}},$locus;
	}
	close IN;
	
	my ($strain) = $file =~ /(\S+)_diamond.txt/;
	open OUT,">$strain\_RM_new.tsv";
	foreach my $gene (keys %l2e)
	{
		my $score = 0;
		if(scalar(@{$l2e{$gene}})==1)
		{
			my $e = @{$l2e{$gene}}[0];
			print OUT $gene."\t".$e."\t".$RM{$e}."\n";
		}
		else #more than one hits
		{
			my $flag = 0;
			my $bit = 0;
			my $score = 0;
			my $re;
			foreach my $e(@{$l2e{$gene}})
			{
				if($e=~/^M\.\S+/)
				{
					my ($RE) = $e =~ /M\.(\S+)/;
					if(exists $e2l{$RE})
					{
						$flag = 1;
						foreach my $l (@{$e2l{$RE}})
						{
							$score = $hash{$gene}{$e} + $hash{$l}{$RE};
							#print OUT $gene."\t".$e."\t".$hash{$gene}{$e}."\t".$l."\t".$RE."\t".$hash{$l}{$RE}."\t".$RM{$e}."\n";
							if($bit < $score)
							{$bit = $score;$re = $gene."\t".$e."\t".$hash{$gene}{$e}."\t".$l."\t".$RE."\t".$hash{$l}{$RE}."\t".$RM{$e}."\n";}
						}
					}
				}
				if($e!~/\./)
				{
					my $MT = "M.".$e;
					if(exists $e2l{$MT})
					{
						$flag = 1;
						my $bit = 0;
						my $score = 0;
						my $re;
						foreach my $l (@{$e2l{$MT}})
						{
							$score = $hash{$gene}{$e} + $hash{$l}{$MT};
							#print OUT $gene."\t".$e."\t".$hash{$gene}{$e}."\t".$l."\t".$MT."\t".$hash{$l}{$MT}."\t".$RM{$e}."\n";
							if($bit < $score)
							{$bit = $score;$re = $gene."\t".$e."\t".$hash{$gene}{$e}."\t".$l."\t".$MT."\t".$hash{$l}{$MT}."\t".$RM{$e}."\n";}
						}
					}
				}
				if($e=~/^S\.\S+/) #Type I R-M system has S,R,M subunits
				{
					my ($group) = $e =~ /S\.(\S+)/;
					my $MT = "M.".$group;
					my $RE = $group;
					if(exists $e2l{$MT} && exists $e2l{$RE})
					{
						$flag = 1;
						my $bit = 0;
						my $score = 0;
						my $re;
						foreach my $l1 (@{$e2l{$RE}})
						{
							foreach my $l2 (@{$e2l{$MT}})
							{
							  $score = $hash{$gene}{$e} + $hash{$l2}{$MT} + $hash{$l1}{$RE};
							  #print OUT $gene."\t".$e."\t".$hash{$gene}{$e}."\t".$l2."\t".$MT."\t".$hash{$l2}{$MT}."\t".$l1."\t".$RE."\t".$hash{$l1}{$RE}."\t".$RM{$e}."\n";
							if($bit < $score)
							{$bit = $score;$re = $gene."\t".$e."\t".$hash{$gene}{$e}."\t".$l2."\t".$MT."\t".$hash{$l2}{$MT}."\t".$l1."\t".$RE."\t".$hash{$l1}{$RE}."\t".$RM{$e}."\n";}
							}
						}
					}
				}
			}
			print OUT $re;
			if($flag == 0)
			{
				my $e = @{$l2e{$gene}}[0];
				print OUT $gene."\t".$e."\t".$RM{$e}."\n";
			}
		}
	}
	close OUT;
}
