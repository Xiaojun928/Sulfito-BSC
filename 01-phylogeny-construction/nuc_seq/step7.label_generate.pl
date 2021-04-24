#!/usr/bin/perl -w

my %gen;
open IN,"../gnm_clade_list.txt";
while(<IN>)
{
 chomp;
 my @a = split(/\t/,$_);
 $gen{$a[0]} = $a[1];
 $hbt{$a[0]} = $a[2];
}
close IN;


system "ls *.label > list.txt";
open IN,"list.txt" or die "can't open0 $!";
while (my $file=<IN>)
{
  chomp $file;
  my ($name) = $file =~ /(\S+).label/;
  open OUT,">$name.labels";
  open IN1,"$file" || die "can't open1 $!";
  while(<IN1>)
  {
     chomp;
     my @a = split(/\_/);
     my $gnm = join ("\_",@a[0..2]);
     print OUT $_."\t".$gen{$gnm}."\t".$hbt{$gnm}."\n";
  }
  close IN1;
  
  close OUT;
}
close IN;
system "rm list.txt";

