#!/usr/bin/perl -w
use strict;

# extract 100 bp flanking sequences in both side of gaps, and output into './Gap/XXX/cycle_1/index.fa'.
open(IN1,"./Genome/genome.fa");
my %genome;
my $chr;
while(<IN1>){
	chomp;
	if($_=~/^\>/){
		$chr = $_;
		$chr=~s/\>//;
	}else{
		$genome{$chr} = $_;
	}

}
close IN1;

system("mkdir ./Gap");
open(IN,"gap.list");
while(<IN>){
	chomp;
	my @temp = split(/\t+|\s+/,$_);
	system("mkdir ./Gap/$temp[0]_$temp[1]_$temp[2]");
	system("mkdir ./Gap/$temp[0]_$temp[1]_$temp[2]/cycle_1");
	open(OUT,">./Gap/$temp[0]_$temp[1]_$temp[2]/cycle_1/index.fa");
	my $seq5 = substr($genome{$temp[0]},$temp[1]-100-1,100);
	my $seq3 = substr($genome{$temp[0]},$temp[2],100);
	print OUT ">5end $temp[0]_$temp[1]_$temp[2]\n$seq5\n";
	print OUT ">3end $temp[0]_$temp[1]_$temp[2]\n$seq3\n";
	close OUT;
}
close IN;
