#!/usr/bin/perl -w
use strict;

open(IN,"index.match.sam");
open(OUT1,">5end_match_reads");
open(OUT2,">3end_match_reads");
while(<IN>){
	chomp;
	if($_=~/^@/){next;}
	my @tmp = split(/\t+|\s+/,$_); 
	if($tmp[2]=~m/5end/){
		print OUT1 ">$tmp[0]\n$tmp[9]\n";
	}elsif($tmp[2]=~m/3end/){
		print OUT2 ">$tmp[0]\n$tmp[9]\n";
	}
}
close IN;
close OUT1;
close OUT2;