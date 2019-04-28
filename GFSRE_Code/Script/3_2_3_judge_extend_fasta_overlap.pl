#!/usr/bin/perl -w
use strict;

system("cat 5end_extend.fa 3end_extend.fa > judge.fa");
system("/home/ktwang/Code/cap3/CAP3/cap3 judge.fa > judge.fa.out");
system("rm *ace *links *qual *info *singlets");
my @args = stat("judge.fa.cap.contigs");
my $size = $args[7];
if($size > 0){
	open(IN,"judge.fa.cap.contigs");
	open(OUT,">Overlap.fa");
	print OUT ">assembly_sequence\n";
	my $bool = 0;
	while(<IN>){
		chomp;
		if($_=~/^>/){
			$bool++;
			next;
		}
		if($bool == 1){
			print OUT "$_";
		}
	}
	print OUT "\n";
	close IN;
	close OUT;
}
