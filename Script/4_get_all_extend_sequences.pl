#!/usr/bin/perl -w
use strict;

# Summarize the extending sequences of completely and partilly filled gap.
my @folders = glob("./Gap/*");
open(OUT1,">Completely_filled.sequence");
print OUT1 "Gap\tCompletely_filled_sequences\n";
open(OUT2,">Partilly_filled.sequence");
print OUT2 "Gap\t5'_Partially_filled_sequences\t3'_Partially_filled_sequences\n";
foreach my $folder(@folders){
	if(-e "$folder/gap_filled.fa"){
		open(IN,"$folder/gap_filled.fa");
		my @dir = split(/\//,$folder);
		my $gap = $dir[@dir-1];
		my $n = 0;
		while(<IN>){
			chomp;
			$n++;
			if($n == 2){
				print OUT1 "$gap\t$_\n";
			}
		}
		close IN;
	}elsif(-e "$folder/gap_shorten.fa"){
		open(IN,"$folder/gap_shorten.fa");
		my @dir = split(/\//,$folder);
		my $gap = $dir[@dir-1];
		my $n = 0;
		my $end5;
		my $end3;
		while(<IN>){
			chomp;
			$n++;
			if($n == 2){
				$end5 = $_;
			}elsif($n == 4){
				$end3 = $_;
			}
		}
		close IN;
		if($end5 || $end3){
			print OUT2 "$gap\t$end5\t$end3\n";
		}
	}
}
