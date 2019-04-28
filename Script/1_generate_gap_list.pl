#!/usr/bin/perl -w
use strict;

# Generate gaps list in reference genome.
open(IN,"./Genome/genome.fa");
open(OUT,">gap.list");
my $chr;
my $gt = 0;
my $s = 0;
my $e = 0;
while(<IN>){
	chomp;
	if($_ =~/^>/){
		$_ =~s/\>//;
		$chr = $_;
		$gt = 0;
		$s = 0;
		$e = 0;
	}else{
		my $len=length($_);
		print "$chr\t$len\n";
		my @array;
		for(my $i=0;$i<length($_);$i++){
			$array[$i] = substr($_,$i,1);
		}
		my $bool = 0;
		for(my $j=0;$j<@array;$j++){
			if(($bool == 0)&&($array[$j] eq "N")){
				my $s = $j + 1;
				print OUT "$chr\t$s";
				$bool = 1;
			}elsif(($bool == 1)&&($array[$j] ne "N")){
				my $e = $j;
				print OUT "\t$e\tN\n";
				$bool = 0;
			}
		}
	}
}
close IN;
close OUT;