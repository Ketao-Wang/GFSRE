#!/usr/bin/perl -w
use strict;

# merge 5' and 3' flanking sequences of all gaps for alignment
my $cycle = $ARGV[0];
my $dir = $ARGV[1];		# gap directory
open(OUT,">$dir/../bwa_tmp/merge_cycle_$cycle.index.fa");		# temp folder for BWA

my @folders = glob("$dir/*");
foreach my $folder(@folders){
	my @tmp = split(/\//,$folder);
	$folder = $tmp[@tmp-1];
	if(!-e "$dir/$folder/cycle_$cycle/index.fa"){next;}
	open(IN,"$dir/$folder/cycle_$cycle/index.fa");
	while(<IN>){
		chomp;
		my ($tmp) = split(/\s+/,$_);
		if($tmp=~/^>/){
			$tmp=~s/\>//;
			print OUT ">$folder\_$tmp\n";
		}else{
			print OUT "$tmp\n";
		}
	}
	close IN;
}
close OUT;