#!/usr/bin/perl -w
use strict;

# merge flank sequences of all gap for align
my $cycle = $ARGV[0];
my $dir = $ARGV[1];		# gap所在目录
open(OUT,">$dir/../bwa_tmp/merge_cycle_$cycle.index.fa");		# temp folder

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