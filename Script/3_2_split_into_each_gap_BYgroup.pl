#!/usr/bin/perl -w
use strict;

# according to the alignment result, the mapped reads of each gap were assembled for sequence extending.
my $group = $ARGV[2];
open(IN,"./group.list");
my %list;
while(<IN>){
	chomp;
	my @tmp = split(/\t+|\s+/,$_);
	if($tmp[1] eq $group){
		$list{$tmp[0]} = 1;
	}
}
close IN;

my $cycle = $ARGV[0];
my $dir = $ARGV[1];		# gap directory

my @folders = glob("$dir/*");
foreach my $folder(@folders){
	my @tmp = split(/\//,$folder);
	$folder = $tmp[@tmp-1];
	if(!exists $list{$folder}){next;}
	if(!-e "$dir/$folder/cycle_$cycle"){next;}
	chdir("$dir/$folder/cycle_$cycle");
	system("grep $folder $dir/../bwa_tmp/index.sam | awk '\$3!~/\*/' > index.match.sam");
	my @args = stat("index.match.sam");
	my $size = $args[7];					
	if($size > 100000000){		# filter too big Sam 
		system("rm index.match.sam");
		next;
	}
	### de novo assembly
	system("perl $dir/../Script/3_2_1_generate_for_cap3.pl");
	system("perl $dir/../Script/3_2_2_get_extend_contig.pl");			# get contig by assembly
	### obtain the extend sequence from 5' and 3' sides
	if((-e "5end_extend.fa")&&(-e "3end_extend.fa")){
		my $idx5 = `sed -n 2p 5end_index.fa`;
		my $idx3 = `sed -n 4p 3end_index.fa`;
		my $ext5 = `sed -n 2p 5end_extend.fa`;
		my $ext3 = `sed -n 2p 3end_extend.fa`;
		if(($idx5 eq $ext5)&&($idx3 eq $ext3)){		# extend sequence same as index sequence in two sides
			next;
		}elsif(($idx5 eq $ext5)&&($idx3 ne $ext3)){			# extend sequence same as index sequence in 5' side
			system("echo -e \">5end\n\"  > result.fa");
			system("cat 3end_extend.fa >> result.fa");
		}elsif(($idx5 ne $ext5)&&($idx3 eq $ext3)){			# extend sequence same as index sequence in 3' side
			system("cat 5end_extend.fa > result.fa");
			system("echo -e \">3end\n\"  >> result.fa");
		}else{
			system("cat 5end_extend.fa 3end_extend.fa > result.fa");
		}
	}elsif((!-e "5end_extend.fa")&&(!-e "3end_extend.fa")){		# no extend sequences
		next;
	}elsif((!-e "5end_extend.fa")&&(-e "3end_extend.fa")){		# only extend sequences in 3' side
		my $idx3 = `sed -n 4p 3end_index.fa`;
		my $ext3 = `sed -n 2p 3end_extend.fa`;
		if($idx3 eq $ext3){next;}		# extend sequence same as index sequence
		my $m = $cycle - 1;
		if($m==0){
			system("cp 5end_index.fa 5end_extend.fa");
		}else{
			system("cp $dir/$folder/cycle_$m/5end_extend.fa 5end_extend.fa");
		}
		system("echo -e \">5end\n\"  > result.fa");
		system("cat 3end_extend.fa >> result.fa");
	}elsif((-e "5end_extend.fa")&&(!-e "3end_extend.fa")){		# only extend sequences in 5' side
		my $idx5 = `sed -n 4p 5end_index.fa`;
		my $ext5 = `sed -n 2p 5end_extend.fa`;
		if($idx5 eq $ext5){next;}		# extend sequence same as index sequence
		my $m = $cycle - 1;
		if($m==0){
			system("cp 3end_index.fa 3end_extend.fa");
		}else{
			system("cp $dir/$folder/cycle_$m/3end_extend.fa 3end_extend.fa");
		}
		system("cat 5end_extend.fa > result.fa");
		system("echo -e \">3end\n\"  >> result.fa");
	}
	### judge the extend sequence from 5' and 3' sides whether exists overlap
	system("perl $dir/../Script/3_2_3_judge_extend_fasta_overlap.pl");
	if(!-e "Overlap.fa"){
		my $n = $cycle + 1;
		system("mkdir $dir/$folder/cycle_$n");
		system("cat result.fa > $dir/$folder/cycle_$n/index.fa");		# if not overlap, make extend sequences as index sequence for next cycle
	}
}