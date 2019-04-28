#!/usr/bin/perl -w
use strict;

my $dir = "./Gap";		# gap directory
my $cycle_cutoff = 30;		# cycle number limit
system("mkdir ./bwa_tmp")

# split gaps into 20 group for parallelly filled
my @gap = `cat gap.list`;
my $m = int(@gap / 20);
open(IN,"gap.list");
open(OUT,">group.list");
my $n = 0;
while(<IN>){
	chomp;
	$n++;
	my $t = int(($n - 1) / $m) + 1;
	if($t > 20){$t = 20;}
	my @temp = split(/\t+|\s+/,$_);
	print OUT "$temp[0]_$temp[1]_$temp[2]\t$t\n";
}
close IN;
close OUT;

# de novo assembly repeatly
for(my $i = 1; $i <= $cycle_cutoff; $i++){	# 30 cycle
	system("perl ./Script/3_1_merge_index_seq.pl $i $dir");
	system("bwa index ./bwa_tmp/merge_cycle_$i.index.fa");
	system("bwa mem -M -t 24 ./bwa_tmp/merge_cycle_$i.index.fa ./Reads/reads.fq > ./bwa_tmp/index.sam 2>./bwa_tmp/bwa.log");
	for(my $j=1;$j<=20;$j++){	# 20 group
		system("nohup perl ./Script/3_2_split_into_each_gap_BYgroup.pl $i $dir $j &")
	}
	### judge if the assembly is complete according to the folder size
	my $size = 0;
	my $t = time();
	LABLE:while(2>1){
		my $t1 = time();
		if($t1 - $t > 300){		# more than 5min
			$t = $t1;
			my $size1 = `du --max-depth=1 | grep gap | cut -f 1`;	# gap folder size
			$size1 =~s/\r|\n//;
			if($size1 == $size){	# gap folder size not change
				print OUT "cycle_$i\t$t\t$size1\t\t\tNo change!\n";
				last LABLE;		# finish this cycle
			}else{		# gap folder size keeping increase
				print OUT "cycle_$i\t$t\t$size1\t\t\tChanging!\n";
				$size = $size1;	
			}
		}
	}
}

# merge sequences of all cycle
my @folders = glob("$dir/*");
foreach my $folder(@folders){
	chdir("$folder");
	system("perl ./Script/3_3_merge_all_contig.pl");
}
