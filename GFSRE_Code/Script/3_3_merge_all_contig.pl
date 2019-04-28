#!/usr/bin/perl -w
use strict;
use File::Basename qw(&basename &dirname);
use Cwd;

my $rawfile=basename(cwd);
my $dirname=dirname(cwd);

my $overlap;
my %save;
my @folders = glob("cycle*");
for(my $i=1;$i<=@folders;$i++){
	chdir("$dirname/$rawfile/cycle_$i");
	my $tag;
	if(-e "5end_final_contig.fa"){
		open(IN,"5end_final_contig.fa");
		while(<IN>){
			chomp;
			my ($tmp) = split(/\t+|\s+/,$_);
			if($tmp=~m/\>/){
				$tmp =~s/\>//;
				$tag = $tmp;
			}else{
				$save{$tag}{$i} = $_;
			}
		}
		close IN;
	}
	if(-e "3end_final_contig.fa"){
		open(IN,"3end_final_contig.fa");
		while(<IN>){
			chomp;
			my ($tmp) = split(/\t+|\s+/,$_);
			if($tmp=~m/\>/){
				$tmp =~s/\>//;
				$tag = $tmp;
			}else{
				$save{$tag}{$i} = $_;
			}
		}
		close IN;
	}
	if(-e "Overlap.fa"){
		open(IN,"Overlap.fa");
		while(<IN>){
			chomp;
			if($_=~/^>/){next;}
			$overlap = $_;
		}
		close IN;
	}
}

chdir("$dirname/$rawfile");
open(OUT1,">5end_all_contig");
foreach my $k1(sort{$a <=> $b} keys %{$save{"5end"}}){
	print OUT1 ">5end_cycle$k1\n$save{\"5end\"}{$k1}\n";
}
close OUT1;
open(OUT2,">3end_all_contig");
foreach my $k1(sort{$a <=> $b} keys %{$save{"3end"}}){
	print OUT2 ">3end_cycle$k1\n$save{\"3end\"}{$k1}\n";
}
close OUT2;

if($overlap){
	system("cat 5end_all_contig 3end_all_contig > all_contig");
	open(OUT,">>all_contig");
	print OUT ">Overlap_contig\n$overlap\n";
	close OUT;
	system("/home/ktwang/Code/cap3/CAP3/cap3 all_contig > all_contig.out");
	system("blastn -query cycle_1/5end_index.fa -subject all_contig.cap.contigs -outfmt 7 -evalue 0.1 -out 5end.all_contig.blast_out");
	system("blastn -query cycle_1/3end_index.fa -subject all_contig.cap.contigs -outfmt 7 -evalue 0.1 -out 3end.all_contig.blast_out");
	my $a1 = `sed -n 6p 5end.all_contig.blast_out`;
	my $a2 = `sed -n 6p 3end.all_contig.blast_out`;
	my @s = split(/\t+|\s+/,$a1);
	my @e = split(/\t+|\s+/,$a2);
	unless((exists $s[1])&&(exists $e[1])){
		system("echo The end can not match with index sequence! >> Error.message");
		next;
	}
	if($s[1] ne $e[1]){next;}
	open(IN,"all_contig.cap.contigs");
	my $t;
	my %seq;
	while(<IN>){
		chomp;
		my ($tmp) = split(/\t+|\s+/,$_);
		if($tmp=~m/\>/){
			$tmp =~s/\>//;
			$t = $tmp;
		}else{
			$seq{$t} .= $_;
		}
	}
	close IN;
	my $fill = substr($seq{$s[1]},$s[9],$e[8]-$s[9]-1);
	open(OUT,">gap_filled.fa");
	print OUT ">gap_filled\n$fill\n";
	close OUT;
}else{
	my $seq5end;
	my $seq3end;
	my $n5 = `grep '>' 5end_all_contig | wc -l`;
	my $n3 = `grep '>' 3end_all_contig | wc -l`;
	if($n5 == 0){
		$seq5end = "";
	}else{
		if($n5 == 1){
			system("blastn -query cycle_1/5end_index.fa -subject 5end_all_contig -outfmt 7 -evalue 0.1 -out 5end.5end_all_contig.blast_out");
		}else{
			system("/home/ktwang/Code/cap3/CAP3/cap3 5end_all_contig > 5end_all_contig.out");
			system("blastn -query cycle_1/5end_index.fa -subject 5end_all_contig.cap.contigs -outfmt 7 -evalue 0.1 -out 5end.5end_all_contig.blast_out");
		}
		my $a = `sed -n 6p 5end.5end_all_contig.blast_out`;
		my @s = split(/\t+|\s+/,$a);
		if(!exists $s[1]){
			system("echo The end can not match with index sequence! >> Error.message");
		}
		if($n5 == 1){
			open(IN,"5end_all_contig");
		}else{
			open(IN,"5end_all_contig.cap.contigs");
		}
		my $t5;
		my %seq5;
		while(<IN>){
			chomp;
			my ($tmp) = split(/\t+|\s+/,$_);
			if($tmp=~m/\>/){
				$tmp =~s/\>//;
				$t5 = $tmp;
			}else{
				$seq5{$t5} .= $_;
			}
		}
		close IN;
		$seq5end = substr($seq5{$s[1]},$s[9],length($seq5{$s[1]})-$s[9]);
	}
	if($n3 == 0){
		$seq3end = "";
	}else{
		if($n3 == 1){
			system("blastn -query cycle_1/3end_index.fa -subject 3end_all_contig -outfmt 7 -evalue 0.1 -out 3end.3end_all_contig.blast_out");
		}else{
			system("/home/ktwang/Code/cap3/CAP3/cap3 3end_all_contig > 3end_all_contig.out");
			system("blastn -query cycle_1/3end_index.fa -subject 3end_all_contig.cap.contigs -outfmt 7 -evalue 0.1 -out 3end.3end_all_contig.blast_out");
		}
		my $a = `sed -n 6p 3end.3end_all_contig.blast_out`;
		my @s = split(/\t+|\s+/,$a);
		if(!exists $s[1]){
			system("echo The end can not match with index sequence! >> Error.message");
		}
		if($n3 == 1){
			open(IN,"3end_all_contig");
		}else{
			open(IN,"3end_all_contig.cap.contigs");
		}
		my $t3;
		my %seq3;
		while(<IN>){
			chomp;
			my ($tmp) = split(/\t+|\s+/,$_);
			if($tmp=~m/\>/){
				$tmp =~s/\>//;
				$t3 = $tmp;
			}else{
				$seq3{$t3} .= $_;
			}
		}
		close IN;
		$seq3end = substr($seq3{$s[1]},0,$s[8]-1);
	}	
	open(OUT,">gap_shorten.fa");
	print OUT ">gap_5end\n$seq5end\n>gap_3end\n$seq3end\n";
	close OUT;	
}
system("rm *ace *links *qual *info *singlets");