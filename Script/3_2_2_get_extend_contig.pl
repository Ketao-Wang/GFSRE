#!/usr/bin/perl -w
use strict;

# assembled the mapped reads in 5’ and 3’respectively.
my @array;
system("sed -n 1,2p index.fa > 5end_index.fa");
system("sed -n 3,4p index.fa > 3end_index.fa");  
@array = `cat 5end_match_reads`;
if(@array>2000){		# only keep 1000 reads
	`cp 5end_match_reads 5end_match_reads2`;
	`head -2000 5end_match_reads2 > 5end_match_reads`;
}
system("cap3 5end_match_reads > 5end_match_reads.out");
@array = `cat 3end_match_reads`;
if(@array>2000){		# only keep 1000 reads
	`cp 3end_match_reads 3end_match_reads2`;
	`head -2000 3end_match_reads2 > 3end_match_reads`;
}
system("cap3 3end_match_reads > 3end_match_reads.out");
system("rm *ace *links *qual *info *singlets");

my @tags = qw/5end 3end/;
foreach my $tag(@tags){
	my $index = `sed -n 2p $tag\_index.fa`;
	system("blastn -query $tag\_index.fa -subject $tag\_match_reads.cap.contigs -outfmt 7 -evalue 0.1 -out $tag\_match_reads.cap.contigs.blast_out");		# align with index seq after firstly assembly
	my $contig_num;
	if(!-e "$tag\_match_reads.cap.contigs.blast_out"){
		$contig_num = 0;
	}else{
		$contig_num = 0;
		open(IN,"$tag\_match_reads.cap.contigs.blast_out");
		while(<IN>){
			chomp;
			if($_=~/^#/){next;}
			my @tmp = split(/\t+|\s+/,$_);
			if(($tmp[7]-$tmp[6]+1)> 0.95*length($index)){$contig_num++;}	# 95% similarity
		}
		close IN;
	}
	if($contig_num == 0){		# 1. no contig was same as index seq after firstly assembly
		next;
	}elsif($contig_num == 1){		# 2. only one contig was same as index seq after firstly assembly
		open(IN,"$tag\_match_reads.cap.contigs.blast_out");
		my $l = 0;
		while(<IN>){
			chomp;
			$l++;
			if($l == 6){
				my @tmp = split(/\t+|\s+/,$_);
				if(($tmp[7]-$tmp[6]+1) > 0.95*length($index)){		# 95% similarity
					open(IN2,"$tag\_match_reads.cap.contigs");
					open(OUT,">$tag\_final_contig.fa");	
					my $bool = 0;
					while(<IN2>){
						chomp;
						if($_=~/^>/){
							$_=~s/>//;
							if($_ eq $tmp[1]){
								$bool = 1;
								print OUT ">$tag\n";
								next;
							}else{
								$bool = 0;
							}
						}
						if($bool == 1){
							print OUT "$_";
						}
					}
					print OUT "\n";
					close IN2;
					close OUT;
					open(OUT,">$tag\_extend.fa");
					my $contig = `sed -n 2p $tag\_final_contig.fa`;
					$contig=~s/\r+|\n+//;
					my $index_seq;
					if($tag eq "5end"){
						$index_seq = substr($contig,$tmp[9]-50,length($contig)-$tmp[9]+50);
					}elsif($tag eq "3end"){
						$index_seq = substr($contig,0,$tmp[8]+50);
					}
					print OUT ">$tag\n$index_seq\n";
					close OUT;
				}
			}
		}
		close IN;
	}else{		# 3. more than one contig were same as index seq after firstly assembly, and then perform secondly assembly 
		@array = `cat $tag\_match_reads.cap.contigs`;
		if(@array>2000){
			`cp $tag\_match_reads.cap.contigs $tag\_match_reads.cap.contigs2`;
			`head -2000 $tag\_match_reads.cap.contigs2 > $tag\_match_reads.cap.contigs`;
		}
		system("cap3 $tag\_match_reads.cap.contigs > $tag\_match_reads.cap.contigs.out");
		system("rm *ace *links *qual *info *singlets");
		system("blastn -query $tag\_index.fa -subject $tag\_match_reads.cap.contigs.cap.contigs -outfmt 7 -evalue 0.1 -out $tag\_match_reads.cap.contigs.cap.contigs.blast_out");		# align with index seq after secondly assembly
		my $contig_num2;
		if(!-e "$tag\_match_reads.cap.contigs.cap.contigs.blast_out"){
			$contig_num2 = 0;
		}else{
			$contig_num2 = 0;
			open(IN,"$tag\_match_reads.cap.contigs.cap.contigs.blast_out");
			while(<IN>){
				chomp;
				if($_=~/^#/){next;}
				my @tmp = split(/\t+|\s+/,$_);
				if(($tmp[7]-$tmp[6]+1)> 0.95*length($index)){$contig_num2++;}
			}
			close IN;
		}
		if($contig_num2 == 0){		# 3.1. no contig was same as index seq after secondly assembly
			open(IN,"$tag\_match_reads.cap.contigs");
			my $t;
			my %seq;
			while(<IN>){
				chomp;
				if($_=~/^>/){
					$t = $_;
					$t=~s/\>//;
				}else{
					$seq{$t} .= $_;
				}
			}
			close IN;
			my %len;
			my %ctg;
			foreach my $k1(sort keys %seq){
				my $length = length($seq{$k1});
				$ctg{$seq{$k1}} = $k1;
				$len{$length}{$k1} = 1;	
			}
			### use the best match contig as extend sequence
			my $tmp = `sed -n 6p $tag\_match_reads.cap.contigs.blast_out`;
			my @tmps = split(/\t+|\s+/,$tmp);
			if($tmps[7]-$tmps[6]+1>80){
				open(OUT,">$tag\_final_contig.fa");
				print OUT ">$tag\n$seq{$tmps[1]}\n";
				close OUT;
				open(OUT,">$tag\_extend.fa");
				my $index_seq;
				if($tag eq "5end"){
					$index_seq = substr($seq{$tmps[1]},$tmps[9]-50,length($seq{$tmps[1]})-$tmps[9]+50);
				}elsif($tag eq "3end"){
					$index_seq = substr($seq{$tmps[1]},0,$tmps[8]+50);
				}
				print OUT ">$tag\n$index_seq\n";
				close OUT;
			}
		}elsif($contig_num2 == 1){		# 3.2. only one contig was same as index seq after secondly assembly
			open(IN,"$tag\_match_reads.cap.contigs.cap.contigs.blast_out");
			my $l = 0;
			while(<IN>){
				chomp;
				$l++;
				if($l == 6){
					my @tmp = split(/\t+|\s+/,$_);
					if(($tmp[7]-$tmp[6]+1) > 0.95*length($index)){		# 95% similarity
						open(IN2,"$tag\_match_reads.cap.contigs.cap.contigs");
						open(OUT,">$tag\_final_contig.fa");
						my $bool = 0;
						while(<IN2>){
							chomp;
							if($_=~/^>/){
								$_=~s/>//;
								if($_ eq $tmp[1]){
									$bool = 1;
									print OUT ">$tag\n";
									next;
								}else{
									$bool = 0;
								}
							}
							if($bool == 1){
								print OUT "$_";
							}
						}
						print OUT "\n";
						close IN2;
						close OUT;
						open(OUT,">$tag\_extend.fa");
						my $contig = `sed -n 2p $tag\_final_contig.fa`;
						$contig=~s/\r+|\n+//;
						my $index_seq;
						if($tag eq "5end"){
							$index_seq = substr($contig,$tmp[9]-50,length($contig)-$tmp[9]+50);
						}elsif($tag eq "3end"){
							$index_seq = substr($contig,0,$tmp[8]+50);
						}
						print OUT ">$tag\n$index_seq\n";
						close OUT;
					}
				}
			}
			close IN;
		}else{		# 3.3. more than one contig were same as index seq after secondly assembly
			open(IN,"$tag\_match_reads.cap.contigs.cap.contigs");
			my $t;
			my %seq;
			while(<IN>){
				chomp;
				if($_=~/^>/){
					$t = $_;
					$t=~s/\>//;
				}else{
					$seq{$t} .= $_;
				}
			}
			close IN;
			### use the best match contig as extend sequence
			my $tmp = `sed -n 6p $tag\_match_reads.cap.contigs.cap.contigs.blast_out`;
			my @tmps = split(/\t+|\s+/,$tmp);
			if($tmps[7]-$tmps[6]+1>80){
				open(OUT,">$tag\_final_contig.fa");
				print OUT ">$tag\n$seq{$tmps[1]}\n";
				close OUT;
				open(OUT,">$tag\_extend.fa");
				my $index_seq;
				if($tag eq "5end"){
					$index_seq = substr($seq{$tmps[1]},$tmps[9]-50,length($seq{$tmps[1]})-$tmps[9]+50);
				}elsif($tag eq "3end"){
					$index_seq = substr($seq{$tmps[1]},0,$tmps[8]+50);
				}
				print OUT ">$tag\n$index_seq\n";
				close OUT;
			}
		}
	}
}
