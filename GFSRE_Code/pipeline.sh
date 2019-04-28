# 0.Collect reads for gap filling (./Reads/reads.fq) and prepare the genome sequence (./Genome/genome.fa).

# 1.generate gaps list in reference genome (delete gaps in chr.0 ).
perl 1_generate_gap_list.pl
# 2.extract 100 bp flank sequences in both side of gaps, and output into 'cycle_1' folder of each gap.
mkdir ./Gap
perl 2_extract_flank_sequence_of_each_gap.pl
# 3.fill gap
perl 3_repeatedly_fill_gap.pl
