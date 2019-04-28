# GFSRE
## What is GFSRE ?
GFSRE (gap filling with single-read extension) performs high quality gap closure of the reference genome using single-end resequencing data. 

## Installation instructions
The package is written by perl scripts and no need to install. 
### Requirements
Burrows-Wheeler Aligner  
CAP3 Sequence Assembly Program  
BLAST+  

## Usage
### 0. Prepare reads for gap filling and the reference genome sequence.
#### Run with a simple example
>mkdir ./Genome  
>cp ./example/example_genome.fa ./Genome/genome.fa  
>mkdir ./Reads  
>cp ./example/example_reads.fq ./Reads/reads.fq  

Make sure the reads data is in './Reads/reads.fq', and the genome sequence is in './Genome/genome.fa'.

### 1. Generate gaps list in reference genome.
>perl ./Script/1_generate_gap_list.pl

### 2. Extract 100 bp flank sequences in both side of gaps as the index of extending
>perl ./Script/2_extract_flank_sequence_of_each_gap.pl  

The flank sequences are output into files './Gap/XXX/cycle_1/index.fa'.

### 3. Fill gaps with single-end reads.
>perl ./Script/3_repeatedly_fill_gap.pl

The reiterative strategy continues to close the gap sequences in step-by-step cycles. Considering the time requirement for the analysis, the maximum number of cycles was 30.  

### 4. Summarize the extending sequences of completely and partilly filled gap.
>perl ./Script/4_get_all_extend_sequences.pl  

Gap sequences for completely and partilly filled gaps was put in the 'Completely_filled.sequence' and 'Partilly_filled.sequence', respectively.
