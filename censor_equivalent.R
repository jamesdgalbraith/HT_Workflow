#!/usr/bin/env Rscript

# According to the supp. material CENSOR wu-blast by default uses parameters:
# blastn DB QU -gi hspmax=0 gspmax=0 B=100000000 V=0 gapE2=0.001 -warnings W 7 Q 30 R 6 S2 112 gapS2 225 S 225 X 225 gapX 450 -matrix 20p<CG>g
# where "In actual use, <CG> refers to G+C content of the query sequence. Censor incorporates substitution matrices that are optimized for sequences of different G+C content [3]. For examples matrix 20p39g is used for nucleotide sequences with average G+C content of 39%. The invocation of these matrices is done transparently to the user."


# According to running censor.ncbi with no input by default CENSOR.ncbi uses parameters:
# blastall -p blastn -b 10000000 -Y 1000000 -e 0.00002 -a 12 -r 3 -q -4 -y 80 -X 130 -Z 150 -W 10 -J T -m 8 -i QU -d DB -F 'm D'

# Here I have "translated" the legacy NCBI blast options into NCBI blast+ options and used the gap open and gap extend options from the wublast options (these two, -G and -E respectively, were not defined in the CENSOR documentation however are required to run).


suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

query_path <- "~/Analysis/Snake/HT_Workflow/Proto2-Snek/Proto2-Snek.fasta"
genome <- "~/Analysis/Genomes/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta"
repeat_database <- "~/Analysis/Snake/HT_Workflow/Proto2-Snek/LINEs.fa"



# blast_1 <- read.table(text=system(paste0("blastn -dbsize 1000000 -max_target_seqs 1000000 -evalue 0.00002 -num_threads 4 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -outfmt 6 -dust yes -gapopen 30 -gapextend 6 -query ", query_path, " -db ", genome), intern = TRUE), col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# blast_1 <- read.table("~/Analysis/Snake/HT_Workflow/Proto2-Snek/tmp_blast.out", col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))


blast_1 <- blast_1 %>%
  as_tibble() %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "-", "+")) %>%
  dplyr::mutate(start = base::ifelse(sstart < send, sstart, send)) %>%
  dplyr::mutate(end = base::ifelse(sstart < send, send, sstart)) %>%
  dplyr::rename(seqnames = sseqid)

blast_1 <- blast_1 %>%
  filter(length >= 250)

blast_1 <- blast_1 %>%
  plyranges::as_granges() %>%
  plyranges::reduce_ranges_directed()

genome_seq <- Biostrings::readDNAStringSet(filepath = genome)

blast_1_seqs <- Biostrings::getSeq(genome_seq, blast_1)

paste0("blastn -dbsize 1000000 -evalue 0.00002 -num_threads 4 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -outfmt 6 -dust yes-gapopen 30 -gapextend 6  -query ", initial_hits, " -db ", repeat_database)