#!/usr/bin/env Rscript

# According to the supp. material CENSOR wu-blast by default uses parameters:
# blastn DB QU -gi hspmax=0 gspmax=0 B=100000000 V=0 gapE2=0.001 -warnings W 7 Q 30 R 6 S2 112 gapS2 225 S 225 X 225 gapX 450 -matrix 20p<CG>g
# where "In actual use, <CG> refers to G+C content of the query sequence. Censor incorporates substitution matrices that are optimized for sequences of different G+C content [3]. For examples matrix 20p39g is used for nucleotide sequences with average G+C content of 39%. The invocation of these matrices is done transparently to the user."


# According to running censor.ncbi with no input by default CENSOR.ncbi uses parameters:
# blastall -p blastn -b 10000000 -Y 1000000 -e 0.00002 -a 12 -r 3 -q -4 -y 80 -X 130 -Z 150 -W 10 -J T -m 8 -i QU -d DB -F 'm D'

# Here I have "translated" the legacy NCBI blast options into NCBI blast+ options and used the gap open and gap extend options from the wublast options (these two, -G and -E respectively, were not defined in the CENSOR documentation however are required to run).

# load packages
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))
suppressMessages(library(argparse))

parser <- ArgumentParser(description = "Rough reciprocal search for repeats. Sequences must have \".fasta\" extension")

parser$add_argument("-q", "--query", help = "File path to query sequence")
parser$add_argument("-g", "--genome", help = "File path to genome being searched")
parser$add_argument("-rd", "--repeat_database", help = "File path to database containing query repeat and all other repeats.")

args <- parser$parse_args()

# set variables
# query_path <- "~/Analysis/Snake/HT_Workflow/RTE-Snek/RTE-Snek.fasta"
# genome <- "~/Analysis/Genomes/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta"
# repeat_database <- "~/Analysis/Snake/HT_Workflow/Proto2-Snek/LINEs.fa"
query_path <- args$query
genome <- args$genome
repeat_database <- args$repeat_database

# create names of query and subject based on query path
query_name <- sub(".*\\/", "", query_path)
query_name <- sub("\\.fasta", "", query_name)
subject_name <- sub(".*\\/", "", genome)
subject_name <- sub("\\.fasta", "", subject_name)

# ensure blast db of genome exists
if(!file.exists(paste0(genome, ".nhr")) | !file.exists(paste0(genome, ".nin")) | !file.exists(paste0(genome, ".nsq"))){
  system(paste0("makeblastdb -dbtype nucl -in ", genome, " -out ", genome))
  }

# ensure working directory exists
if(!dir.exists("./working")){dir.create("./working")}

# run blast / read in blast table
system(paste0("blastn -dbsize 1000000 -max_target_seqs 1000000 -evalue 0.00002 -num_threads 4 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -outfmt 6 -dust yes -gapopen 30 -gapextend 6 -query ", query_path, " -db ", genome, " -out working/", query_name, "_in_", subject_name, ".out"))

blast_1 <- read_tsv(paste0("working/", query_name, "_in_", subject_name, ".out"), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# manipulate bl;ast data into useable format
blast_1 <- blast_1 %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "+", "-")) %>%
  dplyr::mutate(start = base::ifelse(sstart < send, sstart, send)) %>%
  dplyr::mutate(end = base::ifelse(sstart < send, send, sstart)) %>%
  dplyr::rename(seqnames = sseqid)

# filter out small hits
blast_1 <- blast_1 %>%
  filter(length >= 250)

# flatten ranges
blast_ranges <- blast_1 %>%
  plyranges::as_granges() %>%
  plyranges::reduce_ranges_directed() %>%
  GenomicRanges::reduce(min.gapwidth = 500L, ignore.strand=FALSE)

# convert blast ranges back to table for naming
blast_names <- blast_ranges %>%
  as_tibble() %>%
  dplyr::mutate(names = paste0(seqnames, ":", start, "-", end, "(", strand, ")"), names = as.character(names)) %>%
  select(seqnames, names, width)

# read in genome, rename scaffolds, empty garbage in case "leftovers" remain
genome_seq <- Biostrings::readDNAStringSet(filepath = genome)
names(genome_seq) <- sub(" .*", "", names(genome_seq))
gc(verbose = F)

# get sequences
blast_1_seqs <- Biostrings::getSeq(genome_seq, blast_ranges)

# names sequences
names(blast_1_seqs) <- blast_names$names

# write sequences to file
Biostrings::writeXStringSet(x = blast_1_seqs, filepath = paste0("working/", query_name, "_in_", subject_name, ".fasta"))

# reciprocal search - blast sequences against database which includes all LINEs in RepBase and LINEs of interest
system(paste0("blastn -dbsize 1000000 -max_target_seqs 1000000 -evalue 0.00002 -num_threads 4 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -outfmt 6 -dust yes -gapopen 30 -gapextend 6 -query working/", query_name, "_in_", subject_name, ".fasta", " -db ", repeat_database, " -out working/", query_name, "_in_", subject_name, "_2.out"))

# read in reciprocal search                      
blast_2 <- read_tsv(file = paste0("working/", query_name, "_in_", subject_name, "_2.out"), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# correct for strandness, rename where necessary
blast_2 <- blast_2 %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "+", "-")) %>%
  dplyr::mutate(start = base::ifelse(sstart < send, sstart, send)) %>%
  dplyr::mutate(end = base::ifelse(sstart < send, send, sstart)) %>%
  dplyr::rename(names = qseqid)

# select top hits based on coverage and bitscore
top_hits <- blast_2 %>%
  inner_join(blast_names) %>%
  mutate(coverage = (qend - qstart + 1) / width, pident = coverage * pident) %>%
  dplyr::group_by(names, sseqid) %>%
  summarise(sum_bitscore = sum(bitscore), sum_coverage = sum(coverage), sum_pident =  sum(pident)/sum(coverage)) %>% # calculate total coverage and sum bitscore for hits
  filter(sum_coverage > 0.5) %>% # coverage must be over 50%
  dplyr::top_n(n = 1, wt = sum_bitscore) %>% # select highest bitscore
  ungroup() %>%
  tidyr::extract("names", c("seqnames", "coords"), "(.*):([^_]+)$") %>%
  tidyr::separate(coords, into = c("coords", "strand"), sep = "\\(") %>%
  tidyr::separate(coords, into = c("start", "end"), sep = "-") %>%
  mutate(strand = sub("\\)", "", strand), start = as.double(start), end = as.double(end), width = end - start + 1) %>%
  arrange(seqnames, start)

# write to file
top_hits %>%
  select(seqnames, start, end, sseqid, sum_pident, strand) %>%
  readr::write_tsv(paste0("working/reciprocal_", query_name, "_in_", subject_name, ".bed"), col_names = F)