#!/usr/bin/env Rscript

## Script  for creating alignements for manual annotation of repeats
## TO DO: CREATE VERSION FOR MULTIPLE SEQUENCES

library(argparse)

parser <- ArgumentParser(description = "Script to extract flanking sequences around repeats")

parser$add_argument("-f5", "--flank5", default = 2000, type = "double", help = "Length of 5' flanking sequence")
parser$add_argument("-f3", "--flank3", default = 2000, type = "double", help = "Length of 3' flanking sequence")
parser$add_argument("-qr", "--query_repeat", default = NULL, help = "Path to query repeat")
parser$add_argument("-gn", "--genome", default = NULL, help = "Genome to search for repeat")
parser$add_argument("-cv", "--coverage", default = 80, type = "double", help = "Minimum coverage by query repeat")
parser$add_argument("-pi", "--pident", default = 95, type = "double", help = "Minimum percent identity to consensus")
parser$add_argument("-of", "--out_folder", default = ".", help = "Folder for output")

args <- parser$parse_args()

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set and read in genome
if(!file.exists(paste0(args$genome, ".fai"))){system(paste0("samtools faidx ", args$genome))}

genome_fai <- read_tsv(paste0(args$genome, ".fai"), col_names = c("seqnames", "scaffold_length", "x3", "x4", "x5")) %>%
  dplyr::select(1:2)

# set and read in repeats
repeat_seq <- Biostrings::readDNAStringSet(filepath = args$query_repeat)
repeat_name <- names(repeat_seq)

print("Blasting")

# blast, read blast, rearrange
blast_out <- read.table(text=system(paste0("blastn -query ", args$query_repeat, " -db ", args$genome, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length")) %>%
  as_tibble() %>%
  dplyr::mutate(seqnames = as.character(seqnames), sstart = as.double(sstart), send = as.double(send), length = as.double(length)) %>%
  dplyr::arrange(-bitscore) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send))

print("Filtering")

# filter based on coverage and identity, extend and adjust
bed <- blast_out %>%
  inner_join(genome_fai) %>%
  dplyr::filter(length * 100 / width(repeat_seq) > args$coverage) %>%
  dplyr::filter(pident >= args$pident) %>%
  arrange(-length) %>%
  mutate(start = case_when(strand == "+" ~ start - args$flank5, strand == "-" ~ start - args$flank3),
         end = case_when(strand == "-" ~ end + args$flank5, strand == "+" ~ end + args$flank3),
         start = case_when(start <= 1 ~ 1, start > 1 ~ start),
         end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end),
         name = paste0(seqnames, ":", start, "-", end, "(", strand, ")")) %>%
  dplyr::select(seqnames, start, end, strand, name)

if(nrow(bed) < 2){stop("Too few hits, try adjusting pident and coverage")}

print("Snipping")
print(paste0("hits = ", nrow(bed)))
if(nrow(as_tibble(bed))>30){
  bed <- bed %>% dplyr::slice(1:30)
}

print(head(bed))

print("Ranging")

bed <- plyranges::as_granges(bed)

head(bed)

print("Reading")

# read in genome
genome_seq <- Biostrings::readDNAStringSet(filepath = args$genome)

print("Renaming")

# rename  genome seqs
names(genome_seq) <- sub(" .*", "", names(genome_seq))

print(head(names(genome_seq)))
print(head(bed))

print("Getting seq")

# get seqs
seqs <- Biostrings::getSeq(genome_seq, bed)

print("Naming seq")

# name seqs
names(seqs) <- GenomicRanges::elementMetadata(bed)[["name"]]

print("Merging with query")

# merge with transcribed repeat
seqs <- c(repeat_seq, seqs)

print("Writing")

# write seqs to files
Biostrings::writeXStringSet(x = seqs, filepath = paste0(args$out_folder, "/", repeat_name, "_hits.fa"))

print("Aligning")

# perform multiple alignment
system(paste0("mafft --adjustdirection ", args$out_folder, "/", repeat_name, "_hits.fa > ", args$out_folder, "/", repeat_name, "_hits_aligned.fa"))

