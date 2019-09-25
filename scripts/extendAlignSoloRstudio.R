#!/usr/bin/env Rscript

## Script  for creating alignements for manual annotation of repeats
## TO DO: CREATE VERSION FOR MULTIPLE SEQUENCES

# library(argparse)

# parser <- ArgumentParser(description = "Script to extract flanking sequences around repeats")
# 
# parser$add_argument("-f5", "--flank5", default = 2000, type = "double", help = "Length of 5' flanking sequence")
# parser$add_argument("-f3", "--flank3", default = 2000, type = "double", help = "Length of 3' flanking sequence")
# parser$add_argument("-qr", "--query_repeat", default = NULL, help = "Path to query repeat")
# parser$add_argument("-gn", "--genome", default = NULL, help = "Genome to search for repeat")
# parser$add_argument("-cv", "--coverage", default = 80, type = "double", help = "Minimum coverage by query repeat")
# parser$add_argument("-pi", "--pident", default = 95, type = "double", help = "Minimum percent identity to consensus")
# parser$add_argument("-of", "--out_folder", default = ".", help = "Folder for output")
# 
# args <- parser$parse_args()

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set output folder
out_folder <- "~/Analysis/Snake/HT_Workflow/RTE-Snek/curation/"
if(!dir.exists(out_folder)){dir.create(out_folder, recursive = T)}

# set and read in genome and index
genome_path <- "~/Analysis/Genomes/Lethenteron_camtschaticum/LetJap1.0.fasta"

genome_name <- sub(paste0(".*\\/"), "", genome_path)
genome_name <- sub(".fasta", "", genome_name)

if(!file.exists(paste0(genome_path, ".fai"))){system(paste0("samtools faidx ", genome_path))}

genome_fai <- read_tsv(paste0(genome_path, ".fai"), col_names = c("seqnames", "scaffold_length", "x3", "x4", "x5")) %>%
  dplyr::select(1:2)

genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()

names(genome_seq) <- sub(" .*", "", names(genome_seq))

# set flanks, query and coverage
flank5 <- 1000
flank3 <- 100
query_repeat <- "~/Analysis/Snake/HT_Workflow/RTE-Snek/curation/RTE-letJap_3.fasta"
coverage <- 60
pident <- 90

# set and read in repeats
repeat_seq <- Biostrings::readDNAStringSet(filepath = query_repeat)
names(repeat_seq) <- sub(" .*", "", names(repeat_seq))
repeat_name <- names(repeat_seq)

# blast, read blast, rearrange
# For sequence from other species
# blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 4 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -dust yes -gapopen 30 -gapextend 6 -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length"))

# For sequence from same species
blast_out <- read.table(text=system(paste0("blastn -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length"))

blast_out <- blast_out %>%
  as_tibble() %>%
  dplyr::mutate(seqnames = as.character(seqnames), sstart = as.double(sstart), send = as.double(send), length = as.double(length)) %>%
  dplyr::arrange(-bitscore) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send))

# filter based on coverage and identity, extend and adjust
bed <- blast_out %>%
  inner_join(genome_fai) %>%
  # dplyr::filter(pident >= 90, length > 1000) %>%
  dplyr::filter(length >= 1000) %>%
  arrange(-bitscore) %>%
  mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
         end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3),
         start = case_when(start <= 1 ~ 1, start > 1 ~ start),
         end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end),
         name = paste0(seqnames, ":", start, "-", end, "(", strand, ")")) %>%
  dplyr::select(seqnames, start, end, strand, name)

if(nrow(as_tibble(bed))>20){
  bed <- bed %>% 
    dplyr::slice(1:20)
}

bed <- plyranges::as_granges(bed)

# get seqs
seqs <- Biostrings::getSeq(genome_seq, bed)

# name seqs
names(seqs) <- GenomicRanges::elementMetadata(bed)[["name"]]

# merge with transcribed repeat
seqs <- c(repeat_seq, seqs)

# write seqs to files
Biostrings::writeXStringSet(x = seqs, filepath = paste0(out_folder, "/", repeat_name, "_hits.fa"))

# perform multiple alignment
system(paste0("mafft --localpair --thread 4 ", out_folder, "/", repeat_name, "_hits.fa > ", out_folder, "/", repeat_name, "_hits_in_", genome_name, "_aligned.fa"))

