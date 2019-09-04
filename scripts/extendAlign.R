#!/usr/bin/env Rscript

## Script  for creating alignements for manual annotation of repeats
## TO DO: CREATE VERSION FOR MULTIPLE SEQUENCES

library(argparse)
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))
parser <- ArgumentParser(description = "Script to extract flanking sequences around repeats")

parser$add_argument("-f5", "--flank5", default = 2000, type = "double", help = "Length of 5' flanking sequence")
parser$add_argument("-f3", "--flank3", default = 2000, type = "double", help = "Length of 3' flanking sequence")
parser$add_argument("-qr", "--query_repeat", default = NULL, help = "Path to query repeats")
parser$add_argument("-gn", "--genome", default = NULL, help = "Genome to search for repeat")
parser$add_argument("-cv", "--coverage", default = 80, type = "double", help = "Minimum coverage by query repeat")
parser$add_argument("-pi", "--pident", default = 95, type = "double", help = "Minimum percent identity to consensus")
parser$add_argument("-of", "--out_folder", default = "output/", help = "Folder for output")

args <- parser$parse_args()

# set and read in genome
if(!file.exists(paste0(args$genome, ".fai"))){system(paste0("samtools faidx ", args$genome))}
if(!dir.exists(args$out_folder)){dir.create(args$out_folder)}

genome_fai <- read_tsv(paste0(args$genome, ".fai"), col_names = c("seqnames", "scaffold_length", "x3", "x4", "x5")) %>%
  dplyr::select(1:2)
genome_seq <- Biostrings::readDNAStringSet(filepath = args$genome)

# set and read in repeats
repeat_seq <- Biostrings::readDNAStringSet(filepath = args$query_repeat)
repeat_names <- names(repeat_seq)

# run script over all repeats
purrr::map(.x = repeat_names, ~{
  # Select repeat
  working_repeat <- repeat_seq[.]
  
  # Get repeat name
  working_name <- names(working_repeat)
  
  # Write repeat to file
  Biostrings::writeXStringSet(x = working_repeat, filepath = paste0(args$out_folder, "/", working_name, ".fa"))
  
  # Search for repeat, rearrange output
  blast_out <- read.table(text=system(paste0("blastn -query ", args$out_folder, "/", working_name, ".fa -db ", args$genome, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length")) %>%
    as_tibble() %>%
    dplyr::mutate(seqnames = as.character(seqnames), sstart = as.double(sstart), send = as.double(send), length = as.double(length)) %>%
    dplyr::arrange(-bitscore) %>%
    mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
           start = case_when(sstart < send ~ sstart, send < sstart ~ send),
           end = case_when(sstart > send ~ sstart, send > sstart ~ send))
  
  # Remove repeat file
  file.remove(paste0(args$out_folder, "/", working_name, ".fa"))
  
  # Filter search table based on coverage and identity, extend and adjust
  bed <- blast_out %>%
    inner_join(genome_fai) %>%
    dplyr::filter(length * 100 / width(working_repeat) > args$coverage) %>%
    dplyr::filter(pident >= args$pident) %>%
    arrange(-length) %>%
    mutate(start = case_when(strand == "+" ~ start - args$flank5, strand == "-" ~ start - args$flank3),
           end = case_when(strand == "-" ~ end + args$flank5, strand == "+" ~ end + args$flank3),
           start = case_when(start <= 1 ~ 1, start > 1 ~ start),
           end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end),
           name = paste0(seqnames, ":", start, "-", end, "(", strand, ")")) %>%
    dplyr::select(seqnames, start, end, strand, name)

  # If too big, reduce size of hits table
  if(nrow(as_tibble(bed))>30){
    bed <- bed %>% dplyr::slice(1:30)
  }

  # Convert hit table to bed
  bed <- plyranges::as_granges(bed)

  # Get seqs
  seqs <- Biostrings::getSeq(genome_seq, bed)

  # Name seqs
  names(seqs) <- GenomicRanges::elementMetadata(bed)[["name"]]
  
  # Combine seqs with consensus repeat
  seqs <- c(working_repeat, seqs)

  # Write seqs to files
  Biostrings::writeXStringSet(x = seqs, filepath = paste0(args$out_folder, "/", working_name, "_hits.fa"))

  # Perform multiple alignment
  system(paste0("mafft --adjustdirection ", args$out_folder, "/", working_name, "_hits.fa > ", args$out_folder, "/", working_name, "_hits_aligned.fa"))
  
  # Remove unaligned multifasta
  file.remove(paste0(args$out_folder, "/", working_name, "_hits.fa"))
  
})