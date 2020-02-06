#!/usr/bin/env Rscript

## Script  for creating initial alignements for manual annotation of repeats

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

species_hits <- read_tsv("~/all_genomes.txt", col_names = c("clade_name", "species_name", "genome_name"))
query_repeat <- "~/HT_Workflow/new_method/Rex1-Snek_1/Rex1-Snek_1.fasta"

# set output folder
out_folder <- "~/HT_Workflow/new_method/Rex1-Snek_1/"

if(!dir.exists(out_folder)){dir.create(out_folder, recursive = T)}
if(!dir.exists(paste0(out_folder, "/hits"))){dir.create(paste0(out_folder, "/hits"), recursive = T)}
if(!dir.exists(paste0(out_folder, "/aligned"))){dir.create(paste0(out_folder, "/aligned"), recursive = T)}
if(!dir.exists(paste0(out_folder, "/best_hits"))){dir.create(paste0(out_folder, "/best_hits"), recursive = T)}

for(i in 1:nrow(species_hits)){
  clade <- species_hits$clade_name[i]
  species_name <- species_hits$species_name[i]
  genome_name <- species_hits$genome_name[i]
  print(species_name)
  
  # set and read in genome and index
  genome_path <- paste0("~/Genomes/", clade, "/", species_name, "/", genome_name)
  
  if(!file.exists(paste0(genome_path, ".nsq"))){system(paste0("makeblastdb -dbtype nucl -in ", genome_path))}
  
  blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 12 -word_size 7 -dust yes -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length")) %>% as_tibble()
  
  # if no reuslts skip to next species
  if(nrow(blast_out) < 1){
    next
  }
  
  # set and read in repeats
  repeat_seq <- Biostrings::readDNAStringSet(filepath = query_repeat)
  names(repeat_seq) <- sub(" .*", "", names(repeat_seq))
  repeat_name <- names(repeat_seq)
  
  # create fasta index if necessary
  if(!file.exists(paste0(genome_path, ".fai"))){system(paste0("samtools faidx ", genome_path))}
  
  # read in index
  genome_fai <- read_tsv(paste0(genome_path, ".fai"), col_names = c("seqnames", "scaffold_length", "x3", "x4", "x5")) %>%
    dplyr::select(1:2)
  
  # read in genome and garbage collect
  genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
  gc()
  
  # name genome sequences to fit blast output
  names(genome_seq) <- sub(" .*", "", names(genome_seq))
  
  # extract best hit > 1000bp
  initial_bed <- blast_out %>%
    dplyr::arrange(-pident) %>% dplyr::filter(length > 1000) %>% dplyr::slice(1) %>%
    mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
           start = case_when(sstart < send ~ sstart, send < sstart ~ send),
           end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
    select(seqnames, start, end, strand) %>%
    plyranges::as_granges()
  
  initial_seq <- Biostrings::getSeq(genome_seq, initial_bed)
  
  # name seqs
  names(initial_seq) <- paste0(seqnames(initial_bed), ":", ranges(initial_bed), "(", strand(initial_bed), ")")
  
  # write seqs to files
  Biostrings::writeXStringSet(x = initial_seq, filepath = paste0(out_folder, "/best_hits/", repeat_name, "_best_hit_in_", species_name, ".fa"))
  
  # set flanks
  flank5 <- 2000
  flank3 <- 2000
  
  # search for best hit
  blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 12 -word_size 7 -dust yes -query ", out_folder, "/best_hits/", repeat_name, "_best_hit_in_", species_name, ".fa -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length slen\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "slen")) %>% as_tibble()
  
  second_bed <- blast_out %>%
    dplyr::filter(pident > 94, length > 800) %>%
    dplyr::arrange(-pident) %>%
    dplyr::slice(1:20) %>%
    mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
           start = case_when(sstart < send ~ sstart, send < sstart ~ send),
           end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
    mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
           end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3),
           start = case_when(start <= 1 ~ 1, start > 1 ~ start),
           end = case_when(end > as.double(slen) ~ as.double(slen), end <= as.double(slen) ~ end)
           ) %>%
    select(seqnames, start, end, strand) %>%
    plyranges::as_granges()
  
  # get seqs
  second_seqs <- Biostrings::getSeq(genome_seq, second_bed)
  
  # name seqs
  names(second_seqs) <- paste0(seqnames(second_bed), ":", ranges(second_bed), "(", strand(second_bed), ")")
  
  # write seqs to files
  Biostrings::writeXStringSet(x = second_seqs, filepath = paste0(out_folder, "/hits/", species_name, "_hits.fa"))
  
  # perform multiple alignment
  system(paste0("mafft --localpair --maxiterate 10 --thread 12 ", out_folder, "/hits/", species_name, "_hits.fa > ", out_folder, "/aligned/", species_name, "_aligned.fa"))
  
}
