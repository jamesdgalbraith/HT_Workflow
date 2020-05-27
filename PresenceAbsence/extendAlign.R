#!/usr/bin/env Rscript

## Script  for creating initial alignements for manual annotation of repeats

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

query_repeat_name <- "RTE-Snek_2"

species_hits <- read_tsv(paste0("../new_method/", query_repeat_name, "/species_info.tsv"), col_names = c("clade_name", "species_name", "genome_name"))

# set output folder
out_folder <- paste0("../new_method/", query_repeat_name, "/curation")
if(!dir.exists(out_folder)){dir.create(out_folder, recursive = T)}

for(i in 1:nrow(species_hits)){

  clade <- species_hits$clade_name[i]
  species_name <- species_hits$species_name[i]
  genome_name <- species_hits$genome_name[i]
  print(species_name)
  
  # set and read in genome and index
  genome_path <- paste0("~/Genomes/", clade, "/", species_name, "/", genome_name)
  
  query_repeat <- paste0("~/HT_Workflow/", query_repeat_name, "/", query_repeat_name, ".fasta")
  
  if(!file.exists(paste0(genome_path, ".nsq"))){system(paste0("makeblastdb -dbtype nucl -in ", genome_path))}
  
  print("blasting")
  
  blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 12 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 7 -dust yes -gapopen 30 -gapextend 6 -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length"))
  
  print("blasted")
  
  blast_out <- as_tibble(blast_out) %>% dplyr::arrange(-bitscore) %>% dplyr::filter(length > 1000)
  
  # create fasta index if necessary
  if(!file.exists(paste0(genome_path, ".fai"))){system(paste0("samtools faidx ", genome_path))}
  
  # read in index
  genome_fai <- read_tsv(paste0(genome_path, ".fai"), col_names = c("seqnames", "scaffold_length", "x3", "x4", "x5")) %>%
    dplyr::select(1:2)
  
  print("reading genome")
  
  # read in genome and garbage collect
  genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
  names(genome_seq) <- sub(" .*", "", names(genome_seq))
  gc()
  
  print("read")
  
  # set flanks, query and coverage
  flank5 <- 0
  flank3 <- 0
  
  # set and read in repeats
  repeat_seq <- Biostrings::readDNAStringSet(filepath = query_repeat)
  names(repeat_seq) <- sub(" .*", "", names(repeat_seq))
  repeat_name <- names(repeat_seq)
  
  # manipulate blast out
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
    dplyr::filter(length > 1000) %>%
    arrange(-length) %>%
    mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
           end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3),
           start = case_when(start <= 1 ~ 1, start > 1 ~ start),
           end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end)) %>%
    dplyr::select(seqnames, start, end, strand, bitscore)
  
  bed_tbl <- as_tibble(bed) %>%
    mutate(X4 = ".", X5 = ".") %>%
    dplyr::select(seqnames, start, end, X4, X5, strand)
  
  write_tsv(x = bed_tbl, path = paste0("../new_method/", query_repeat_name, "/curation/", species_name, ".bed"))
  
  # if less than two hits move on to next sequence
  if(nrow(as_tibble(bed_tbl))<2){
    next
  }

  if(nrow(bed_tbl)>30){
    bed <- bed %>%
      arrange(-bitscore) %>%
      dplyr::slice(1:30)
  }
  
  bed <- bed %>%
    mutate(name = paste0(seqnames, ":", start, "-", end, "(", strand, ")")) %>%
    as_granges()
  
  # get seqs
  seqs <- Biostrings::getSeq(genome_seq, bed)
  
  # name seqs
  names(seqs) <- bed$name
  
  # merge with transcribed repeat
  seqs <- c(repeat_seq, seqs)
  
  # write seqs to files
  Biostrings::writeXStringSet(x = seqs, filepath = paste0(out_folder, "/", repeat_name, "_", species_name, "_hits.fa"))
  
  if(!dir.exists(paste0(out_folder, "/aligned"))){dir.create(paste0(out_folder, "/aligned"))}

  # perform multiple alignment
  system(paste0("mafft --localpair --maxiterate 10 --thread 12 ", out_folder, "/", repeat_name, "_hits.fa > ", out_folder, "/aligned/", repeat_name, "_", species_name, "_aligned.fa"))
  # 
}
